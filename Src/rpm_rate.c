/*
 * rpm_rate.c
 *
 * Time derivative of the electrical RPM (eRPM rate, i.e. rotor acceleration), estimated from the
 * electrical-revolution period e_com_time [us] for the extended KISS telemetry frame.
 *
 * Measurement
 *   e_com_time T is the duration of the last electrical revolution (sum of the last six commutation intervals,
 *   see main.c), so every sample yields the speed averaged over that revolution:
 *     z = 60 * 10^6 / T  [eRPM]
 *   Samples arrive every dt = 1.05 ms from the 1 kHz block of tenKhzRoutine().
 *
 * Savitzky-Golay derivative: quadratic fit over the last N = 2m + 1 = 17 samples, differentiated at the centre
 *   1. Fit samples: The samples get the index j \in [-m,...,m] counted from the centre sample: j = 0 is m samples old,
 *      j = m is the newest. A quadratic is fitted to them by least squares:
 *        p(j) = c0 + c1*j + c2*j^2,   minimising  sum_j (z_j - p(j))^2
 *      The window is symmetric about j = 0, so the normal equations decouple. The odd coefficient depends only on
 *      the odd moment of the data:
 *        c1 = sum_j j*z_j / sum_j j^2
 *      while c0 and c2 follow from sum_j z_j and sum_j j^2*z_j alone.
 *   2. Differentiate at the centre. dp/dj = c1 + 2*c2*j, so the slope at j = 0 is c1 [eRPM per sample]. c0 and c2
 *      do not enter it, which is why they are not computed. Per second:
 *        rate = c1 / dt  [eRPM/s]
 *   Together this is a fixed weighted sum of the samples, weights j / (dt * sum j^2): a linear-phase FIR filter.
 *   Its delay is exactly m samples (8 * 1.05 = 8.4 ms) at all frequencies, plus the averaging over one electrical
 *   revolution in e_com_time.
 *
 * Integer arithmetic only. Called from the main loop with the e_com_time samples taken at the fixed rate of the
 * 1 kHz block in tenKhzRoutine(), which keeps the divisions out of the interrupt.
 */

#include "rpm_rate.h"
#include "targets.h"

// The samples come from the 1 kHz block of tenKhzRoutine(), which runs every (PID_LOOP_DIVIDER + 1) interrupts,
// so dt = RPM_RATE_TICKS_PER_UPDATE / LOOP_FREQUENCY_HZ = 21 / 20000 s = 1.05 ms at 20 kHz
#define RPM_RATE_TICKS_PER_UPDATE (PID_LOOP_DIVIDER + 1)
// Update frequency of 1 / dt in Hz, rounded (952 at 20 kHz, -0.04 %)
#define RPM_RATE_UPDATE_HZ ((LOOP_FREQUENCY_HZ + RPM_RATE_TICKS_PER_UPDATE / 2) / RPM_RATE_TICKS_PER_UPDATE)

// Window of N = 2m + 1 samples, j \in [-m,...,m]
#define RPM_RATE_SG_HALF_WINDOW 8
#define RPM_RATE_SG_WINDOW (2 * RPM_RATE_SG_HALF_WINDOW + 1)
// Scaling factor: sum_j j^2 over j \in [-m,...,m] = m(m + 1)(2m + 1) / 3, 408 for m = 8
#define RPM_RATE_SG_SUM_J2 \
    (RPM_RATE_SG_HALF_WINDOW * (RPM_RATE_SG_HALF_WINDOW + 1) * (2 * RPM_RATE_SG_HALF_WINDOW + 1) / 3)

// Sanity check for corrupted periods, single outliers are replaced by the previous sample, a persistent jump (e.g. after a re-sync)
// restarts the window on the measurement. The bound also keeps the fit sum below in int32. Zero-crossing jitter
// alone moves the samples by ~1300 eRPM RMS (peaks ~4800) at 13 000 rpm, so the gate is wider than the 4096 eRPM
// of the former filter.
#ifndef RPM_RATE_OUTLIER_Q4
#define RPM_RATE_OUTLIER_Q4 131072L // 8192 eRPM
#endif
#define RPM_RATE_MAX_OUTLIERS 8     // consecutive outliers tolerated before the restart
// Rate limit, also keeps c1 * RPM_RATE_UPDATE_HZ in int32
#define RPM_RATE_MAX_ERPM_RATE 4000000L // eRPM/s
#define RPM_RATE_MAX_C1_Q8 (RPM_RATE_MAX_ERPM_RATE * 256 / RPM_RATE_UPDATE_HZ)

// Ring buffer of the last N measurements z in 1/16 eRPM. window_q4[oldest] is j = -m and is overwritten next.
static int32_t window_q4[RPM_RATE_SG_WINDOW];
static uint8_t oldest;
static int32_t newest_q4;                // last sample entered, for the outlier gate
static uint8_t outlier_count;            // consecutive replaced samples
static int16_t erpm_rate_telem;          // eRPM/s / 100, as sent in the telemetry frame

// Measurement z in Q4: 16 * 60 * 10^6 / T = 960e6 / T for T in us per electrical revolution (T > 0). 960e6 fits in
// uint32 and the quotient in int32 for any T >= 1. This is one software division, the G0 has no hardware divider.
static int32_t periodToErpmQ4(int e_com_time_us)
{
    return (int32_t)(960000000UL / (uint32_t)e_com_time_us); // 60e6 us per minute * 16
}

// Restart with every sample of the window on the current measurement, i.e. at a rate of zero.
// Called by main.c while the motor is not settled, and below after a persistent run of outliers.
void rpmRateReset(int e_com_time_us)
{
    int32_t z_q4 = (e_com_time_us > 0) ? periodToErpmQ4(e_com_time_us) : 0;
    for (uint8_t i = 0; i < RPM_RATE_SG_WINDOW; i++) {
        window_q4[i] = z_q4;
    }
    oldest = 0;
    newest_q4 = z_q4;
    outlier_count = 0;
    erpm_rate_telem = 0;
}

// One estimation step per 1 kHz sample: enter the sample, fit the window, differentiate at its centre
void rpmRateUpdate(int e_com_time_us)
{
    if (e_com_time_us <= 0) {
        return;
    }

    // Gate, see RPM_RATE_OUTLIER_Q4
    int32_t z_q4 = periodToErpmQ4(e_com_time_us);
    int32_t jump_q4 = z_q4 - newest_q4;
    if (jump_q4 > RPM_RATE_OUTLIER_Q4 || jump_q4 < -RPM_RATE_OUTLIER_Q4) {
        if (++outlier_count > RPM_RATE_MAX_OUTLIERS) {
            rpmRateReset(e_com_time_us);
            return;
        }
        z_q4 = newest_q4;
    } else {
        outlier_count = 0;
    }

    // Enter the sample in place of the oldest one, which moves the window on by one sample
    window_q4[oldest] = z_q4;
    newest_q4 = z_q4;
    oldest = (oldest + 1 == RPM_RATE_SG_WINDOW) ? 0 : oldest + 1;

    // Step 1: fit p(j) = c0 + c1*j + c2*j^2 to the window, of which only c1 = sum_j j*z_j / sum_j j^2 is needed.
    // Since sum_j j = 0, summing j * (z_j - z_0) gives the same result, and the gate bounds every term:
    // |z_j - z_0| <= 8 * 131072, so |sum| <= 72 * 8 * 131072 < 2^27 and the * 16 below stays in int32.
    uint8_t centre = oldest + RPM_RATE_SG_HALF_WINDOW;
    if (centre >= RPM_RATE_SG_WINDOW) {
        centre -= RPM_RATE_SG_WINDOW;
    }
    const int32_t z_centre_q4 = window_q4[centre];
    int32_t sum_jz_q4 = 0;
    uint8_t i = oldest; // j = -m
    // Kept as a loop: -O3 would unroll it to ~850 bytes, which overflows the 27 KB flash of the E230 targets
#pragma GCC unroll 1
    for (int32_t j = -RPM_RATE_SG_HALF_WINDOW; j <= RPM_RATE_SG_HALF_WINDOW; j++) {
        sum_jz_q4 += j * (window_q4[i] - z_centre_q4);
        i = (i + 1 == RPM_RATE_SG_WINDOW) ? 0 : i + 1;
    }
    // c1 in 1/256 eRPM per sample (the sum is in 1/16 eRPM, times 16)
    int32_t c1_q8 = sum_jz_q4 * 16 / RPM_RATE_SG_SUM_J2;

    // Step 2: differentiate p at the centre, p'(0) = c1 per sample, and convert to per second: rate = c1 / dt
    if (c1_q8 > RPM_RATE_MAX_C1_Q8) {
        c1_q8 = RPM_RATE_MAX_C1_Q8;
    }
    if (c1_q8 < -RPM_RATE_MAX_C1_Q8) {
        c1_q8 = -RPM_RATE_MAX_C1_Q8;
    }
    int32_t erpm_rate = c1_q8 * RPM_RATE_UPDATE_HZ / 256; // eRPM/s

    // Telemetry value rate / 100, rounded to nearest (truncation read up to 99 eRPM/s toward zero), saturated to
    // the int16 range of the frame
    int32_t telem = (erpm_rate >= 0 ? erpm_rate + 50 : erpm_rate - 50) / 100;
    if (telem > INT16_MAX) {
        telem = INT16_MAX;
    }
    if (telem < INT16_MIN) {
        telem = INT16_MIN;
    }
    erpm_rate_telem = (int16_t)telem;
}

// Latest estimate in units of 100 eRPM/s, held between updates and 0 while the motor is not settled.
// It describes the rotor half a window (8.4 ms) before the newest sample.
// makeTelemPackage() sends it as a big-endian two's-complement int16.
int16_t getEstimatedRPMRateForTelemetry(void)
{
    return erpm_rate_telem;
}
