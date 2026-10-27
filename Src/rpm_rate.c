/*
 * rpm_rate.c
 *
 * Time derivative of the electrical RPM (eRPM rate, i.e. rotor acceleration), estimated from the
 * electrical-revolution period e_com_time [us] for the extended KISS telemetry frame.
 *
 * Measurement
 *   e_com_time T is the duration of the last electrical revolution (sum of the last six commutation intervals,
 *   see main.c), so every sample yields the speed averaged over that revolution:
 *     z = 60e6 / T  [eRPM]
 *
 * Model and filter
 *   Constant-rate model with x = eRPM and v = eRPM rate [eRPM/s], sampled every dt:
 *     x[k+1] = x[k] + v[k]*dt,   v[k+1] = v[k]
 *   tracked by an alpha-beta filter, which has the structure of the steady-state (fixed-gain) Kalman filter
 *   of this model:
 *     predict: x- = x + v*dt
 *     correct: r = z - x-,  x = x- + alpha*r,  v = v + (beta/dt)*r
 *
 * Gain design (pole placement, not noise covariances)
 *   The estimation error evolves with the characteristic polynomial  z^2 - (2 - alpha - beta)*z + (1 - alpha).
 *   alpha = 1 - p^2 and beta = (1 - p)^2 turn it into (z - p)^2, a critically damped double pole at
 *     p = exp(-2*pi*f*dt).
 *   The rate output is then exactly a backward difference followed by two identical unit-gain exponential
 *   moving averages:
 *     V(z) = (1 - z^-1) / dt * ((1 - p) / (1 - p*z^-1))^2 * Z(z)
 *   whose continuous-time equivalent is s / (1 + s/(2*pi*f))^2. So f is the pole frequency, not the -3 dB
 *   bandwidth. For white acceleration noise the Kalman-optimal beta at this alpha would be 2*(1 - p)^2
 *   (damping ~0.71). The critically damped choice halves it: a monotone rate estimate without overshoot,
 *   for slightly more delay.
 *
 * Resulting behaviour for f = 30 Hz and dt = 1.05 ms
 *   - delay of a slowly varying rate:          (2*p/(1 - p) + 1/2) * dt = 10.1 ms
 *   - gain relative to the exact derivative:   ~ 1 / (1 + (F/f)^2) at frequency F, -3 dB at 19.3 Hz
 *   - response to an acceleration step a:      a * (1 - (1 + t/tau)*exp(-t/tau)), tau = 5.3 ms,
 *                                              10-90 % in 17.8 ms, no overshoot, no steady-state error
 *   - white measurement noise sigma_z [eRPM]:  sigma_v = sigma_z * (1 - p)^1.5 / dt * sqrt(2 / (1 + p)^3)
 *                                              = 41.7/s * sigma_z (a raw difference gives 1347/s * sigma_z)
 *
 * Integer arithmetic only. Called from the main loop with the e_com_time samples taken at the fixed rate of the
 * 1 kHz block in tenKhzRoutine(), which keeps the divisions out of the interrupt.
 */

#include "rpm_rate.h"
#include "targets.h"

// The samples come from the 1 kHz block of tenKhzRoutine(), which runs every (PID_LOOP_DIVIDER + 1) interrupts,
// so dt = RPM_RATE_TICKS_PER_UPDATE / LOOP_FREQUENCY_HZ = 21 / 20000 s = 1.05 ms at 20 kHz
#define RPM_RATE_TICKS_PER_UPDATE (PID_LOOP_DIVIDER + 1)

// Gains in fixed point, alpha in Q11 (x 2048) and beta in Q16 (x 65536), for f = 30 Hz at dt = 1.05 ms.
// They are tied to this dt: with another dt the poles stay at p, so f scales with 1.05 ms / dt.
// Other choices of f (alpha / beta, delay): 15 Hz: 368 / 582, 20.7 ms; 20 Hz: 475 / 1001, 15.4 ms;
// 25 Hz: 575 / 1515, 12.2 ms; 40 Hz: 840 / 3526, 7.5 ms. The delay scales with 1/f, the noise with f^1.5.
#ifndef RPM_RATE_ALPHA_Q11
#define RPM_RATE_ALPHA_Q11 670 // alpha = 670 / 2048 = 0.3271 (design value 0.3269)
#endif
#ifndef RPM_RATE_BETA_Q16
#define RPM_RATE_BETA_Q16 2113 // beta = 2113 / 65536 = 0.03224
#endif

// beta/dt in Q10, applied to the Q4 residual: beta * LOOP_FREQUENCY_HZ / (16 * RPM_RATE_TICKS_PER_UPDATE)
// The velocity correction is (beta/dt) * r with r = residual_q4 / 16 eRPM, i.e. beta / (16 * dt) eRPM/s per Q4 count
// (= 1.919 for the gains above, 1965 in Q10). Computed at compile time in 64 bit to avoid overflow.
#define RPM_RATE_KV_Q10 ((int32_t)((RPM_RATE_BETA_Q16 * 1024ULL * LOOP_FREQUENCY_HZ) / (65536ULL * 16ULL * RPM_RATE_TICKS_PER_UPDATE)))

// A rotor cannot move 4096 eRPM away from the prediction within one update, so such a residual is a corrupted
// period. Single outliers are skipped, a persistent jump (e.g. after a re-sync) restarts the filter on the measurement.
// A sudden acceleration step a causes a peak residual of 2.26 * a * dt, so real motion trips this gate only above
// ~1.7e6 eRPM/s. Measurement noise stays below it (simulated sample-to-sample jumps < 1600 eRPM up to 13 300 rpm).
#define RPM_RATE_OUTLIER_Q4 65536L   // 4096 eRPM
#define RPM_RATE_MAX_OUTLIERS 8      // consecutive outliers tolerated before the restart
#define RPM_RATE_MAX_ERPM_RATE 4000000L // eRPM/s, keeps erpm_rate * 16 * RPM_RATE_TICKS_PER_UPDATE in int32 (< 1.35e9)

// Filter state: x = erpm_q4 / 16 [eRPM] and v = erpm_rate [eRPM/s]. x carries 4 fractional bits because the
// prediction step v*dt is small (1 eRPM per update at v = 950 eRPM/s) and would otherwise be lost to truncation.
static int32_t erpm_q4;                  // x, filtered eRPM in 1/16 eRPM
static int32_t erpm_rate;                // v, eRPM/s
static uint8_t outlier_count;            // consecutive gated (skipped) updates
static int16_t erpm_rate_telem;          // eRPM/s / 100, as sent in the telemetry frame

// Measurement z in Q4: 16 * 60e6 / T = 960e6 / T for T in us per electrical revolution (T > 0). 960e6 fits in
// uint32 and the quotient in int32 for any T >= 1. This is one software division, the G0 has no hardware divider.
static int32_t periodToErpmQ4(int e_com_time_us)
{
    return (int32_t)(960000000UL / (uint32_t)e_com_time_us); // 60e6 us per minute * 16
}

// Restart at rest on the current measurement (x = z, v = 0), so the next update only sees the motion of one dt
// instead of a jump from a stale state. v then settles to the true rate within ~30 ms (98 % of an acceleration step).
// Called by main.c while the motor is not settled, and below after a persistent run of outliers.
void rpmRateReset(int e_com_time_us)
{
    erpm_q4 = (e_com_time_us > 0) ? periodToErpmQ4(e_com_time_us) : 0;
    erpm_rate = 0;
    outlier_count = 0;
    erpm_rate_telem = 0;
}

// One filter step per 1 kHz sample, see the equations at the top of this file.
void rpmRateUpdate(int e_com_time_us)
{
    // Defensive: main.c only passes valid periods and calls rpmRateReset() for the rest
    if (e_com_time_us <= 0) {
        return;
    }

    // Predict x- = x + v*dt in Q4, where 16*v*dt = v * (16 * RPM_RATE_TICKS_PER_UPDATE) / LOOP_FREQUENCY_HZ.
    // The division truncates toward zero, which biases v by at most (1/16 eRPM) / dt = 60 eRPM/s (< 1 telemetry count).
    int32_t predicted_q4 = erpm_q4 + (erpm_rate * (16 * RPM_RATE_TICKS_PER_UPDATE)) / LOOP_FREQUENCY_HZ;
    // Innovation r = z - x- in Q4
    int32_t residual_q4 = periodToErpmQ4(e_com_time_us) - predicted_q4;
    // Gate: an implausible residual is not used. The filter coasts on its prediction with v unchanged and restarts on
    // the measurement once more than RPM_RATE_MAX_OUTLIERS samples in a row were rejected.
    if (residual_q4 > RPM_RATE_OUTLIER_Q4 || residual_q4 < -RPM_RATE_OUTLIER_Q4) {
        if (++outlier_count > RPM_RATE_MAX_OUTLIERS) {
            rpmRateReset(e_com_time_us);
        } else {
            erpm_q4 = predicted_q4;
        }
        return;
    }
    outlier_count = 0;

    // Correct x = x- + alpha*r (Q4 * Q11 / 2^11 = Q4) and v = v + (beta/dt)*r (Q4 count * Q10 gain / 2^10 = eRPM/s).
    // The gate bounds |r| <= 65536, so both products fit in int32. Dividing instead of shifting truncates toward zero
    // for either sign, so positive and negative residuals are treated symmetrically.
    erpm_q4 = predicted_q4 + (residual_q4 * RPM_RATE_ALPHA_Q11) / 2048;
    erpm_rate += (residual_q4 * RPM_RATE_KV_Q10) / 1024;
    if (erpm_rate > RPM_RATE_MAX_ERPM_RATE) {
        erpm_rate = RPM_RATE_MAX_ERPM_RATE;
    }
    if (erpm_rate < -RPM_RATE_MAX_ERPM_RATE) {
        erpm_rate = -RPM_RATE_MAX_ERPM_RATE;
    }

    // Telemetry value v / 100 (truncated toward zero), saturated to the int16 range of the frame (+-3.28e6 eRPM/s)
    int32_t telem = erpm_rate / 100;
    if (telem > INT16_MAX) {
        telem = INT16_MAX;
    }
    if (telem < INT16_MIN) {
        telem = INT16_MIN;
    }
    erpm_rate_telem = (int16_t)telem;
}

// Latest estimate in units of 100 eRPM/s, held between updates and 0 while the motor is not settled.
// makeTelemPackage() sends it as a big-endian two's-complement int16.
int16_t rpmRateTelemetry(void)
{
    return erpm_rate_telem;
}
