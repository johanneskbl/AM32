/*
 * zc_speed.c
 *
 * Created on: 2026-11-09
 * Rotor speed from the zero-crossing (ZC) times that the commutation routines log in zc_log (main.c), and the
 * crossing processing shared with the eRPM rate fit (rpm_rate.c).
 *
 * Why the crossings and not e_com_time
 *   e_com_time averages over one electrical revolution and is smoothed by the commutation_interval filter.
 *   But this loses information since the timestamps are differenced. Detecting the zero crossings directly
 *   preserves the full information content including timing jitter. At the same delay a far less noisy
 *   estimate is possible.
 *
 * 1. Crossing processing, shared by the speed estimate and the rate fit
 *    zcSpeedProcess() takes the new records from zc_log in every main-loop pass.
 *    Take timestamps t_k [us] from the free-running UTILITY_TIMER, unwrapped to a 32 bit integer.
 *    Note, this is not equal to the sum of thiszctime since INTERVAL_TIMER is reset at every crossing and
 *    loses a tick about every fifth time.
 *    A run of consecutive crossings ends when the log was overrun, at a change of direction, at a step that does not
 *    follow, at an interval outside (0, 22.5 ms] and when the motor stops. Both estimates start over then.
 *    It does not end when the detection changes between comparator and polling (~1430 rpm up, ~1140 rpm down
 *    with 7 pole pairs, higher with bi_direction):
 *    polling finds the crossings 0.1-0.2 ms later, which shows as a blip of ~600-1200 rpm/s for one rate window, while
 *    a restart would hold the rate for ~90 ms right where the rotor spins down or up.
 *    Sanity bound to correct timestamp without using memory:
 *    t_k is replaced by the midpoint of its neighbours when
 *      |t_k - (t_{k-1} + t_{k+1}) / 2| > interval / 4,   interval = (t_{k+3} - t_{k-3}) / 6
 *    It reads only the neighbouring times and doesn't use stored estimates, so it cannot lock out. The crossing three
 *    back is checked, corrected in place in the history and handed to the rate fit (rpm_rate.c).
 *
 * 2. Speed estimate, the measurement of the RPM filter (rpm_filter.c) for the speed PID and the telemetry RPM
 *    Mean electrical period over m revolutions, from the newest six crossings and the six crossings 6m before them:
 *      D = sum_{i=0..5} (t_{k-i} - t_{k-i-6m}),   eRPM = 60e6 / (D / 6m) = 360e6 m / D
 *    In simple terms, the formula averages over the zero-crossing intervals. Written in terms of the individual
 *    crossing intervals τ, D is a weighted sum with trapezoidal weights:
 *      weight  1 2 3 4 5 6 6 6 ... 6 6 5 4 3 2 1      (6m + 5 intervals, sum 36m)
 *              newest ─────────────────── oldest
 *    So the result is the mean electrical period over the last m revolutions. Mathematically that is two moving
 *    averages in a row, a linear-phase FIR filter with no recursive state:
 *      1. 6 crossings (make one electrical revolution):
 *         Each difference starts and ends on the same comparator input and edge, so the per-phase timing offsets
 *         cancel exactly, and the jitter of the end crossings averages down by sqrt(6).
 *      2. 6m crossings (m electrical revolutions):
 *         The longer the baseline, the smaller the jitter relative to the period.
 *    A moving average lags a ramp by half its span, so m is reduced until the baseline spans at most 2.5 ms: with
 *    7 pole pairs m = 1 below ~6900 rpm, 2 up to ~10 300 rpm, 3 above. The delay is (6m + 5) / 2 crossings, 1.2-1.8 ms
 *    above 6000 rpm and 8.7 ms at 900 rpm. The smoothing is left to the alpha-beta filter (15 Hz), which follows a
 *    ramp without lag; rpm_filter.c compares the averaging spans and filter frequencies.
 *    The newest three crossings enter before their check, to keep the delay short.
 *    The filter gets measured_rpm while a run has fewer than 12 crossings (start-up, after a break) and on targets
 *    without the microsecond timer (NXP).
 *
 * Integer arithmetic only: ~150 cycles per crossing and one division per pass with new crossings, in the main loop.
 */

#include "zc_speed.h"
#include "rpm_rate.h"

// Crossing history, a power of two: holds the 6 * ZC_SPEED_MAX_REVS + 6 crossings of the longest baseline. Within the
// 2.5 ms span m is at most 4 up to ~17 000 rpm with 7 pole pairs.
#define ZC_SPEED_HISTORY 32
#define ZC_SPEED_MASK (ZC_SPEED_HISTORY - 1)
#define ZC_SPEED_MAX_REVS ((ZC_SPEED_HISTORY - 6) / 6) // 4 electrical revolutions
#if (ZC_SPEED_HISTORY & ZC_SPEED_MASK) != 0 || ZC_SPEED_MAX_REVS < 1
#error "zc_speed.c: ZC_SPEED_HISTORY must be a power of two of at least 16"
#endif
// Longest baseline of the speed estimate. Its delay is about half of it; the alpha-beta filter does the smoothing.
#define ZC_SPEED_MAX_SPAN_US 2500UL
// Fewest crossings in a run for an estimate: the newest six and the six one revolution before them
#define ZC_SPEED_MIN_CROSSINGS 12
// Records taken per read, so one pass has a bounded run time (the rest waits in zc_log)
#define ZC_SPEED_READ_BATCH 16
#define ZC_SPEED_MAX_READS 2 // batches per pass: 32 records cover 3.4 ms at 13 300 rpm

static uint32_t zc_log_tail; // records taken from zc_log
// Crossing times [us] of the current run, crossing[cross_newest] newest; checked except for the newest three
static uint32_t crossing[ZC_SPEED_HISTORY];
static uint32_t cross_newest;
static uint32_t crossings;    // crossings in the current run, saturates at ZC_SPEED_HISTORY
static uint32_t clock_us;     // unwrapped UTILITY_TIMER time of the newest crossing
static uint32_t last_flags;   // flags of the newest crossing: step, direction
// Estimated eRPM, 0 without an estimate. Written by the main loop, read by the 1 kHz interrupt for the RPM filter.
static volatile uint32_t erpm_estimate;

// Step that follows step s in the direction of rotation: next_step_of[s] forward, next_step_of[s + 8] reverse
static const uint8_t next_step_of[16] = { 0, 2, 3, 4, 5, 6, 1, 0, 0, 6, 1, 2, 3, 4, 5, 0 };

// The crossings no longer belong to one run: start both estimates over
static void endRun(void)
{
    crossings = 0;
    erpm_estimate = 0;
    rpmRateRestart();
}

// Take up to max records in order from zc_log, without masking interrupts. Records the ISR overwrote before they
// were read end the run, so crossings on both sides of the loss are never treated as consecutive.
static uint8_t readCrossings(uint32_t* records, uint8_t max)
{
    uint32_t pending = zc_log_head - zc_log_tail;
    if (pending > ZC_LOG_SIZE) {
        endRun(); // continue with the oldest record still in the log
        zc_log_tail = zc_log_head - ZC_LOG_SIZE;
        pending = ZC_LOG_SIZE;
    }
    if (pending > max) {
        pending = max;
    }
    // Kept as a loop: -O3 would unroll it (~400 bytes)
#pragma GCC unroll 1
    for (uint32_t i = 0; i < pending; i++) {
        records[i] = zc_log[(zc_log_tail + i) & (ZC_LOG_SIZE - 1)];
    }
    // The ISR may have run during the copy. The slot of record zc_log_tail is reused by record zc_log_tail + SIZE,
    // after which the head is past zc_log_tail + SIZE: then the copy may hold newer records. Drop it; the next call
    // finds the overrun above.
    if (zc_log_head - zc_log_tail > ZC_LOG_SIZE) {
        return 0;
    }
    zc_log_tail += pending;
    return (uint8_t)pending;
}

// One crossing record: continuity, unwrapped time, the check of the crossing three back, which goes to the rate fit.
// Called for every crossing, so it stays out of line and small.
__attribute__((noinline)) static void addCrossing(uint32_t record)
{
    uint32_t flags = record & 0xFF;
    uint32_t time_us = record >> 16;
    uint32_t previous = last_flags;
    last_flags = flags;
    uint32_t interval = (time_us - clock_us) & 0xFFFF;
    uint32_t next_step = next_step_of[(previous & ZC_LOG_STEP_MASK) | ((flags & ZC_LOG_FORWARD) ^ ZC_LOG_FORWARD)];
    // Sanity checks: change of direction, missing next step,
    // timestamp outside interval (0, 22.5 ms] (interval - 1 wraps for 0)
    uint32_t run_break = ((flags ^ previous) & ZC_LOG_FORWARD) | ((flags & ZC_LOG_STEP_MASK) ^ next_step);
    if (run_break || crossings == 0 || interval - 1 >= ZC_MAX_INTERVAL_US) {
        endRun();
        clock_us = time_us; // a new run: only differences within it matter
    } else {
        clock_us += interval;
    }
    if (flags & ZC_LOG_STARTUP) {
        rpmRateRestart(); // the rate fit leaves out the start-up, the speed PID needs the estimate there already
    }

    cross_newest = (cross_newest + 1) & ZC_SPEED_MASK;
    crossing[cross_newest] = clock_us;
    if (crossings < ZC_SPEED_HISTORY) {
        crossings++;
    }
    if (crossings < 7) {
        return;
    }

    // Correction check for t_k: Compare with timestamps against its neighbours, in integers:
    //   |2 t_k - t_{k-1} - t_{k+1}| > (t_{k+3} - t_{k-3}) / 12
    uint32_t checked = (cross_newest - 3) & ZC_SPEED_MASK;
    uint32_t t_next = crossing[(cross_newest - 2) & ZC_SPEED_MASK];
    uint32_t t = crossing[checked];
    uint32_t t_previous = crossing[(cross_newest - 4) & ZC_SPEED_MASK];
    int32_t deviation = (int32_t)(t - t_previous) - (int32_t)(t_next - t);
    if (deviation < 0) {
        deviation = -deviation;
    }
    if (12 * deviation > (int32_t)(clock_us - crossing[(cross_newest - 6) & ZC_SPEED_MASK])) {
        t = t_previous + (t_next - t_previous) / 2;
        crossing[checked] = t; // in place, so the speed estimate and the later checks use the corrected time
    }
    rpmRateAddCrossing(t);
}

// Speed estimate from the history, see the top
static void updateSpeed(uint8_t pole_pairs)
{
    if (crossings < ZC_SPEED_MIN_CROSSINGS) {
        erpm_estimate = 0;
        return;
    }
    uint32_t m = (crossings - 6) / 6; // revolutions the run covers, at most ZC_SPEED_MAX_REVS
    if (m > pole_pairs && pole_pairs > 0) {
        m = pole_pairs;
    }
    uint32_t newest = crossing[cross_newest];
    while (m > 1 && newest - crossing[(cross_newest - 6 * m) & ZC_SPEED_MASK] > ZC_SPEED_MAX_SPAN_US) {
        m--;
    }
    // Each difference is at most 6 intervals of 22.5 ms in a run, so D < 2^20 and 360e6 m + D / 2 < 2^32 for m <= 4.
    // Every interval of a run is at least 1 us, so D > 0.
    uint32_t d = 0;
    for (uint32_t i = 0; i < 6; i++) {
        d += crossing[(cross_newest - i) & ZC_SPEED_MASK] - crossing[(cross_newest - i - 6 * m) & ZC_SPEED_MASK];
    }
    erpm_estimate = (360000000UL * m + d / 2) / d; // rounded
}

// Every main-loop pass: take the new crossings, then update the speed estimate once
void zcSpeedProcess(uint8_t pole_pairs, uint8_t running)
{
    if (!running) {
        // Stopped: no estimate, and the next crossing starts a new run
        zc_log_tail = zc_log_head;
        if (crossings != 0) {
            endRun();
        }
        return;
    }
    uint32_t records[ZC_SPEED_READ_BATCH];
    uint8_t taken = 0;
    for (uint8_t batch = 0; batch < ZC_SPEED_MAX_READS; batch++) {
        uint8_t count = readCrossings(records, ZC_SPEED_READ_BATCH);
        for (uint8_t i = 0; i < count; i++) {
            addCrossing(records[i]);
        }
        taken |= count;
        if (count < ZC_SPEED_READ_BATCH) {
            break;
        }
    }
    if (taken) {
        updateSpeed(pole_pairs);
    }
}

// Mechanical RPM for the RPM filter, rounded, or fallback_rpm without an estimate. One read of the estimate, which the
// main loop may update at any time.
int32_t getZcSpeedRPM(uint8_t pole_pairs, int32_t fallback_rpm)
{
    uint32_t erpm = erpm_estimate;
    if (erpm == 0 || pole_pairs == 0) {
        return fallback_rpm;
    }
    return (int32_t)((erpm + pole_pairs / 2) / pole_pairs);
}
