/*
 * rpm_rate.c
 *
 * Time derivative of the electrical RPM (eRPM rate, i.e. rotor acceleration) for the extended KISS telemetry frame,
 * fitted to the zero-crossing (ZC) times that the commutation ISRs record (zc_capture.c).
 *
 * Why the crossings and not e_com_time
 *   e_com_time averages over one electrical revolution and is smoothed by the commutation_interval filter.
 *   But this loses information since the timestamps are differenced. Detecting the zero crossings directly
 *   preserves the full information content including timing jitter. At the same delay a far less noisy
 *   estimate is possible.
 *
 * 1. Timestamp crossings through zcCaptureRead()
 *    Take timestamps t_k [us] from the free-running UTILITY_TIMER, unwrapped to a 32 bit integer.
 *    Note, this is not equal to the sum of thiszctime since INTERVAL_TIMER is reset at every crossing and
 *    loses a tick about every fifth time.
 *    Sanity bound to correct timestamp without using memory:
 *    t_k is replaced by the midpoint of its neighbours when
 *      |t_k - (t_{k-1} + t_{k+1}) / 2| > interval / 4,   interval = (t_{k+3} - t_{k-3}) / 6
 *    It reads only the raw neighbouring times and doesn't use stored data, so it cannot lock out.
 *    The window restarts at a record flagged GAP (ring overrun) or STARTUP, at a change of direction or detection
 *    mode (polling finds the crossing ~0.1 ms later), at a step that does not follow, and at an interval outside
 *    (0, 22.5 ms].
 * 2. Average electrical revolutions
 *    Sample down the crossings: P_j = sum of 6 consecutive checked crossing times which represent one electrical
 *    revolution with each comparator input, so the per-phase timing offsets cancel out.
 *    Importantly, for the sum of times,
 *      P_j = t_1 + t_2 + ... + t_6,
 *    every crossing enters once, so all six crossings contribute and their jitter averages down by √6,
 *    consecutive up and down movemet is canceled out. While for the period e_com_time,
 *      T = (t_2-t_1) + ... + (t_7-t_6) = t_7 − t_1,
 *    the intervals telescope, so the four crossings in between cancel exactly.
 *    Only the two boundary crossings are left, and their jitter adds.
 *    The window of this Savgol filter holds as many points needed to span the last 42 ms. Its sums
 *      S0 = sum y, S1 = sum u*y, S2 = sum u^2*y
 *    run along with it, with u = 0..n-1 and y = P - P_oldest counted from the oldest point,
 *    which therefore adds nothing. A new point is added in O(1). Dropping the oldest moves the origin
 *    to the next point (u -> u - 1, y -> y - d):
 *      S2' = S2 - 2 S1 + S0 - d * sum (u-1)^2,  S1' = S1 - S0 - d * sum (u-1),  S0' = S0 - n' d
 *    with the sums over u = 1..n' of the points that remain.
 * 3. Apply Savgol filter by fitting a polynomial
 *    Only runs when the telemetry frame asks for the rate and a new point arrived since the last time
 *    For the N = 2m + 1 newest points, calculate the weighted sums around the centre:
 *      w = u - u_c,
 *      x = y - y_c,
 *    and fit a quadratic polynomial
 *      x(w) = a + b w + c w^2
 *    by least squares.
 * 4. Derivate to calculate the eRPM rate
 *    The period at the window centre T = b / 6 [us], its change per revolution dT = c / 3 [us], are needed 
 *    to calculate the rate of change of the eRPM
 *      rate = d(eRPM)/dt = -60e12 * dT / T^3 = -eRPM^2 * (dT / T) / 60  [eRPM/s]
 *    The window is symmetric, so sum w = sum w^3 = 0 and the sums of w^2, w^4 have closed forms. With
 *      S_x = S0 - N y_c,  S_wx = S1 - u_c S0,  S_w2x = S2 - 2 u_c S1 + u_c^2 S0 - y_c sum w^2
 *    this becomes
 *      T = 2 S_wx / (N(N^2-1)),  E = 12 S_w2x - (N^2-1) S_x,  rate = -eRPM^2 * E / (24 (N^2-4) S_wx)
 *    It describes the rotor acceleration at the window centre, with 42 ms window coverage this means a delay of ~23 ms.
 * 5. Sanity check
 *    A fit is anomalous when the window is not full, when T is outside 120 us..135 ms, when |rate| exceeds
 *    RPM_RATE_MAX_ERPM_RATE, or when it differs from the last output by more than RPM_RATE_MAX_STEP per 1 kHz update
 *    since the last accepted value. Then the last output is kept for this frame and the check repeats with the next.
 *    The allowed difference grows with every held update, so a real change of the rate is accepted after a few updates
 *    and an anomaly cannot freeze the output.
 *
 * Integer arithmetic only, no loop over the window. rpmRateUpdate() takes the crossings at the 1 kHz block of
 * tenKhzRoutine() (~150 cycles per crossing), getEstimatedRPMRateForTelemetry() fits only when a frame is built
 * (~1700 cycles), both in the main loop.
 */

#include "rpm_rate.h"
#include "zc_capture.h"

// Fitting window
#define RPM_RATE_WINDOW_US 42000L
// Revolution points stored
#define RPM_RATE_MAX_REVS 96
// Fewest points for a fit
#define RPM_RATE_MIN_REVS 5
// Crossing interval limit: the stall restart in main.c fires at INTERVAL_TIMER_COUNT > 45000 (0.5 us)
#define RPM_RATE_MAX_INTERVAL_US 22500
// Plausible electrical period limit at the window centre, Q4 us: 120 us (500 000 eRPM, keeps eRPM << 12 in uint32) up to
// six of the longest intervals
#define RPM_RATE_MIN_PERIOD_Q4 (120UL * 16)
#define RPM_RATE_MAX_PERIOD_Q4 (6UL * RPM_RATE_MAX_INTERVAL_US * 16)
// Plausible rate limit: 2 000 000 eRPM/s = 286 000 rpm/s with 7 pole pairs
#ifndef RPM_RATE_MAX_ERPM_RATE
#define RPM_RATE_MAX_ERPM_RATE 2000000L
#endif
// Plausible change limit of the output per update (1.05 ms) since the last accepted value
#ifndef RPM_RATE_MAX_STEP
#define RPM_RATE_MAX_STEP 150000L // eRPM/s
#endif
#define RPM_RATE_MAX_HOLDS 255
#if RPM_RATE_MAX_HOLDS * RPM_RATE_MAX_STEP < 2 * RPM_RATE_MAX_ERPM_RATE
#error "rpm_rate.c: after the longest hold every plausible rate must be accepted again"
#endif
#if RPM_RATE_MAX_HOLDS * RPM_RATE_MAX_STEP > 2147483647 || RPM_RATE_MAX_ERPM_RATE > 1000000000
#error "rpm_rate.c: the allowed change and the rate difference must fit int32"
#endif
// Quotient bits of divideU64(). |rate| <= RPM_RATE_MAX_ERPM_RATE needs q = |rate| * 24 (N^2 - 4) / eRPM below
// 3.2e6 < 2^23 for every window (largest at N = 95 and 137 000 eRPM), so larger quotients are not plausible.
#define RPM_RATE_QUOTIENT_BITS 23
// Records taken per read, so one update has a bounded run time (the rest waits in the capture ring)
#define RPM_RATE_READ_BATCH 16
#define RPM_RATE_MAX_READS 2 // batches per update: 32 records cover 3.4 ms at 13 300 rpm

// Last crossing times [us] in a ring of 8, crossing[cross_newest] newest, and how many belong to the current run
static uint32_t crossing[8];
static uint32_t cross_newest;
static uint32_t crossings;
static uint32_t clock_us;     // unwrapped UTILITY_TIMER time of the newest crossing
static uint32_t last_flags;   // flags of the newest crossing: step, direction, detection mode
// Revolution points: the sum of the six checked crossing times collected so far, and the window of completed points
static uint32_t rev_sum;
static uint32_t rev_crossings;
static uint32_t rev[RPM_RATE_MAX_REVS];
static uint8_t rev_oldest;    // index of the oldest point in rev[], u = 0
static uint8_t revs;          // points in the window
static uint8_t window_full;   // the window spans 42 ms: a point has left it since the restart
static uint8_t new_point;     // a point arrived since the last fit
static int64_t s0, s1, s2;    // window sums of y, u*y and u^2*y, see the top
// Output
static int32_t erpm_rate;     // last output, eRPM/s
static uint8_t holds;         // updates since the last accepted fit (age of the output)
static int16_t erpm_rate_telem; // eRPM/s / 100, as sent in the telemetry frame

static void restartWindow(void)
{
    crossings = 0;
    rev_sum = 0;
    rev_crossings = 0;
    revs = 0;
    s0 = 0;
    s1 = 0;
    s2 = 0;
    window_full = 0;
    new_point = 0;
}

// Helper function for integer arithmetic to reduce memory usage
// x / 3 for an exact multiple of 3, by the inverse of 3 mod 2^32: one multiplication, while the G0 would call the
// division routine even for a constant
static uint32_t divideExactBy3(uint32_t x)
{
    return x * 0xAAAAAAABUL;
}

// sum of v^2 for v = 0..k-1 = (k-1) k (2k-1) / 6, exact for k <= 96
static uint32_t sumOfSquares(uint32_t k)
{
    return (k == 0) ? 0 : divideExactBy3((k - 1) * k / 2 * (2 * k - 1));
}

static uint8_t revIndex(uint8_t u)
{
    uint32_t i = (uint32_t)rev_oldest + u;
    return (uint8_t)((i >= RPM_RATE_MAX_REVS) ? i - RPM_RATE_MAX_REVS : i);
}

// num / d by shift and subtract for a quotient below 2^23 (num < d * 2^23, checked by the caller, d < 2^30). The libgcc
// 64-bit division would link ~5 KB including the C++ unwinder; this loop has 23 fixed iterations.
static uint32_t divideU64(uint64_t num, uint32_t d)
{
    uint32_t remainder = (uint32_t)(num >> RPM_RATE_QUOTIENT_BITS);   // < d
    uint32_t low = (uint32_t)num << (32 - RPM_RATE_QUOTIENT_BITS);    // the remaining bits, top-aligned
    uint32_t quotient = 0;
    for (uint8_t bit = 0; bit < RPM_RATE_QUOTIENT_BITS; bit++) {
        remainder = (remainder << 1) | (low >> 31); // < 2d < 2^31
        low <<= 1;
        quotient <<= 1;
        if (remainder >= d) {
            remainder -= d;
            quotient |= 1;
        }
    }
    return quotient;
}

// Drop the oldest point and move the origin of the sums to the next one, see the top. All points lie within
// 6 * 42 ms (252 000) of each other, so S0 < 2^25, S1 < 2^31, S2 < 2^37 and every product below fits int64.
static void dropOldest(void)
{
    uint8_t next = revIndex(1);
    uint32_t d = rev[next] - rev[rev_oldest];
    revs--;
    uint32_t n = revs;
    // d * sum (u-1) = d n (n-1) / 2 < 252 000 * 4465 and d * n fit 32 bit; only d * sum (u-1)^2 needs 64
    s2 = s2 - (s1 << 1) + s0 - (int64_t)d * sumOfSquares(n);
    s1 = s1 - s0 - (int64_t)(d * (n * (n - 1) / 2));
    s0 = s0 - (int64_t)(d * n);
    rev_oldest = next;
}

// Step that follows step s in the direction of rotation: next_step_of[s] forward, next_step_of[s + 8] reverse
static const uint8_t next_step_of[16] = { 0, 2, 3, 4, 5, 6, 1, 0, 0, 6, 1, 2, 3, 4, 5, 0 };

// Enter a completed revolution point, after dropping the points it pushes out of the window or the ring
__attribute__((noinline)) static void addRevolution(uint32_t point)
{
    while (revs > 0 && (revs == RPM_RATE_MAX_REVS || point - rev[rev_oldest] > 6UL * RPM_RATE_WINDOW_US)) {
        dropOldest(); // points are sums of six times
        window_full = 1;
    }
    if (revs == 0) {
        rev_oldest = 0;
        s0 = 0;
        s1 = 0;
        s2 = 0;
    }
    uint32_t u = revs;
    uint32_t y = (revs == 0) ? 0 : point - rev[rev_oldest]; // at most 252 000, so u^2 * y < 2^32
    rev[revIndex((uint8_t)u)] = point;
    s0 += y;
    s1 += u * y;
    s2 += u * u * y;
    revs++;
    new_point = 1;
}

// One crossing record: continuity, unwrapped time, the check of the crossing three back, revolution sum.
// Called for every crossing, so it stays out of line and small; the rarer work is in addRevolution().
__attribute__((noinline)) static void addCrossing(uint32_t record)
{
    uint32_t flags = ZC_CAPTURE_FLAGS(record);
    uint32_t previous = last_flags;
    last_flags = flags;
    uint32_t interval = (ZC_CAPTURE_TIME_US(record) - clock_us) & 0xFFFF;
    uint32_t next_step = next_step_of[(previous & ZC_CAPTURE_STEP_MASK) | ((flags & ZC_CAPTURE_FORWARD) ^ ZC_CAPTURE_FORWARD)];
    // Sanity checks: lost records, start-up, change of direction or detection mode, missing next step,
    // timestamp outside interval (0, 22.5 ms] (interval - 1 wraps for 0)
    uint32_t run_break = (flags & (ZC_CAPTURE_GAP | ZC_CAPTURE_STARTUP))
        | ((flags ^ previous) & (ZC_CAPTURE_FORWARD | ZC_CAPTURE_POLLING)) | ((flags & ZC_CAPTURE_STEP_MASK) ^ next_step);
    if (run_break || crossings == 0 || interval - 1 >= RPM_RATE_MAX_INTERVAL_US) {
        restartWindow();
        clock_us = ZC_CAPTURE_TIME_US(record); // a new run: only differences within it matter
    } else {
        clock_us += interval;
    }

    cross_newest = (cross_newest + 1) & 7;
    crossing[cross_newest] = clock_us;
    if (crossings < 7) {
        crossings++;
        if (crossings < 7) {
            return;
        }
    }

    // Correction check for t_k: Compare with timestamps against its neighbours, in integers:
    //   |2 t_k - t_{k-1} - t_{k+1}| > (t_{k+3} - t_{k-3}) / 12
    uint32_t t_next = crossing[(cross_newest - 2) & 7];
    uint32_t t = crossing[(cross_newest - 3) & 7];
    uint32_t t_previous = crossing[(cross_newest - 4) & 7];
    int32_t deviation = (int32_t)(t - t_previous) - (int32_t)(t_next - t);
    if (deviation < 0) {
        deviation = -deviation;
    }
    if (12 * deviation > (int32_t)(clock_us - crossing[(cross_newest - 6) & 7])) {
        t = t_previous + (t_next - t_previous) / 2;
    }
    rev_sum += t;
    // Add a revolution for six detected and validated crossings
    if (++rev_crossings == 6) {
        addRevolution(rev_sum);
        rev_sum = 0;
        rev_crossings = 0;
    }
}

// Quadratic fit over the window
// Note: The least squares problem is solved in a closed-form
static uint8_t fitRate(int32_t* rate)
{
    if (!window_full || revs < RPM_RATE_MIN_REVS) {
        return 0;
    }
    int32_t n = revs;
    int32_t u_low = 0;
    if ((n & 1) == 0) {
        n--;
        u_low = 1;
    }
    int32_t m = (n - 1) / 2;
    uint32_t u_c = (uint32_t)(u_low + m);
    uint32_t y_c = rev[revIndex((uint8_t)u_c)] - rev[rev_oldest];
    uint32_t s_w2 = divideExactBy3((uint32_t)(m * (m + 1) * (2 * m + 1))); // sum of w^2

    // Centred sums: All y >= 0, so S0, S1, S2 >= 0. u_c <= 95, y_c <= 252 000 and S0 < 2^25, so
    // N y_c and u_c S0 fit uint32; the other products need 64 bit.
    int64_t s_x = s0 - (int64_t)((uint32_t)n * y_c);
    int64_t s_wx64 = s1 - (int64_t)(u_c * (uint32_t)s0);
    int64_t s_w2x = s2 - (int64_t)(u_c << 1) * s1 + (int64_t)(u_c * u_c) * s0 - (int64_t)y_c * s_w2;
    if (s_wx64 <= 0 || s_wx64 >= (1LL << 30)) {
        return 0; // time must advance; at most 2256 * 252 000 for points within the window
    }
    uint32_t s_wx = (uint32_t)s_wx64;

    // Period at the centre T = 2 S_wx / (N(N^2 - 1)) in Q4 us. 32 * S_wx would overflow, so the division is split:
    // with S_wx = q * D + r, 32 * S_wx / D = 32 q + 32 r / D, r < D <= 857 280.
    uint32_t d = (uint32_t)(n * (n * n - 1));
    uint32_t period_q4 = (s_wx / d) * 32 + (s_wx % d) * 32 / d;
    if (period_q4 < RPM_RATE_MIN_PERIOD_Q4 || period_q4 > RPM_RATE_MAX_PERIOD_Q4) {
        return 0;
    }
    uint32_t erpm = 960000000UL / period_q4; // 60e6 * 16 / T, below 2^19

    // Derivate polynomial at the center
    // rate = -eRPM^2 * E / (24 (N^2 - 4) S_wx), evaluated as -(|E| * eRPM / S_wx) * (eRPM / (24 (N^2 - 4))) with
    // the sign of E. |S_w2x| < 2^39 and |S_x| < 2^26 even for the worst sums, so |E| < 2^43 and |E| * eRPM < 2^62.
    // 12 S_w2x by shifts, done unsigned because S_w2x can be negative
    int64_t e = (int64_t)(((uint64_t)s_w2x << 3) + ((uint64_t)s_w2x << 2)) - (int64_t)(n * n - 1) * s_x;
    uint64_t product = (uint64_t)((e < 0) ? -e : e) * erpm;
    if ((product >> RPM_RATE_QUOTIENT_BITS) >= s_wx) {
        return 0; // curvature not plausible
    }
    uint32_t q = divideU64(product, s_wx);
    uint32_t k_q12 = (erpm << 12) / (uint32_t)(24 * (n * n - 4)); // eRPM / (24 (N^2 - 4)) in Q12
    uint64_t magnitude = ((uint64_t)q * k_q12 + 2048) >> 12;      // rounded
    if (magnitude > RPM_RATE_MAX_ERPM_RATE) {
        return 0;
    }
    *rate = (e < 0) ? (int32_t)magnitude : -(int32_t)magnitude;
    return 1;
}

// Restart while the motor is not settled: drop the queued crossings, rate 0. The first plausible fit afterwards has
// no earlier value to be compared with, hence the full hold count.
void rpmRateReset(void)
{
    uint32_t records[RPM_RATE_READ_BATCH];
    for (uint8_t batch = 0; batch < RPM_RATE_MAX_READS; batch++) {
        if (zcCaptureRead(records, RPM_RATE_READ_BATCH) < RPM_RATE_READ_BATCH) {
            break;
        }
    }
    restartWindow();
    erpm_rate = 0;
    erpm_rate_telem = 0;
    holds = RPM_RATE_MAX_HOLDS;
}

// Every 1 kHz block while the motor is settled: take the new crossings into the window. The fit itself runs only when
// the telemetry frame asks for the rate through getEstimatedRPMRateForTelemetry()
void rpmRateUpdate(void)
{
    uint32_t records[RPM_RATE_READ_BATCH];
    for (uint8_t batch = 0; batch < RPM_RATE_MAX_READS; batch++) {
        uint8_t count = zcCaptureRead(records, RPM_RATE_READ_BATCH);
        for (uint8_t i = 0; i < count; i++) {
            addCrossing(records[i]);
        }
        if (count < RPM_RATE_READ_BATCH) {
            break;
        }
    }
    if (holds < RPM_RATE_MAX_HOLDS) {
        holds++; // age of the output in updates, for the allowed change
    }
}

// Latest estimate in units of 100 eRPM/s
int16_t getEstimatedRPMRateForTelemetry(void)
{
    if (!new_point) {
        return erpm_rate_telem;
    }
    new_point = 0;
    int32_t rate;
    if (!fitRate(&rate)) {
        return erpm_rate_telem;
    }
    // The allowed change grows with the updates since the last accepted value. A fit always follows an update, so
    // holds >= 1 here.
    int32_t allowed = RPM_RATE_MAX_STEP * holds;
    int32_t change = rate - erpm_rate;
    if (change > allowed || change < -allowed) {
        return erpm_rate_telem;
    }
    erpm_rate = rate;
    holds = 0;

    // Telemetry value rate / 100, rounded to nearest, saturated to the int16 range of the frame
    int32_t telem = (rate >= 0 ? rate + 50 : rate - 50) / 100;
    if (telem > INT16_MAX) {
        telem = INT16_MAX;
    }
    if (telem < INT16_MIN) {
        telem = INT16_MIN;
    }
    erpm_rate_telem = (int16_t)telem;
    return erpm_rate_telem;
}
