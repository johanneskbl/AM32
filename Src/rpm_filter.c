/*
 * rpm_filter.c
 * Created on: 2026-11-06
 *
 * Low-pass filtered mechanical RPM for the speed PID and the telemetry frame:
 * Implementation of the position state x of an alpha-beta tracking filter on the zero-crossing speed (zc_speed.c,
 * the mean speed of the last <= 2.5 ms of crossings), or on measured_rpm while that has no estimate,
 * updated once per 1 kHz speed loop tick (dt = 1.05 ms at 20 kHz).
 *
 * Model and filter
 *   Constant-rate model with the measurement z = zero-crossing speed, the position state x = RPM
 *   and the velocity state v = RPM change per tick:
 *     predict:     x̂⁻_k = x̂_{k−1} + Δt·v̂_{k−1}
 *     innovation:  r_k  = z_k − x̂⁻_k
 *     correct:     x̂_k  = x̂⁻_k + α·r_k
 *                  v̂_k  = v̂_{k−1} + (β/Δt)·r_k
 * 
 *   The velocity state v is only internal: it lets x follow a constant
 *   acceleration without lag, which makes this more accurate than a moving average filter.
 *
 * Gain design
 *   Steady-state Kalman gains (Kalata relation) for white acceleration and white measurement noise:
 *     alpha = 1 - p^2,  beta = 2*(1 - p)^2,  p = exp(-2*pi*f*dt)
 *   The frequency f is the tuning parameter, here f = 15 Hz. The poles are complex with a damping of ~0.71.
 *   The transfer function from z to x follows:
 *     X(z)/Z(z) = (alpha + (beta - alpha)*z^-1) / (1 - (2 - alpha - beta)*z^-1 + (1 - alpha)*z^-2)
 *   Choice of f and of the zero-crossing averaging (closed speed-loop model fitted to the 2026-11-09 ESC3 bench bags,
 *   telemetry RPM noise at normal / PWM-beat targets and RMS error against the true speed in a 5 Hz chirp):
 *     2.5 ms average + 15 Hz:  4.0 / 7.8 RPM, 18.6 RPM   <- selected
 *     2.5 ms average + 20 Hz:  4.8 / 10.1 RPM, 16.5 RPM
 *     5 ms average + 15 Hz:    2.5 / 5.6 RPM, 25.7 RPM
 *     10 ms average alone:     5.4 / 11.9 RPM, 30.7 RPM  (a moving average lags a ramp by half its span)
 *     e_com_time + 25 Hz:      20.4 / 18.4 RPM, 17.1 RPM (the former input)
 *   The short average keeps its delay at ~1.5 ms, and the filter smooths without lagging a ramp. The PWM-beat
 *   targets are where a commutation step spans a whole number of PWM periods (8571 RPM: 4, 11429 RPM: 3, ...): the
 *   detection error then varies at a few hundred Hz, which the filter removes before the telemetry samples it.
 * 
 * Integer arithmetic only (no division), so it runs in the 1 kHz block of tenKhzRoutine() next to the speed PID.
 */

#include "rpm_filter.h"

// Gains in fixed point, alpha in Q11 (x 2048) and beta in Q16 (x 65536), for f = 15 Hz at dt = 1.05 ms
// Other choices of f (alpha / beta): 20 Hz: 475 / 2003, 25 Hz: 575 / 3030. Lower f: less noise, more lag.
// Above ~29 Hz beta reaches 4096, which the residual clamp below no longer keeps in int32.
#ifndef RPM_FILTER_ALPHA_Q11
#define RPM_FILTER_ALPHA_Q11 368 // alpha = 368 / 2048 = 0.1797
#endif
#ifndef RPM_FILTER_BETA_Q16
#define RPM_FILTER_BETA_Q16 1164 // beta = 1164 / 65536 = 0.01776
#endif
#if RPM_FILTER_BETA_Q16 >= 4096
#error "RPM_FILTER_BETA_Q16 must stay below 4096, else residual * beta overflows int32"
#endif

// The rotor cannot move this far from the prediction within one tick (a 150 000 RPM/s step gives a residual of
// ~360 RPM), so a larger residual is a corrupted period. It is clamped instead of skipped, so the PID input never
// freezes. The bound also keeps residual * RPM_FILTER_BETA_Q16 below 2^31.
#define RPM_FILTER_MAX_RESIDUAL_Q8 (2048L * 256) // 2048 RPM
#define RPM_FILTER_MAX_V_Q8 (2048L * 256)        // 2048 RPM per tick, about 2e6 RPM/s

// Filter state in 1/256 RPM: the fractional bits keep the small corrections alpha*r and beta*r
static int32_t x_q8;                    // x, filtered RPM
static int32_t v_q8;                    // v, RPM change per tick
static volatile int32_t published_q8;   // x as last finished, read by the main loop for the telemetry frame

// Start on the measurement at rest (x = z, v = 0), so a later update only sees the motion of one tick
void rpmFilterReset(int32_t rpm)
{
    x_q8 = rpm * 256;
    v_q8 = 0;
    published_q8 = x_q8;
}

// One filter step per 1 kHz tick, see the equations at the top of this file
void rpmFilterUpdate(int32_t rpm)
{
    // Predict x- = x + v
    int32_t predicted_q8 = x_q8 + v_q8;
    // Residual r = z - x-, clamped (see RPM_FILTER_MAX_RESIDUAL_Q8)
    int32_t residual_q8 = rpm * 256 - predicted_q8;
    if (residual_q8 > RPM_FILTER_MAX_RESIDUAL_Q8) {
        residual_q8 = RPM_FILTER_MAX_RESIDUAL_Q8;
    }
    if (residual_q8 < -RPM_FILTER_MAX_RESIDUAL_Q8) {
        residual_q8 = -RPM_FILTER_MAX_RESIDUAL_Q8;
    }

    // Correct x = x- + alpha*r and v = v + beta*r. Dividing instead of shifting truncates toward zero for either
    // sign, so positive and negative residuals are treated symmetrically.
    x_q8 = predicted_q8 + residual_q8 * RPM_FILTER_ALPHA_Q11 / 2048;
    v_q8 += residual_q8 * RPM_FILTER_BETA_Q16 / 65536;
    if (v_q8 > RPM_FILTER_MAX_V_Q8) {
        v_q8 = RPM_FILTER_MAX_V_Q8;
    }
    if (v_q8 < -RPM_FILTER_MAX_V_Q8) {
        v_q8 = -RPM_FILTER_MAX_V_Q8;
    }
    published_q8 = x_q8;
}

// Filtered RPM for the speed PID, rounded to whole RPM
int32_t getFilteredRPM(void)
{
    return (x_q8 + 128) / 256;
}

// Filtered RPM as eRPM / 10 for the extended telemetry frame, rounded and limited to the uint16 field.
// x * pole_pairs / 10 with x in 1/256 RPM: (x_q8 * pole_pairs + 1280) / 2560.
uint16_t getFilteredRPMForTelemetry(uint8_t pole_pairs)
{
    int32_t rpm_q8 = published_q8;  // One read, the 1 kHz interrupt may update it at any time
    if (rpm_q8 <= 0) {
        return 0;
    }
    uint32_t erpm_div10 = ((uint32_t)rpm_q8 * pole_pairs + 1280) / 2560;
    return (erpm_div10 > 65535) ? 65535 : (uint16_t)erpm_div10;
}
