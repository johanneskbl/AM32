/*
 * rpm_filter.c
 * Created on: 2026-11-06
 *
 * Low-pass filtered mechanical RPM for the speed PID and the telemetry frame:
 * Implementation of the position state x of an alpha-beta tracking filter on measured_rpm,
 * updated once per 1 kHz speed loop tick (dt = 1.05 ms at 20 kHz).
 *
 * Model and filter
 *   Constant-rate model with the measurement z = measured_rpm, the position state x = RPM
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
 *   The frequency f is the tuning parameter, here f = 25 Hz. The poles are complex with a damping of ~0.71.
 *   The transfer function from z to x follows:
 *     X(z)/Z(z) = (alpha + (beta - alpha)*z^-1) / (1 - (2 - alpha - beta)*z^-1 + (1 - alpha)*z^-2)
 *   Compared with the critically damped double pole (beta = (1 - p)^2, f = 30 Hz) at the same 2.25 ms lag on
 *   rotor steps (simulation calibrated on the 29 Oct 2026 ESC3 bench data, 9 us zero-crossing jitter):
 *   noise 11.0 instead of 12.6 RPM at 8000 RPM (-13 %) and 36.0 instead of 37.1 RPM at 11000 RPM (-3 %), same step
 *   tracking and speed loop behaviour, frequency-response peak +1.7 dB instead of +1.0 dB.
 * 
 * Integer arithmetic only (no division), so it runs in the 1 kHz block of tenKhzRoutine() next to the speed PID.
 */

#include "rpm_filter.h"

// Gains in fixed point, alpha in Q11 (x 2048) and beta in Q16 (x 65536), for f = 25 Hz at dt = 1.05 ms
// Other choices of f (alpha / beta): 15 Hz: 368 / 1164, 20 Hz: 475 / 2003. Lower f: less noise, more lag.
// Above ~29 Hz beta reaches 4096, which the residual clamp below no longer keeps in int32.
#ifndef RPM_FILTER_ALPHA_Q11
#define RPM_FILTER_ALPHA_Q11 575 // alpha = 575 / 2048 = 0.2808
#endif
#ifndef RPM_FILTER_BETA_Q16
#define RPM_FILTER_BETA_Q16 3030 // beta = 3030 / 65536 = 0.04623
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
