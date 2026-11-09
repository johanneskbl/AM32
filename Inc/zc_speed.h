/*
 * zc_speed.h
 *
 * Created on: 2026-11-09
 * Rotor speed from the zero crossings (ZC) that the commutation routines in main.c log in zc_log: the shared
 * crossing processing for the RPM filter (rpm_filter.c: speed PID input and telemetry RPM) and the eRPM rate fit
 * (rpm_rate.c).
 */

#ifndef ZC_SPEED_H_
#define ZC_SPEED_H_

#include <stdint.h>

// Zero-crossing log: a ring of the last ZC_LOG_SIZE detected crossings, written by interruptRoutine() and
// zcfoundroutine() with interrupts disabled, read by zcSpeedProcess(). One 32-bit record per crossing, written and
// read with a single access:
//   bits 31..16  UTILITY_TIMER count at the crossing [us]. That timer runs freely and is never written, so differences
//                of consecutive records are exact (mod 65536 us), unlike sums of thiszctime: INTERVAL_TIMER is
//                written to 0 at every crossing while its prescaler keeps running, which loses a tick about every
//                fifth crossing.
//   bits  7..0   flags below
// At 13 300 rpm and 7 pole pairs there are 9300 crossings/s, so the main loop may stall for 64 / 9300 = 6.9 ms before
// records are overwritten; zcSpeedProcess() then starts a new run.
#define ZC_LOG_SIZE 64 // power of two
#define ZC_LOG_STEP_MASK 0x07u // commutation step 1..6 at the crossing
#define ZC_LOG_FORWARD 0x08u   // forward == 1
#define ZC_LOG_STARTUP 0x10u   // zero_crosses < 100: start-up, stall restart, desync

extern volatile uint32_t zc_log[ZC_LOG_SIZE];
extern volatile uint32_t zc_log_head; // records written since boot, wraps

// Longest interval between consecutive crossings of one run: the stall restart in main.c fires at
// INTERVAL_TIMER_COUNT > 45000 (0.5 us)
#define ZC_MAX_INTERVAL_US 22500

void zcSpeedProcess(uint8_t pole_pairs, uint8_t running); // main loop, every pass: take the new crossings
int32_t getZcSpeedRPM(uint8_t pole_pairs, int32_t fallback_rpm); // mean speed of the last <= 2.5 ms of crossings

#endif /* ZC_SPEED_H_ */
