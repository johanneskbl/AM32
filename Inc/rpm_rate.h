/*
 * rpm_rate.h
 *
 * Created on: 2026-10-26
 * Estimate of the eRPM rate: time derivative of the electrical RPM for telemetry, fitted to the zero-crossing times
 * that the commutation routines in main.c log in zc_log.
 */

#ifndef RPM_RATE_H_
#define RPM_RATE_H_

#include <stdint.h>

// Zero-crossing log: a ring of the last ZC_LOG_SIZE detected crossings, written by interruptRoutine() and
// zcfoundroutine() with interrupts disabled, read by rpmRateUpdate(). One 32-bit record per crossing, written and
// read with a single access:
//   bits 31..16  UTILITY_TIMER count at the crossing [us]. That timer runs freely and is never written, so differences
//                of consecutive records are exact (mod 65536 us), unlike sums of thiszctime: INTERVAL_TIMER is
//                written to 0 at every crossing while its prescaler keeps running, which loses a tick about every
//                fifth crossing.
//   bits  7..0   flags below
// At 13 300 rpm and 7 pole pairs there are 9300 crossings/s, so the main loop may stall for 64 / 9300 = 6.9 ms before
// records are overwritten; rpmRateUpdate() then restarts the fit.
#define ZC_LOG_SIZE 64 // power of two
#define ZC_LOG_STEP_MASK 0x07u // commutation step 1..6 at the crossing
#define ZC_LOG_FORWARD 0x08u   // forward == 1
#define ZC_LOG_STARTUP 0x10u   // zero_crosses < 100: start-up, stall restart, desync

extern volatile uint32_t zc_log[ZC_LOG_SIZE];
extern volatile uint32_t zc_log_head; // records written since boot, wraps

void rpmRateReset(void);  // motor not settled: drop the logged crossings, rate 0
void rpmRateUpdate(void); // every 1 kHz block while settled: take the new crossings into the fit window
int16_t getEstimatedRPMRateForTelemetry(void); // fits on demand; eRPM rate / 100, so 1 means 100 eRPM/s

#endif /* RPM_RATE_H_ */
