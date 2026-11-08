/*
 * zc_capture.h
 *
 * Created on: 2026-11-08
 * Capture of the raw zero crossings (ZC) of the back-EMF:
 * Use cases:
 *   - input for estimating the eRPM rate in rpm_rate.c
 *   - create snapshots of the crossings for debugging on a test bench, enable by defining USE_ZC_CAPTURE,
 *     i.e, set `#define USE_ZC_CAPTURE` in the block for the specific target (e.g. TMOTOR_G071) in targets.h
 */

#ifndef ZC_CAPTURE_H_
#define ZC_CAPTURE_H_

#include <stdint.h>
#include "targets.h"

// One crossing is represented as a 32-bit , written and read with a single access, so the main loop never sees a half-written
// record:
//   bits 31..16  UTILITY_TIMER count at the crossing [us]. That timer runs freely and is never written, so differences
//                of consecutive records are exact (mod 65536 us). The sum of thiszctime is not: INTERVAL_TIMER is
//                written to 0 at every crossing while its prescaler keeps running, which loses a 0.5 us tick about
//                every fifth crossing, a random walk in the summed time.
//   bits 15..8   PWM timer count at the crossing / 32 (0.5 us at 64 MHz, saturated at 255, STM32 only, else 0), to
//                measure how the detection delay depends on the PWM phase (beats between crossing rate and PWM)
//   bits  7..0   flags below
// No crossings are recorded on NXP (no free-running microsecond timer) and WCH (no PRIMASK) targets.
#define ZC_CAPTURE_STEP_MASK 0x07u // commutation step 1..6 whose floating phase crossed zero
#define ZC_CAPTURE_FORWARD 0x08u   // forward == 1
#define ZC_CAPTURE_POLLING 0x10u   // found by zcfoundroutine() on the 20 kHz polling grid, else by the comparator ISR
#define ZC_CAPTURE_STARTUP 0x20u   // zero_crosses below ZC_CAPTURE_SETTLED_CROSSINGS: start-up, stall restart, desync
#define ZC_CAPTURE_GAP 0x80u       // set by zcCaptureRead(): records were lost right before this one (ring overrun)

#define ZC_CAPTURE_TIME_US(record) ((uint16_t)((record) >> 16))
#define ZC_CAPTURE_PWM_PHASE(record) ((uint8_t)((record) >> 8))
#define ZC_CAPTURE_FLAGS(record) ((uint8_t)(record))

void zcCaptureRecord(uint8_t step, uint8_t forward, uint32_t zero_crosses, uint8_t polling); // ISR side, per crossing
uint8_t zcCaptureRead(uint32_t* records, uint8_t max); // main loop: up to max records in order, returns the count

#ifdef USE_ZC_CAPTURE
// Debug snapshot of consecutive crossings: a copy of the records zcCaptureRead() hands to the rate fit
//
// Bench readout over SWD while the motor runs. openocd reads and writes RAM without halting the core. Never halt it
// while the motor is powered: the phase outputs would stay switched on one commutation step.
//   arm-none-eabi-nm -S obj/debug.elf | grep zc_snapshot        -> address A and size S (hex)
//   openocd -f obj/openocd.cfg -c "mwb 0xA 1" -c "sleep 3000" -c "dump_image zc.bin 0xA 0xS" -c "shutdown"
// state = 2 in the dump means the snapshot is complete. tools/scripts/zc_capture.py automates this and evaluates it.
#ifndef ZC_CAPTURE_SNAPSHOT_LENGTH
#define ZC_CAPTURE_SNAPSHOT_LENGTH 2048 // 8 KB of RAM, 0.39 s at 7500 rpm and 7 pole pairs
#endif
#define ZC_CAPTURE_VERSION 1

typedef struct {
    volatile uint8_t arm;      // offset 0: write 1 to record the next ZC_CAPTURE_SNAPSHOT_LENGTH crossings
    volatile uint8_t state;    // offset 1: 0 idle, 1 recording, 2 complete
    uint8_t max_pending;       // offset 2: most records queued when the main loop read them (> 64: overrun)
    uint8_t version;           // offset 3: ZC_CAPTURE_VERSION
    volatile uint16_t count;   // offset 4: records in entry[]
    uint16_t lost;             // offset 6: records lost to ring overruns while recording (saturates at 65535)
    uint16_t length;           // offset 8: ZC_CAPTURE_SNAPSHOT_LENGTH
    uint16_t pwm_arr;          // offset 10: PWM auto-reload when armed: the PWM period is pwm_arr + 1 timer counts
    uint32_t entry[ZC_CAPTURE_SNAPSHOT_LENGTH]; // offset 12
} zcCaptureSnapshot_t;

extern zcCaptureSnapshot_t zc_snapshot;
#endif /* USE_ZC_CAPTURE */

#endif /* ZC_CAPTURE_H_ */
