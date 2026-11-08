/*
 * zc_capture.c
 *
 * Created on: 2026-11-08
 * Capture of the raw zero crossings, see zc_capture.h for the record format and the bench readout.
 *
 * Why is this needed?
 *  In the main loop, only the latest zero crossing is stored and folded into a filter, so each value is overwritten.
 * 
 *   Existing variable	        What it holds
 *   thiszctime, lastzctime:	Raw intervals of the last two crossings (0.5 µs ticks), updated every crossing
 *   commutation_interval:	    IIR mix of those intervals, updated every step
 *   commutation_intervals[6]:	One IIR value per step and summed into e_com_time, overwritten every revolution
 * 
 *   The rate update runs every 1.05 ms, so reading these from the main loop would see only 1 in 6 to 10 crossings.
 *   This down-sampling is limiting the performance of filters (especially for the RPM rate).
 *   To make use of every crossing, each one has to be stored the moment it happens, inside the ISR.
 *   That is what the zcCaptureRecord() does: it writes one 32-bit word into a 64-entry ring that the
 *   main loop empties.
 * 
 * Data path
 *   ISR:       interruptRoutine() and zcfoundroutine() call zcCaptureRecord() once per detected crossing, after the
 *              interval and commutation timers are set, so the commutation timing does not change. It stores one
 *              32-bit record in a ring buffer and advances the head: ~60 cycles (~1 us at 64 MHz) including the call,
 *              no loop and no division.
 *   Main loop: zcCaptureRead() takes the records in order without masking interrupts and checks afterwards that the
 *              ISR has not overwritten what was read. Lost records are counted and the next record is flagged GAP,
 *              never passed on as consecutive.
 *
 * Writers
 *   Only one record is written at a time: interruptRoutine() runs with interrupts masked, zcfoundroutine() runs in the
 *   20 kHz interrupt while the comparator interrupt is off, and the stall restart in the main loop calls
 *   zcfoundroutine() as well, which the 20 kHz interrupt can preempt. zcCaptureRecord() masks interrupts itself and
 *   restores the caller's state, so this holds for every caller.
 */

#include "zc_capture.h"
#include "functions.h"

#if defined(NXP) || defined(WCH)
// No free-running microsecond timer (NXP) or no PRIMASK (WCH): nothing is recorded, the rate fit stays without input
void zcCaptureRecord(uint8_t step, uint8_t forward, uint32_t zero_crosses, uint8_t polling)
{
}

uint8_t zcCaptureRead(uint32_t* records, uint8_t max)
{
    return 0;
}
#else

// Crossings before this count belong to the start-up, as for the rpm rate in tenKhzRoutine()
#define ZC_CAPTURE_SETTLED_CROSSINGS 100
// PWM timer count to the 8-bit phase field: / 32 keeps auto-reload values up to 8191 (8 kHz PWM at 64 MHz) in range
#define ZC_CAPTURE_PWM_SHIFT 5

// Ring buffer from the ISR to the main loop, read every 1.05 ms by rpm_rate.c. At 13 300 rpm and 7 pole pairs there
// are 9300 crossings/s, so the main loop may stall for 64 / 9300 = 6.9 ms before records are lost.
#define ZC_CAPTURE_RING_SIZE 64
#if (ZC_CAPTURE_RING_SIZE & (ZC_CAPTURE_RING_SIZE - 1)) != 0 || ZC_CAPTURE_RING_SIZE > 255
#error "zc_capture.c: ZC_CAPTURE_RING_SIZE must be a power of two below 256 (max_pending is 8 bit)"
#endif
#define ZC_CAPTURE_RING_MASK (ZC_CAPTURE_RING_SIZE - 1)

static volatile uint32_t ring[ZC_CAPTURE_RING_SIZE];
// Records written since boot, free-running (wraps after 2^32 records, > 5 days at 9300/s). Only the ISR writes it,
// and only after the record is stored, so every record below the head is complete.
static volatile uint32_t ring_head;
static uint32_t ring_tail;  // records the main loop has taken from the ring
static uint8_t gap_pending; // records were lost: flag the next record handed out

#ifdef USE_ZC_CAPTURE
#if ZC_CAPTURE_SNAPSHOT_LENGTH < 1 || ZC_CAPTURE_SNAPSHOT_LENGTH > 65535
#error "zc_capture.c: ZC_CAPTURE_SNAPSHOT_LENGTH must fit the 16-bit count"
#endif
// In .bss, zero at boot: idle until armed. Not static, so the debugger finds it by name.
zcCaptureSnapshot_t zc_snapshot;

// Arming (by the debugger) starts the snapshot with the next record handed out
static void snapshotService(uint32_t pending)
{
    if (zc_snapshot.arm) {
        zc_snapshot.arm = 0;
        zc_snapshot.count = 0;
        zc_snapshot.lost = 0;
        zc_snapshot.max_pending = 0;
        zc_snapshot.version = ZC_CAPTURE_VERSION;
        zc_snapshot.length = ZC_CAPTURE_SNAPSHOT_LENGTH;
#ifdef STMICRO
        zc_snapshot.pwm_arr = (uint16_t)TIM1->ARR;
#endif
        zc_snapshot.state = 1;
    }
    if (zc_snapshot.state == 1 && pending > zc_snapshot.max_pending) {
        zc_snapshot.max_pending = (uint8_t)pending;
    }
}

// Append while recording. The count is checked before the write because the debugger can change it.
static void snapshotAppend(uint32_t record)
{
    if (zc_snapshot.state != 1) {
        return;
    }
    uint16_t count = zc_snapshot.count;
    if (count >= ZC_CAPTURE_SNAPSHOT_LENGTH) {
        zc_snapshot.state = 2;
        return;
    }
    zc_snapshot.entry[count] = record;
    zc_snapshot.count = count + 1;
    if (count + 1 >= ZC_CAPTURE_SNAPSHOT_LENGTH) {
        zc_snapshot.state = 2;
    }
}
#endif /* USE_ZC_CAPTURE */

void zcCaptureRecord(uint8_t step, uint8_t forward, uint32_t zero_crosses, uint8_t polling)
{
    uint32_t primask = __get_PRIMASK();
    __disable_irq();
    // Timers first, at a fixed delay after the caller read thiszctime
    uint32_t time_us = get_timer_us16();
    uint32_t pwm_phase = 0;
#ifdef STMICRO
    pwm_phase = TIM1->CNT >> ZC_CAPTURE_PWM_SHIFT;
    if (pwm_phase > 0xFF) {
        pwm_phase = 0xFF;
    }
#endif
    uint32_t flags = step & ZC_CAPTURE_STEP_MASK;
    if (forward) {
        flags |= ZC_CAPTURE_FORWARD;
    }
    if (polling) {
        flags |= ZC_CAPTURE_POLLING;
    }
    if (zero_crosses < ZC_CAPTURE_SETTLED_CROSSINGS) {
        flags |= ZC_CAPTURE_STARTUP;
    }
    uint32_t record = (time_us << 16) | (pwm_phase << 8) | flags;
    uint32_t head = ring_head;
    ring[head & ZC_CAPTURE_RING_MASK] = record;
    ring_head = head + 1; // publish after the record is stored
    __set_PRIMASK(primask);
}

uint8_t zcCaptureRead(uint32_t* records, uint8_t max)
{
    uint32_t pending = ring_head - ring_tail;
#ifdef USE_ZC_CAPTURE
    snapshotService((pending > 0xFF) ? 0xFF : pending); // first, so a loss found in this read counts in the snapshot
#endif
    if (pending > ZC_CAPTURE_RING_SIZE) {
        // The ISR has overwritten records that were not read yet: continue with the oldest one still in the ring
#ifdef USE_ZC_CAPTURE
        if (zc_snapshot.state == 1) {
            uint32_t total = zc_snapshot.lost + (pending - ZC_CAPTURE_RING_SIZE);
            zc_snapshot.lost = (total > 0xFFFF) ? 0xFFFF : (uint16_t)total;
        }
#endif
        gap_pending = 1;
        ring_tail = ring_head - ZC_CAPTURE_RING_SIZE;
        pending = ZC_CAPTURE_RING_SIZE;
    }
    if (pending > max) {
        pending = max;
    }
    for (uint32_t i = 0; i < pending; i++) {
        records[i] = ring[(ring_tail + i) & ZC_CAPTURE_RING_MASK];
    }
    // The ISR may have run during the copy. The slot of record ring_tail is reused by record ring_tail + RING_SIZE,
    // after which the head is past ring_tail + RING_SIZE: then the copy may hold newer records. Drop it; the next call
    // finds the overrun above and flags the gap.
    if (ring_head - ring_tail > ZC_CAPTURE_RING_SIZE) {
        return 0;
    }
    ring_tail += pending;
    if (pending > 0 && gap_pending) {
        records[0] |= ZC_CAPTURE_GAP;
        gap_pending = 0;
    }
#ifdef USE_ZC_CAPTURE
    for (uint32_t i = 0; i < pending; i++) {
        snapshotAppend(records[i]);
    }
#endif
    return (uint8_t)pending;
}

#endif /* NXP || WCH */
