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

void rpmRateAddCrossing(uint32_t time_us); // zc_speed.c: the next checked crossing time of the current run
void rpmRateRestart(void);                 // zc_speed.c: a new run of crossings, or start-up
void rpmRateTick(uint8_t settled);         // every 1 kHz block: motor settled or not
int16_t getEstimatedRPMRateForTelemetry(void); // fits on demand; eRPM rate / 100, so 1 means 100 eRPM/s

#endif /* RPM_RATE_H_ */
