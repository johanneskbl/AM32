/*
 * rpm_rate.h
 *
 * Created on: 2026-10-26
 * Estimate of the eRPM rate: time derivative of the electrical RPM for telemetry, fitted to the zero-crossing times
 * recorded by zc_capture.c.
 */

#ifndef RPM_RATE_H_
#define RPM_RATE_H_

#include <stdint.h>

void rpmRateReset(void);  // motor not settled: drop the crossings, rate 0
void rpmRateUpdate(void); // every 1 kHz block while settled: take the new crossings into the fit window
int16_t getEstimatedRPMRateForTelemetry(void); // fits on demand; eRPM rate / 100, so 1 means 100 eRPM/s

#endif /* RPM_RATE_H_ */
