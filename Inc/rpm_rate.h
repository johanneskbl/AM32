/*
 * rpm_rate.h
 *
 * Created on: 2026-10-26
 * Estimate of the eRPM rate: time derivative of the electrical RPM for telemetry.
 */

#ifndef RPM_RATE_H_
#define RPM_RATE_H_

#include <stdint.h>

void rpmRateReset(int e_com_time_us);
void rpmRateUpdate(int e_com_time_us);
int16_t rpmRateTelemetry(void); // eRPM rate / 100, so 1 means 100 eRPM/s

#endif /* RPM_RATE_H_ */
