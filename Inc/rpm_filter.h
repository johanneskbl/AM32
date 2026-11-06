/*
 * rpm_filter.h
 *
 * Created on: 2026-11-06
 * Low-pass filtered mechanical RPM for the speed PID and the telemetry frame.
 */

#ifndef RPM_FILTER_H_
#define RPM_FILTER_H_

#include <stdint.h>

void rpmFilterReset(int32_t rpm);
void rpmFilterUpdate(int32_t rpm);
int32_t getFilteredRPM(void);  // Filtered mechanical RPM, rounded
uint16_t getFilteredRPMForTelemetry(uint8_t pole_pairs);  // Filtered eRPM / 10, rounded, for the telemetry frame

#endif /* RPM_FILTER_H_ */
