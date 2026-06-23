/*
 * Copyright 2026 ValeEng
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/* Classical Keplerian element <-> Cartesian state conversions.
 *
 * Inputs are the six classical orbital elements referenced to the
 * inertial frame the caller intends the output to live in: angles in
 * radians, sma in km, mu in km^3/s^2 -- output (r, v) in km, km/s.
 * The caller is responsible for choosing the reference frame (and,
 * for CR3BP, for translating/rotating the resulting state from the
 * primary-centered inertial frame into the synodic frame). */
#ifndef SPODY_KEPLER_H
#define SPODY_KEPLER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Solve Kepler's equation M = E - e*sin(E) for the eccentric anomaly
 * E. Newton iteration, machine-precision convergence in ~5 steps for
 * the elliptical regime (e < 1); capped at 30 as a hard safety. The
 * mean anomaly is wrapped to [-pi, pi] internally. */
double spody_kepler_solve_E(double mean_anom_rad, double ecc);

/* Mean -> true anomaly (radians). Uses spody_kepler_solve_E + the
 * standard half-angle formula. Returns true anomaly in (-pi, pi]. */
double spody_kepler_mean_to_true_anom(double mean_anom_rad, double ecc);

/* Classical elements -> Cartesian state in the same inertial frame
 * the elements reference.
 *
 *   sma_km           : semi-major axis (km), > 0 for elliptical
 *   ecc              : eccentricity, [0, 1) for elliptical
 *   inc_rad          : inclination, [0, pi]
 *   raan_rad         : right ascension of the ascending node, [0, 2pi)
 *   argp_rad         : argument of periapsis, [0, 2pi)
 *   true_anom_rad    : true anomaly nu, any value (wrapped internally)
 *   mu_km3_s2        : gravitational parameter of the reference body
 *   r_out, v_out     : output position (km) and velocity (km/s)
 *
 * Convention:
 *   inc = 0 means orbit lies in the reference frame's xy plane;
 *   raan = 0 means the ascending node coincides with the +x axis;
 *   argp = 0 means periapsis at the ascending node;
 *   true_anom = 0 means object at periapsis. */
void spody_keplerian_to_cartesian(double sma_km, double ecc,
                                  double inc_rad, double raan_rad,
                                  double argp_rad, double true_anom_rad,
                                  double mu_km3_s2,
                                  double r_out[3], double v_out[3]);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_KEPLER_H */
