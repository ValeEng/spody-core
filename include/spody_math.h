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
#ifndef SPODY_MATH_H
#define SPODY_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdio.h>

void spody_transpose_matrix(double in[3][3], double out[3][3]);
void spody_getrotmatrix_XZX(double alpha, double beta, double gamma, double R[3][3]);
void spody_rotate_vector(const double C[3][3], const double v_in[3], double v_out[3]);

/* 3-vector dot and cross products. Shared primitives for frame /
 * basis construction (e.g. the RIC triad in trajectory-diff tooling);
 * out must not alias a or b in spody_cross3. */
double spody_dot3(const double a[3], const double b[3]);
void   spody_cross3(const double a[3], const double b[3], double out[3]);

/* Body-fixed cartesian -> geodetic coordinates on the oblate
 * ellipsoid (semi-major axis a_km, inverse flattening inv_f; pass
 * WGS84_A_KM / WGS84_INV_F from spody_const.h for Earth). Bowring's
 * method with a fixed parametric-latitude refinement: machine-
 * precision latitude/altitude for any Earth-orbit point (verified by
 * forward/inverse round-trip up to GEO heights). Ellipsoid-
 * parameterised on purpose: a future Mars central body reuses it
 * with its own a/f, nothing here is Earth-specific.
 * Outputs: geodetic latitude (rad), longitude (rad, atan2 range),
 * altitude above the ellipsoid (km). Any output may be NULL. */
void spody_bf_to_geodetic(const double r_bf_km[3], double a_km,
                            double inv_f, double *lat_rad,
                            double *lon_rad, double *alt_km);

/* Shape of a body for shadow, impact and altitude: a spheroid of
 * equatorial radius r_eq and polar radius r_pol [km] about the unit
 * spin axis `pole`, expressed in the frame of the vectors it is used
 * with (the integration frame, ICRF axes). r_pol <= 0 or r_pol == r_eq
 * is a sphere of radius r_eq, and every consumer then takes the
 * spherical path of old, bit for bit. r_eq <= 0: no shape. */
typedef struct {
    double r_eq;
    double r_pol;
    double pole[3];
} SpodyBodyShape;

/* Spin axis in ICRF from the IAU WGCCRE linear elements (pck00011
 * BODYnnn_POLE_RA / _POLE_DEC): RA = ra0 + ra1*T, DEC = dec0 + dec1*T,
 * degrees and degrees per Julian century, T counted from J2000 TDB
 * (et in seconds). */
void spody_iau_pole(double ra0_deg, double ra1_deg, double dec0_deg,
                    double dec1_deg, double et, double pole[3]);

/* "Distance" of a point from the body centre measured against the
 * body's shape, so that (result - r_eq) is the altitude above it.
 *   sphere   : |r_rel| exactly (the spherical code path of old);
 *   spheroid : r_eq + geodetic altitude above the spheroid (Bowring,
 *              spody_bf_to_geodetic, in the frame of the pole).
 * r_rel: satellite relative to the body centre, same frame as pole. */
double spody_body_shape_distance(const SpodyBodyShape *shape,
                                 const double r_rel[3]);

#ifdef __cplusplus
}
#endif

#endif // SPODY_MATH_H