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

/* Bracketing binary search over an ascending array xs[0..n-1]:
 * returns i such that xs[i] <= x < xs[i+1], clamped to 0 when
 * x < xs[0] and to n-2 when x >= xs[n-1], so i+1 is always a valid
 * index. n must be >= 2. Shared primitive for every tabulated
 * lookup (density-scale nodes, resampled trajectories, ...); the
 * EOP / space weather parsers keep their own struct-array searches
 * because their records interleave several channels per row. */
size_t spody_bracket_index(const double *xs, size_t n, double x);

/* Piecewise-linear interpolation over ascending nodes (xs, ys),
 * clamped at both ends: x <= xs[0] returns ys[0], x >= xs[n-1]
 * returns ys[n-1]. n >= 1; n == 1 returns ys[0] everywhere. */
double spody_interp_linear(const double *xs, const double *ys,
                           size_t n, double x);

#ifdef __cplusplus
}
#endif

#endif // SPODY_MATH_H