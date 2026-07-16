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
/*
 * Interpolation primitives.
 *
 * Pure-math routines, no dependency on integrator state or any other
 * spody-core module. THE home of everything that evaluates tabulated
 * data: the bracketing lookup, piecewise-linear nodes, cubic Hermite
 * dense output. Future additions (cubic spline, Lagrange, B-spline --
 * e.g. for an SPK Type 9/13 reader) should live in this same header.
 *
 * The Hermite functions clamp the parameter s = (t - t_a) / (t_b - t_a)
 * into [0, 1], so a query slightly outside [t_a, t_b] (e.g. due to
 * floating point rounding on a grid edge) returns the endpoint value
 * rather than extrapolating wildly; the tabulated lookups clamp to the
 * first/last node the same way.
 */
#ifndef SPODY_INTERP_H
#define SPODY_INTERP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

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

/* Scalar cubic Hermite. Returns the C^1 cubic at t that matches
 *   y(t_a) = y_a,   y'(t_a) = dy_a,
 *   y(t_b) = y_b,   y'(t_b) = dy_b.
 * Caller must ensure t_a != t_b. */
double spody_hermite_cubic_1d(double t,
                              double t_a, double t_b,
                              double y_a, double dy_a,
                              double y_b, double dy_b);

/* Vector cubic Hermite: applies the same scalar interpolant component-
 * wise to a state of `dim` doubles. Optionally also returns the
 * derivative at t (pass dy_out = NULL to skip). The y_out and dy_out
 * buffers must hold at least `dim` doubles each. */
void spody_hermite_cubic(double t,
                         double t_a, double t_b,
                         const double *y_a, const double *dy_a,
                         const double *y_b, const double *dy_b,
                         int dim,
                         double *y_out, double *dy_out);

/* ODE dense-output convenience for the standard (r, v) state layout:
 *   y      = [r(3), v(3)]   with v = dr/dt
 *   y_a    = state at t_a
 *   y_b    = state at t_b
 *   y_out  = interpolated state at t, same layout
 *
 * Internally calls spody_hermite_cubic on the position triple using
 * the velocity portion of y as the derivative; the velocity portion
 * of y_out is the analytic derivative of that cubic, so position and
 * velocity are C^1-consistent at every t. */
void spody_hermite_dense_rv6(double t,
                             double t_a, const double y_a[6],
                             double t_b, const double y_b[6],
                             double y_out[6]);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_INTERP_H */
