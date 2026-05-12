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
#ifndef SPODY_SOLVER_H
#define SPODY_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Generic numeric solvers
 *
 * Today: a single root finder (Brent's method) used by the events
 * module to localise crossings inside a Runge-Kutta step (in
 * combination with the integrator dense output). Generic by design:
 * the function under test is a `double (*)(double, void*)`, the
 * args pointer is opaque -- so the solver is reusable for any
 * future scalar root problem (orbital geometry, constraint
 * residuals, fixed-point iterations, etc.).
 * ============================================================ */

#define SPODY_SOLVER_OK              0
#define SPODY_SOLVER_ERR_NULL       -1
#define SPODY_SOLVER_ERR_NOT_BRACKET -2  /* f(lo) and f(hi) have the same sign */
#define SPODY_SOLVER_ERR_MAX_ITER   -3

typedef double (*spody_scalar_fn)(double x, void *args);

/* Brent's method: combines bisection (always-converging) with inverse
 * quadratic interpolation (super-linear convergence in regular zones).
 * Standard reference: Brent, "Algorithms for minimization without
 * derivatives", 1973, ch. 4.
 *
 *   f         : function whose root is sought
 *   args      : opaque payload forwarded to f
 *   x_lo,x_hi : initial bracket. Must satisfy f(x_lo) * f(x_hi) <= 0.
 *   f_lo,f_hi : f evaluated at the endpoints (caller-provided to avoid
 *               re-evaluation when the caller already has them; pass any
 *               value if you don't, set use_provided to 0).
 *   use_provided_endpoints : 1 = use f_lo/f_hi as f(x_lo)/f(x_hi) directly
 *                            0 = the function will compute them
 *   tol       : absolute tolerance on the root (in the units of x).
 *   max_iter  : iteration cap (typical 50-100).
 *   x_root_out: output, the root.
 *
 * Returns SPODY_SOLVER_OK on convergence. */
int spody_solver_brent(spody_scalar_fn f, void *args,
                       double x_lo, double x_hi,
                       double f_lo, double f_hi,
                       int use_provided_endpoints,
                       double tol, int max_iter,
                       double *x_root_out);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_SOLVER_H */
