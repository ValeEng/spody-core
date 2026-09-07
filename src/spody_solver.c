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
#include <float.h>
#include <math.h>
#include "spody_solver.h"

/* Brent's zeroin (Brent 1973, ch. 4; the same formulation Numerical
 * Recipes reprints as zbrent). Three points are carried: b is the best
 * estimate of the root so far, c is the previous best and always lies
 * on the other side of the root from b, so [b, c] is the bracket; a is
 * the point before that. Every iteration tries an inverse quadratic
 * interpolation through (a, b, c) -- the secant through (a, b) when two
 * of the three coincide -- and bisects [b, c] instead unless the
 * interpolated step is well inside the bracket and shrinks faster than
 * the one before last.
 *
 * Two details make it terminate in a handful of evaluations on a
 * smooth residual, and neither is optional:
 *
 *   - convergence is declared on the half-width of [b, c], the live
 *     bracket, not on the width of the original [x_lo, x_hi]. A far
 *     end that never moves -- the normal case when the iterates
 *     approach the root from one side -- does not have to be bisected
 *     down 40 times after b has already found the root;
 *
 *   - a step shorter than tol1 is stretched to tol1 in the direction
 *     of c. Once b sits on the root the next probe lands just across
 *     it, the sign flips, c collapses onto the old b and the bracket
 *     closes in one evaluation instead of a bisection cascade.
 *
 * tol1 = 2*DBL_EPSILON*|b| + tol/2: the caller's tol is an absolute
 * tolerance on x with a machine-precision floor, and the returned b
 * lies within tol1 of a sign change. */
int spody_solver_brent(spody_scalar_fn f, void *args,
                       double x_lo, double x_hi,
                       double f_lo, double f_hi,
                       int use_provided_endpoints,
                       double tol, int max_iter,
                       double *x_root_out)
{
    if (!f || !x_root_out) return SPODY_SOLVER_ERR_NULL;

    double a  = x_lo;
    double b  = x_hi;
    double fa = use_provided_endpoints ? f_lo : f(a, args);
    double fb = use_provided_endpoints ? f_hi : f(b, args);

    if (fa == 0.0) { *x_root_out = a; return SPODY_SOLVER_OK; }
    if (fb == 0.0) { *x_root_out = b; return SPODY_SOLVER_OK; }
    if ((fa > 0.0) == (fb > 0.0)) return SPODY_SOLVER_ERR_NOT_BRACKET;

    double c  = a;
    double fc = fa;
    double d  = b - a;      /* last step taken                          */
    double e  = d;          /* the step before that (controls fallback) */

    for (int iter = 0; iter < max_iter; iter++) {

        /* Keep c on the far side of the root from b. */
        if ((fb > 0.0) == (fc > 0.0)) {
            c = a; fc = fa;
            d = b - a; e = d;
        }
        /* Keep b the best estimate: |f(b)| <= |f(c)|. */
        if (fabs(fc) < fabs(fb)) {
            a = b;   b = c;   c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol1 = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * tol;
        double m    = 0.5 * (c - b);
        if (fabs(m) <= tol1 || fb == 0.0) {
            *x_root_out = b;
            return SPODY_SOLVER_OK;
        }

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            double s = fb / fa;
            double p, q;
            if (a == c) {
                /* two distinct points: secant */
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                /* three distinct points: inverse quadratic interpolation */
                double qq = fa / fc;
                double r  = fb / fc;
                p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;

            /* Take the interpolated step only if it stays inside the
             * bracket with room to spare and is at most half the step
             * before last; otherwise bisect. */
            if (2.0 * p < fmin(3.0 * m * q - fabs(tol1 * q), fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = m; e = m;
            }
        } else {
            d = m; e = m;
        }

        a = b; fa = fb;
        if (fabs(d) > tol1) b += d;
        else                b += (m > 0.0) ? tol1 : -tol1;
        fb = f(b, args);
    }

    *x_root_out = b;
    return SPODY_SOLVER_ERR_MAX_ITER;
}
