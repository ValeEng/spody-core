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
#include <math.h>
#include "spody_solver.h"

/* Brent's root-finding method, classical formulation. The implementation
 * follows the Numerical Recipes / Wikipedia presentation: maintain three
 * approximations (a, b, c) with a current bracket [a, b], and at each
 * step decide between:
 *   - inverse quadratic interpolation (3 distinct values)
 *   - secant method (2 distinct values)
 *   - bisection fallback (whenever the trial step is unsafe)
 * The fallback set of conditions guarantees that |bracket| at least halves
 * every two iterations, so the method is as robust as bisection but
 * usually 5-10 iterations to ~1e-12 in regular regions. */
int spody_solver_brent(spody_scalar_fn f, void *args,
                       double x_lo, double x_hi,
                       double f_lo, double f_hi,
                       int use_provided_endpoints,
                       double tol, int max_iter,
                       double *x_root_out)
{
    if (!f || !x_root_out) return SPODY_SOLVER_ERR_NULL;

    double a = x_lo;
    double b = x_hi;
    double fa = use_provided_endpoints ? f_lo : f(a, args);
    double fb = use_provided_endpoints ? f_hi : f(b, args);

    if (fa == 0.0) { *x_root_out = a; return SPODY_SOLVER_OK; }
    if (fb == 0.0) { *x_root_out = b; return SPODY_SOLVER_OK; }
    if ((fa > 0.0) == (fb > 0.0)) return SPODY_SOLVER_ERR_NOT_BRACKET;

    /* keep |f(a)| >= |f(b)| so b is the best estimate so far */
    if (fabs(fa) < fabs(fb)) {
        double tmp = a; a = b; b = tmp;
        tmp = fa; fa = fb; fb = tmp;
    }

    double c  = a;
    double fc = fa;
    int    mflag = 1;          /* whether the previous step was bisection */
    double d = 0.0;            /* "previous-previous" iterate (used by mflag logic) */

    for (int iter = 0; iter < max_iter; iter++) {

        if (fabs(b - a) < tol) {
            *x_root_out = b;
            return SPODY_SOLVER_OK;
        }

        double s;
        if (fa != fc && fb != fc) {
            /* inverse quadratic interpolation */
            double L1 = (a * fb * fc) / ((fa - fb) * (fa - fc));
            double L2 = (b * fa * fc) / ((fb - fa) * (fb - fc));
            double L3 = (c * fa * fb) / ((fc - fa) * (fc - fb));
            s = L1 + L2 + L3;
        } else {
            /* secant */
            s = b - fb * (b - a) / (fb - fa);
        }

        /* Brent's safety conditions: fall back to bisection if any holds */
        double s_lo = (3.0 * a + b) * 0.25;
        double cond1 = (s < s_lo) ^ (s < b);                        /* not in (s_lo, b) */
        double cond2 = mflag      && fabs(s - b) >= fabs(b - c) * 0.5;
        double cond3 = !mflag     && fabs(s - b) >= fabs(c - d) * 0.5;
        double cond4 = mflag      && fabs(b - c) < tol;
        double cond5 = !mflag     && fabs(c - d) < tol;
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = 0.5 * (a + b);
            mflag = 1;
        } else {
            mflag = 0;
        }

        double fs = f(s, args);
        d  = c;
        c  = b;
        fc = fb;

        if ((fa > 0.0) != (fs > 0.0)) {
            b  = s;
            fb = fs;
        } else {
            a  = s;
            fa = fs;
        }
        if (fabs(fa) < fabs(fb)) {
            double tmp = a; a = b; b = tmp;
            tmp = fa; fa = fb; fb = tmp;
        }

        if (fs == 0.0) {
            *x_root_out = s;
            return SPODY_SOLVER_OK;
        }
    }

    *x_root_out = b;
    return SPODY_SOLVER_ERR_MAX_ITER;
}
