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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spody_integrators.h"
#include "spody_interp.h"   /* spody_hermite_cubic_1d for dense output */

//----- DP45 well-known constants (hardcoded, do not expose in options) ------

#define RKDP45_FACTOR_UPSCALE     5.0   // max growth of h after an accepted step
#define RKDP45_FACTOR_DOWNSCALE   0.1   // min shrinkage of h after a rejected step
#define RKDP45_RETRIES_PER_STEP   500   // safety cap on rejection retries per step

// Dormand-Prince 5(4) Butcher tableau -- "7S" pair from the same paper.
// J. R. Dormand and P. J. Prince, "A family of embedded Runge-Kutta formulae",
// J. Comput. Appl. Math., Vol. 6, 1980.
// Same order as the classical "7M" pair but with a different stage pattern
// that yields a more reliable embedded error estimate; this is the variant
// used by the LRO reference simulation, and the one we standardise on.
//   b_5 = a[6]  (FSAL: first-same-as-last)
//   b_4 = {431/5000, 0, 333/500, -7857/10000, 957/1000, 193/2000, -1/50}
//   e   = b_5 - b_4
static const double rkdp45_a[7][7] = {
    {0.0,         0.0,         0.0,         0.0,         0.0,        0.0,       0.0},
    {2.0/9.0,     0.0,         0.0,         0.0,         0.0,        0.0,       0.0},
    {1.0/12.0,    1.0/4.0,     0.0,         0.0,         0.0,        0.0,       0.0},
    {55.0/324.0, -25.0/108.0,  50.0/81.0,   0.0,         0.0,        0.0,       0.0},
    {83.0/330.0, -13.0/22.0,   61.0/66.0,   9.0/110.0,   0.0,        0.0,       0.0},
    {-19.0/28.0,  9.0/4.0,     1.0/7.0,    -27.0/7.0,    22.0/7.0,   0.0,       0.0},
    {19.0/200.0,  0.0,         3.0/5.0,    -243.0/400.0, 33.0/40.0,  7.0/80.0,  0.0}
};

static const double rkdp45_c[7] = {
    0.0,
    2.0/9.0,
    1.0/3.0,
    5.0/9.0,
    2.0/3.0,
    1.0,
    1.0
};

// error coefficients e = b_5 - b_4
static const double rkdp45_e[7] = {
    (19.0/200.0   - 431.0/5000.0),
     0.0,
    (3.0/5.0      - 333.0/500.0),
    (-243.0/400.0 + 7857.0/10000.0),
    (33.0/40.0    - 957.0/1000.0),
    (7.0/80.0     - 193.0/2000.0),
     1.0/50.0
};

// Classical Runge-Kutta 4th order Butcher tableau
static const double rk4_a[4][4] = {
    {0.0, 0.0, 0.0, 0.0},
    {0.5, 0.0, 0.0, 0.0},
    {0.0, 0.5, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0}
};

static const double rk4_c[4] = {
    0.0,
    0.5,
    0.5,
    1.0
};

static const double rk4_b[4] = {
    1.0/6.0,
    1.0/3.0,
    1.0/3.0,
    1.0/6.0
};

//----- internal helpers --------------------------------------------------

static int alloc_buf(double **p, size_t n) {
    *p = (double*)malloc(n * sizeof(double));
    return (*p) ? 0 : SPODY_INTEG_ERR_NULL;
}

static void free_buf(double **p) {
    if (*p) { free(*p); *p = NULL; }
}

static void zero_all_buffers(IntegratorAllData *integ) {
    integ->k     = NULL;
    integ->y_tmp = NULL;
    integ->y_err = NULL;
    integ->y     = NULL;
    integ->y_old = NULL;
}

//----- defaults ----------------------------------------------------------

void spody_default_integrator_options(spody_integrator_method method, IntegratorOptions *opt) {
    if (!opt) return;
    opt->h_init    = 10.0;
    opt->h_min     = 1e-6;
    opt->h_max     = 600.0;
    opt->abs_tol   = 1e-9;
    opt->rel_tol   = 1e-9;
    opt->max_steps = 0;     // unlimited
    opt->safety    = 0.9;

    switch (method) {
        case SPODY_INTEG_RK4:
            opt->h_init = 10.0;
            break;
        case SPODY_INTEG_RK45:
            opt->h_init = 10.0;
            opt->rel_tol = 1e-9;
            break;
        case SPODY_INTEG_RK78:
            opt->h_init = 30.0;
            opt->rel_tol = 1e-12;
            break;
        case SPODY_INTEG_VERLET:
            opt->h_init = 1.0;
            break;
    }
}

//----- setup / teardown --------------------------------------------------

int spody_setup_integrator(IntegratorAllData *integ,
                           spody_integrator_method method,
                           const IntegratorOptions *opt,
                           int dim,
                           spody_rhs_fn rhs,
                           void *user) {
    if (!integ || !rhs) return SPODY_INTEG_ERR_NULL;
    if (dim <= 0)       return SPODY_INTEG_ERR_DIM;

    zero_all_buffers(integ);

    integ->method   = method;
    integ->dim      = dim;
    integ->rhs      = rhs;
    integ->user     = user;
    integ->pre_step = NULL;   // opt-in; see spody_prestep_fn
    integ->t      = 0.0;
    integ->h_old  = 0.0;
    integ->t_old  = 0.0;

    integ->n_accepted = 0;
    integ->n_rejected = 0;
    integ->n_rhs      = 0;

    if (opt) {
        integ->opt = *opt;
    } else {
        spody_default_integrator_options(method, &integ->opt);
    }
    integ->h = integ->opt.h_init;

    // common allocations
    if (alloc_buf(&integ->y, (size_t)dim))     goto fail;
    if (alloc_buf(&integ->y_tmp, (size_t)dim)) goto fail;
    if (alloc_buf(&integ->y_old, (size_t)dim)) goto fail;   // for dense output

    // method-specific scratch sizing for the stage buffer `k`
    int n_stages = 0;
    int needs_yerr = 0;
    switch (method) {
        case SPODY_INTEG_RK4:    n_stages = 4;  needs_yerr = 0; break;
        case SPODY_INTEG_RK45:   n_stages = 7;  needs_yerr = 1; break;
        case SPODY_INTEG_RK78:   n_stages = 13; needs_yerr = 1; break;
        case SPODY_INTEG_VERLET:
            // Verlet expects state laid out as [r(0..dim/2-1), v(dim/2..dim-1)].
            if ((dim & 1) != 0) goto fail_dim;
            n_stages = 1; needs_yerr = 0;
            break;
        default: goto fail;
    }

    if (alloc_buf(&integ->k, (size_t)(n_stages * dim))) goto fail;
    if (needs_yerr) {
        if (alloc_buf(&integ->y_err, (size_t)dim)) goto fail;
    }

    return SPODY_INTEG_OK;

fail_dim:
    spody_free_integrator(integ);
    return SPODY_INTEG_ERR_DIM;

fail:
    spody_free_integrator(integ);
    return SPODY_INTEG_ERR_NULL;
}

int spody_free_integrator(IntegratorAllData *integ) {
    if (!integ) return SPODY_INTEG_ERR_NULL;
    free_buf(&integ->y);
    free_buf(&integ->y_tmp);
    free_buf(&integ->y_err);
    free_buf(&integ->k);
    free_buf(&integ->y_old);
    return SPODY_INTEG_OK;
}

int spody_set_integrator_state(IntegratorAllData *integ, double t0, const double *y0) {
    if (!integ || !y0) return SPODY_INTEG_ERR_NULL;
    if (!integ->y)     return SPODY_INTEG_ERR_NULL;
    integ->t = t0;
    memcpy(integ->y, y0, (size_t)integ->dim * sizeof(double));
    return SPODY_INTEG_OK;
}

//----- RKDP45 (Dormand-Prince 5(4), adaptive, GMAT-style step control) -----

static int step_rkdp45(IntegratorAllData *integ) {

    //A family of embedded Runge-Kutta formulae
    //J. R. Dormand and P. J. Prince

    const int dim   = integ->dim;
    double *state   = integ->y;
    double *temp    = integ->y_tmp;     // intermediate state for stage evaluation
    double *k       = integ->k;          // flat buffer: k[stage*dim + i]
    double *rel_err = integ->y_err;      // accumulator: e[j] * k[j]

    double clock = integ->t;
    double temp_clock = clock;
    int returnNumber;
    int steps = 0;

    // Snapshot of the state at the start of the step. Needed for dense
    // output (y(theta) = y_old + h_old * sum_i b_i(theta) * k_i).
    memcpy(integ->y_old, state, sizeof(double) * (size_t)dim);

    do {

        for (int j = 0; j < 7; j++) {

            memcpy(temp, state, sizeof(double) * (size_t)dim);

            #if DEBUG_INTEGRATORS == 1
            printf("[rkdp45_01][K%d] K%d evaluation\n", j+1, j+1);
            #endif

            for (int i = 0; i < dim; i++) {
                for (int x = 0; x < j; x++) {

                    temp[i] += rkdp45_a[j][x] * k[x*dim + i];

                    #if DEBUG_INTEGRATORS == 1
                    printf("[rkdp45_03][K%d] vector position %d | column %d :  a[%d][%d]--->temp %.21f\n", j+1, i, x, j, x, temp[i]);
                    #endif
                }
                #if DEBUG_INTEGRATORS == 1
                printf("[rkdp45_02][K%d] vector position %d | temp tot:  %.21f\n", j+1, i, temp[i]);
                #endif
            }

            temp_clock = clock + integ->h * rkdp45_c[j];
            returnNumber = integ->rhs(temp_clock, temp, k + j*dim, integ->user);
            integ->n_rhs++;
            if (returnNumber != 0) return SPODY_INTEG_ERR_RHS;

            for (int i = 0; i < dim; i++) {
                k[j*dim + i] *= integ->h;   // save the K factor already multiplied by h
            }

            #if DEBUG_INTEGRATORS == 1
            for (int i = 0; i < dim; i++) {
                printf("[rkdp45_04][K%d] vector position %d : %.21f | sorting \n", j+1, i, k[j*dim + i]);
            }
            #endif

        }

        // accumulate the embedded-pair error estimate: rel_err[i] = sum_j e[j] * k[j][i]
        for (int i = 0; i < dim; i++) {
            rel_err[i] = 0.0;
            for (int j = 0; j < 7; j++) {
                rel_err[i] += rkdp45_e[j] * k[j*dim + i];

                #if DEBUG_INTEGRATORS == 1
                printf("[rkdp45_05] vector position %d | e[%d] factor | relative error : %.21f \n", i, j, rel_err[i]);
                #endif
            }
        }

        #if DEBUG_INTEGRATORS == 1
        printf("[rkdp45_06] relative error : ");
        for (int i = 0; i < dim; i++) printf("%.21f | ", rel_err[i]);
        printf("\n");
        #endif

        // ---> error estimation RSS Step from GMAT | Thank you!
        // Iterate over the state in 3-component sub-vectors and take the worst block.
        double err = 0.0;

        for (int i = 0; i + 2 < dim; i += 3) {

            double dx = temp[i]   - state[i];
            double dy = temp[i+1] - state[i+1];
            double dz = temp[i+2] - state[i+2];

            double mag_delta_state = dx*dx + dy*dy + dz*dz;
            double err_i = rel_err[i]*rel_err[i] + rel_err[i+1]*rel_err[i+1] + rel_err[i+2]*rel_err[i+2];

            #if DEBUG_INTEGRATORS == 1
            printf("[rkdp45_07] i = %d | delta state %.21f , %.21f , %.21f | mag delta %.21f | mag rel err %.21f \n", i, dx, dy, dz, mag_delta_state, err_i);
            #endif

            if (mag_delta_state > 0.1) {
                err_i = sqrt(err_i / mag_delta_state);
            } else {
                err_i = sqrt(err_i);
            }
            if (err_i > err) err = err_i;

            #if DEBUG_INTEGRATORS == 1
            printf("[rkdp45_08] sqrt(mag rel err/mag delta) %.21f | if equal to (mag rel err) no division occurred\n", err);
            #endif
        }

        steps++;

        #if DEBUG_INTEGRATORS == 1
        printf("[rkdp45_08] Iteration : %d | Estimated Error : %.2e | Tolerance : %.2e \n", steps, err, integ->opt.rel_tol);
        #endif

        if (err < integ->opt.rel_tol) {

            // from GMAT
            double incPower = 1.0 / 5.0;
            double scale = integ->opt.safety * pow(integ->opt.rel_tol / err, incPower); // safety factor = 0.9 (default)
            if (scale > RKDP45_FACTOR_UPSCALE) scale = RKDP45_FACTOR_UPSCALE;

            memcpy(state, temp, sizeof(double) * (size_t)dim);
            integ->h_old = integ->h;
            integ->t_old = clock;
            integ->t     = clock + integ->h_old;

            integ->h *= fmin(RKDP45_FACTOR_UPSCALE, scale);
            if (integ->h > integ->opt.h_max) integ->h = integ->opt.h_max;

            integ->n_accepted++;

            #if DEBUG_INTEGRATORS == 1
            printf("[RKDP45_09][YES] iter : %d | start clock : %.6f | old time step : %.6f | new time step : %.6f | err : %.2e | scale : %.6f\n\n", steps, clock, integ->h_old, integ->h, err, scale);
            #endif

            return SPODY_INTEG_OK;

        } else {

            integ->n_rejected++;

            // from GMAT
            double decPower = 1.0 / 4.0;
            double scale = integ->opt.safety * pow(integ->opt.rel_tol / err, decPower); // safety factor = 0.9 (default)
            if (scale < RKDP45_FACTOR_DOWNSCALE) scale = RKDP45_FACTOR_DOWNSCALE;

            double h_prev = integ->h;
            integ->h *= fmax(RKDP45_FACTOR_DOWNSCALE, scale);
            if (fabs(integ->h) < integ->opt.h_min) {
                #if DEBUG_INTEGRATORS == 1
                printf("[RKDP45_09][STEP_TOO_SMALL] reached h_min : %.6e\n", integ->opt.h_min);
                #endif
                return SPODY_INTEG_ERR_STEP_TOO_SMALL;
            }

            #if DEBUG_INTEGRATORS == 1
            printf("[RKDP45_09][NO] iter : %d | start clock : %.6f | old time step : %.6f | new time step : %.6f | err : %.2e | scale : %.6f\n\n", steps, clock, h_prev, integ->h, err, scale);
            #else
            (void)h_prev;
            #endif
        }

    } while (steps < RKDP45_RETRIES_PER_STEP);

    printf("[RKDP45_05][MAX_ITER] !Iteration limit reached! ---> Bad constraints\n");
    return SPODY_INTEG_ERR_MAX_STEPS;
}

//----- RK4 (classical Runge-Kutta 4th order, fixed step) ----------------

static int step_rk4(IntegratorAllData *integ) {

    const int dim = integ->dim;
    double *state = integ->y;
    double *temp  = integ->y_tmp;
    double *k     = integ->k;          // flat buffer: k[stage*dim + i]

    double clock = integ->t;
    double temp_clock = clock;
    int returnNumber;

    // snapshot for dense output bookkeeping (RK4 dense_eval not yet implemented)
    memcpy(integ->y_old, state, sizeof(double) * (size_t)dim);

    for (int j = 0; j < 4; j++) {

        memcpy(temp, state, sizeof(double) * (size_t)dim);

        #if DEBUG_INTEGRATORS == 1
        printf("[rk4_01][K%d] K%d evaluation\n", j+1, j+1);
        #endif

        for (int i = 0; i < dim; i++) {
            for (int x = 0; x < j; x++) {

                temp[i] += rk4_a[j][x] * k[x*dim + i];

                #if DEBUG_INTEGRATORS == 1
                printf("[rk4_03][K%d] vector position %d | column %d :  a[%d][%d]--->temp %.21f\n", j+1, i, x, j, x, temp[i]);
                #endif
            }
            #if DEBUG_INTEGRATORS == 1
            printf("[rk4_02][K%d] vector position %d | temp tot:  %.21f\n", j+1, i, temp[i]);
            #endif
        }

        temp_clock = clock + integ->h * rk4_c[j];
        returnNumber = integ->rhs(temp_clock, temp, k + j*dim, integ->user);
        integ->n_rhs++;
        if (returnNumber != 0) return SPODY_INTEG_ERR_RHS;

        for (int i = 0; i < dim; i++) {
            k[j*dim + i] *= integ->h;   // save the K factor already multiplied by h
        }

        #if DEBUG_INTEGRATORS == 1
        for (int i = 0; i < dim; i++) {
            printf("[rk4_04][K%d] vector position %d : %.21f | sorting \n", j+1, i, k[j*dim + i]);
        }
        #endif
    }

    // y_{n+1} = y_n + sum_j b[j] * k[j]
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 4; j++) {
            state[i] += rk4_b[j] * k[j*dim + i];
        }
    }

    integ->h_old = integ->h;
    integ->t = clock + integ->h_old;

    /* Fixed step: every step is accepted by construction, so
     * n_rejected stays at zero for this method. */
    integ->n_accepted++;

    #if DEBUG_INTEGRATORS == 1
    printf("[RK4_05] start clock : %.6f | time step : %.6f | new clock : %.6f\n\n", clock, integ->h_old, integ->t);
    #endif

    return SPODY_INTEG_OK;
}

//----- other methods (placeholders) -------------------------------------
// To be implemented following the same style of step_rkdp45 / step_rk4.

static int step_rk78(IntegratorAllData *integ)   { (void)integ; return SPODY_INTEG_ERR_NULL; }
static int step_verlet(IntegratorAllData *integ) { (void)integ; return SPODY_INTEG_ERR_NULL; }

//----- public step / drive ----------------------------------------------

int spody_propagate_onestep(IntegratorAllData *integ) {
    if (!integ || !integ->rhs || !integ->y) return SPODY_INTEG_ERR_NULL;
    /* One call per step, before any stage: the hook's whole purpose is
     * that whatever it tunes stays FIXED for every stage below,
     * retries included. Placed here rather than inside each step_*()
     * so all methods get it from one site. */
    if (integ->pre_step
            && integ->pre_step(integ->t, integ->y, integ->h,
                               integ->user) != 0) {
        return SPODY_INTEG_ERR_RHS;
    }
    switch (integ->method) {
        case SPODY_INTEG_RK4:    return step_rk4(integ);
        case SPODY_INTEG_RK45:   return step_rkdp45(integ);
        case SPODY_INTEG_RK78:   return step_rk78(integ);
        case SPODY_INTEG_VERLET: return step_verlet(integ);
    }
    return SPODY_INTEG_ERR_NULL;
}

int spody_propagate_untilend(IntegratorAllData *integ, double t_end) {
    if (!integ) return SPODY_INTEG_ERR_NULL;

    size_t n_steps = 0;
    const size_t cap = integ->opt.max_steps;

    while (1) {
        double remaining = t_end - integ->t;
        if (remaining == 0.0) return SPODY_INTEG_OK;

        // align step direction with target
        if ((remaining > 0.0 && integ->h <= 0.0) ||
            (remaining < 0.0 && integ->h >= 0.0)) {
            integ->h = (remaining > 0.0 ? integ->opt.h_init : -integ->opt.h_init);
        }

        // clip last step so we land exactly on t_end
        if (fabs(integ->h) > fabs(remaining)) integ->h = remaining;

        int returnNumber = spody_propagate_onestep(integ);
        if (returnNumber != SPODY_INTEG_OK) return returnNumber;

        n_steps++;
        if (cap && n_steps >= cap) return SPODY_INTEG_ERR_MAX_STEPS;
    }
}

/* ============================================================
 * Dense output: cubic Hermite C^1 via FSAL endpoint derivatives
 *
 * For RKDP45 (Dormand-Prince 5(4)) the integrator advance here uses the
 * "7S" Butcher tableau (see top of file -- chosen because its embedded
 * error estimate is the more reliable of the two DP5 variants). The
 * classical 5th-order DOPRI5 dense-output coefficients from Hairer-
 * Wanner Table 6.1 are tied to the "7M" k-stage pattern instead, so
 * applying that P-matrix to the 7S stages produces an interpolant that
 * matches at the endpoints (by construction of any reasonable formula)
 * but drifts off the integrated trajectory mid-step. The bug is
 * harmless for trajectory plotting (errors stay small) but it breaks
 * downstream localisation -- spody_event_check_refined runs Brent on
 * the dense curve, so a non-integrator-consistent interpolant gives a
 * wrong t_trigger.
 *
 * Cubic Hermite sidesteps the issue entirely: it depends only on the
 * endpoint values and the endpoint derivatives, both of which are
 * integrator-consistent for any Butcher tableau. The accuracy drops
 * one order vs the DOPRI5 standard formula (4th order on y vs 5th),
 * still strictly higher than the method's local truncation error;
 * over a 30 s LRO step this is sub-microsecond on a time-localised
 * event surface, far below any physically meaningful threshold.
 *
 * FSAL note: for RKDP45 the last stage is "first same as last", so
 *     k_1 = h_old * f(t_old, y_old)
 *     k_7 = h_old * f(t_new, y_new)
 * are exactly the endpoint derivatives (times h_old) we need. No
 * extra RHS evaluation is performed.
 *
 * The Hermite basis itself lives in spody_interp.{h,c}; we only stage
 * the inputs here.
 * ============================================================ */
int spody_dense_eval(const IntegratorAllData *integ, double theta, double *y_out) {
    if (!integ || !y_out) return SPODY_INTEG_ERR_NULL;
    if (integ->method != SPODY_INTEG_RK45) return SPODY_INTEG_ERR_NULL;
    if (!integ->y_old || !integ->y || !integ->k) return SPODY_INTEG_ERR_NULL;

    if (theta < 0.0) theta = 0.0;
    if (theta > 1.0) theta = 1.0;

    /* spody_hermite_cubic_1d expects derivatives in physical units
     * (per-second), but k-stages are stored pre-multiplied by h_old.
     * Convert once and pass the scalars to the library helper component
     * by component -- the basis computation is cheap and the call
     * pattern stays dim-agnostic without needing a temporary buffer. */
    const int    dim   = integ->dim;
    const double h_old = integ->h_old;
    const double t_q   = integ->t_old + theta * h_old;
    const double inv_h = (h_old != 0.0) ? 1.0 / h_old : 0.0;
    const double *k1   = &integ->k[0 * dim];   /* h_old * f(t_old, y_old)  */
    const double *k7   = &integ->k[6 * dim];   /* h_old * f(t_new, y_new)  */

    for (int i = 0; i < dim; i++) {
        y_out[i] = spody_hermite_cubic_1d(
                        t_q,
                        integ->t_old, integ->t,
                        integ->y_old[i], k1[i] * inv_h,
                        integ->y[i],     k7[i] * inv_h);
    }
    return SPODY_INTEG_OK;
}
