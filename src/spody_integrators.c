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
#include "spody_interp.h"   /* Hermite interpolants for dense output */
#include "spody_io.h"

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

//----- DOP853 constants (hardcoded, do not expose in options) ---------------
// Step-size control with the defaults of dop853.f (E. Hairer, G. Wanner,
// version of October 11, 2009): h_new = h * fac,
// fac = safety * (tol/err)^(1/8) limited to [0.333, 6], no growth on the
// step after a rejection. E. Hairer, S. P. Norsett, G. Wanner, "Solving
// Ordinary Differential Equations I. Nonstiff Problems", 2nd ed.,
// Springer, 1993, Sect. II.4.

#define DOP853_STAGES             12    // stage evaluations per attempt (stage 1 is the FSAL derivative)
#define DOP853_FACTOR_UPSCALE     6.0   // max growth of h after an accepted step
#define DOP853_FACTOR_DOWNSCALE   0.333 // min shrinkage of h after a rejected step
#define DOP853_RETRIES_PER_STEP   500   // safety cap on rejection retries per step

// Dormand-Prince 8(5,3) tableau (DOP853): explicit Runge-Kutta method of
// order 8 with embedded error estimators of order 5 and 3 (Hairer,
// Norsett, Wanner 1993, Sect. II.10). Coefficients as in dop853.f.
//   k_i   = h f(t + c_i h, y + sum_j a_ij k_j)
//   y_new = y + sum_i b_i k_i
//   err5  = sum_i e5_i k_i
//   err3  = sum_i (b_i - bhh_i) k_i
// f(t + h, y_new) is the first stage of the next step (FSAL).
static const double dop853_c[DOP853_STAGES] = {
    0.0,
    0.526001519587677318785587544488e-01,
    0.789002279381515978178381316732e-01,
    0.118350341907227396726757197510,
    0.281649658092772603273242802490,
    0.333333333333333333333333333333,
    0.25,
    0.307692307692307692307692307692,
    0.651282051282051282051282051282,
    0.6,
    0.857142857142857142857142857142,
    1.0,
};

static const double dop853_a[DOP853_STAGES][DOP853_STAGES] = {
    {0.0},
    { /* stage 2 */
        [ 0] = 5.26001519587677318785587544488e-2,
    },
    { /* stage 3 */
        [ 0] = 1.97250569845378994544595329183e-2,
        [ 1] = 5.91751709536136983633785987549e-2,
    },
    { /* stage 4 */
        [ 0] = 2.95875854768068491816892993775e-2,
        [ 2] = 8.87627564304205475450678981324e-2,
    },
    { /* stage 5 */
        [ 0] = 2.41365134159266685502369798665e-1,
        [ 2] = -8.84549479328286085344864962717e-1,
        [ 3] = 9.24834003261792003115737966543e-1,
    },
    { /* stage 6 */
        [ 0] = 3.7037037037037037037037037037e-2,
        [ 3] = 1.70828608729473871279604482173e-1,
        [ 4] = 1.25467687566822425016691814123e-1,
    },
    { /* stage 7 */
        [ 0] = 3.7109375e-2,
        [ 3] = 1.70252211019544039314978060272e-1,
        [ 4] = 6.02165389804559606850219397283e-2,
        [ 5] = -1.7578125e-2,
    },
    { /* stage 8 */
        [ 0] = 3.70920001185047927108779319836e-2,
        [ 3] = 1.70383925712239993810214054705e-1,
        [ 4] = 1.07262030446373284651809199168e-1,
        [ 5] = -1.53194377486244017527936158236e-2,
        [ 6] = 8.27378916381402288758473766002e-3,
    },
    { /* stage 9 */
        [ 0] = 6.24110958716075717114429577812e-1,
        [ 3] = -3.36089262944694129406857109825,
        [ 4] = -8.68219346841726006818189891453e-1,
        [ 5] = 2.75920996994467083049415600797e1,
        [ 6] = 2.01540675504778934086186788979e1,
        [ 7] = -4.34898841810699588477366255144e1,
    },
    { /* stage 10 */
        [ 0] = 4.77662536438264365890433908527e-1,
        [ 3] = -2.48811461997166764192642586468,
        [ 4] = -5.90290826836842996371446475743e-1,
        [ 5] = 2.12300514481811942347288949897e1,
        [ 6] = 1.52792336328824235832596922938e1,
        [ 7] = -3.32882109689848629194453265587e1,
        [ 8] = -2.03312017085086261358222928593e-2,
    },
    { /* stage 11 */
        [ 0] = -9.3714243008598732571704021658e-1,
        [ 3] = 5.18637242884406370830023853209,
        [ 4] = 1.09143734899672957818500254654,
        [ 5] = -8.14978701074692612513997267357,
        [ 6] = -1.85200656599969598641566180701e1,
        [ 7] = 2.27394870993505042818970056734e1,
        [ 8] = 2.49360555267965238987089396762,
        [ 9] = -3.0467644718982195003823669022,
    },
    { /* stage 12 */
        [ 0] = 2.27331014751653820792359768449,
        [ 3] = -1.05344954667372501984066689879e1,
        [ 4] = -2.00087205822486249909675718444,
        [ 5] = -1.79589318631187989172765950534e1,
        [ 6] = 2.79488845294199600508499808837e1,
        [ 7] = -2.85899827713502369474065508674,
        [ 8] = -8.87285693353062954433549289258,
        [ 9] = 1.23605671757943030647266201528e1,
        [10] = 6.43392746015763530355970484046e-1,
    },
};

/* b: 8th-order weights (= A[12] in the extended tableau) */
static const double dop853_b[DOP853_STAGES] = {
    [ 0] = 5.42937341165687622380535766363e-2,
    [ 5] = 4.45031289275240888144113950566,
    [ 6] = 1.89151789931450038304281599044,
    [ 7] = -5.8012039600105847814672114227,
    [ 8] = 3.1116436695781989440891606237e-1,
    [ 9] = -1.52160949662516078556178806805e-1,
    [10] = 2.01365400804030348374776537501e-1,
    [11] = 4.47106157277725905176885569043e-2,
};

static const double dop853_e5[DOP853_STAGES] = {
    [ 0] = 0.1312004499419488073250102996e-1,
    [ 5] = -0.1225156446376204440720569753e+1,
    [ 6] = -0.4957589496572501915214079952,
    [ 7] = 0.1664377182454986536961530415e+1,
    [ 8] = -0.3503288487499736816886487290,
    [ 9] = 0.3341791187130174790297318841,
    [10] = 0.8192320648511571246570742613e-1,
    [11] = -0.2235530786388629525884427845e-1,
};

static const double dop853_bhh[DOP853_STAGES] = {
    [ 0] = 0.244094488188976377952755905512,
    [ 8] = 0.733846688281611857341361741547,
    [11] = 0.220588235294117647058823529412e-1,
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
    integ->f_now = NULL;
    integ->f_new = NULL;
}

//----- defaults ----------------------------------------------------------

void spody_default_integrator_options(spody_integrator_method method, IntegratorOptions *opt) {
    if (!opt) return;
    opt->h_init    = 10.0;
    opt->h_min     = 1e-6;
    opt->h_max     = 600.0;
    opt->rel_tol   = 1e-9;
    opt->max_steps = 0;     // unlimited
    opt->safety    = 0.9;

    switch (method) {
        case SPODY_INTEG_RK4:
            opt->h_init = 10.0;
            break;
        case SPODY_INTEG_RK45:
        case SPODY_INTEG_DOP853:
            opt->h_init = 10.0;
            opt->rel_tol = 1e-9;
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

    integ->method = method;
    integ->dim    = dim;
    integ->rhs    = rhs;
    integ->user   = user;
    integ->t      = 0.0;
    integ->h_old  = 0.0;
    integ->t_old  = 0.0;

    integ->n_accepted = 0;
    integ->n_rejected = 0;
    integ->n_rhs      = 0;
    integ->fsal_valid = 0;

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
        case SPODY_INTEG_DOP853: n_stages = DOP853_STAGES; needs_yerr = 0; break;
        default: goto fail;
    }

    if (alloc_buf(&integ->k, (size_t)(n_stages * dim))) goto fail;
    if (needs_yerr) {
        if (alloc_buf(&integ->y_err, (size_t)dim)) goto fail;
    }
    /* FSAL derivative buffers: the DP5(4) 7S tableau and DOP853 have
     * the property (the derivative at the step's end is the next
     * step's first stage). RK4 does not. */
    if (method == SPODY_INTEG_RK45 || method == SPODY_INTEG_DOP853) {
        if (alloc_buf(&integ->f_now, (size_t)dim)) goto fail;
        if (alloc_buf(&integ->f_new, (size_t)dim)) goto fail;
    }

    return SPODY_INTEG_OK;

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
    free_buf(&integ->f_now);
    free_buf(&integ->f_new);
    return SPODY_INTEG_OK;
}

int spody_set_integrator_state(IntegratorAllData *integ, double t0, const double *y0) {
    if (!integ || !y0) return SPODY_INTEG_ERR_NULL;
    if (!integ->y)     return SPODY_INTEG_ERR_NULL;
    integ->t = t0;
    memcpy(integ->y, y0, (size_t)integ->dim * sizeof(double));
    integ->fsal_valid = 0;
    return SPODY_INTEG_OK;
}

void spody_integrator_invalidate_fsal(IntegratorAllData *integ) {
    if (integ) integ->fsal_valid = 0;
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

            if (j == 0) {
                /* Stage 1 is f(t, y). It is already in f_now when the
                 * previous step's last stage produced it (FSAL) or
                 * when a rejected attempt evaluated it: the state has
                 * not moved, so it is evaluated at most once per step. */
                if (!integ->fsal_valid) {
                    returnNumber = integ->rhs(temp_clock, temp, integ->f_now, integ->user);
                    integ->n_rhs++;
                    if (returnNumber != 0) return SPODY_INTEG_ERR_RHS;
                    integ->fsal_valid = 1;
                }
                for (int i = 0; i < dim; i++) {
                    k[i] = integ->f_now[i] * integ->h;
                }
            } else if (j == 6) {
                /* Stage 7 sits at (t+h, y_new) because a[6] == b: keep
                 * the raw derivative, it becomes f_now if the step is
                 * accepted. */
                returnNumber = integ->rhs(temp_clock, temp, integ->f_new, integ->user);
                integ->n_rhs++;
                if (returnNumber != 0) return SPODY_INTEG_ERR_RHS;
                for (int i = 0; i < dim; i++) {
                    k[6*dim + i] = integ->f_new[i] * integ->h;
                }
            } else {
                returnNumber = integ->rhs(temp_clock, temp, k + j*dim, integ->user);
                integ->n_rhs++;
                if (returnNumber != 0) return SPODY_INTEG_ERR_RHS;
                for (int i = 0; i < dim; i++) {
                    k[j*dim + i] *= integ->h;   // save the K factor already multiplied by h
                }
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
            /* FSAL hand-over: the derivative at the new state is the
             * one stage 7 just computed. Swap the buffers, no copy. */
            {
                double *swap  = integ->f_now;
                integ->f_now  = integ->f_new;
                integ->f_new  = swap;
            }
            integ->fsal_valid = 1;
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

    spody_log_eprintf("[RKDP45_05][MAX_ITER] !Iteration limit reached! ---> Bad constraints\n");
    return SPODY_INTEG_ERR_MAX_STEPS;
}

//----- DOP853 (Dormand-Prince 8(5,3), adaptive) --------------------------

/*
 * One adaptive DOP853 step; same contract as step_rkdp45.
 *
 * Error: per 3-component block, with the RSS-step scaling of
 * step_rkdp45 (err5, err3 divided by the RSS of the block's change over
 * the step when its square exceeds 0.1),
 *     err = err5^2 / sqrt(err5^2 + 0.01 err3^2)
 * (Hairer, Norsett, Wanner 1993, Sect. II.10). The step is accepted
 * when the worst block is below rel_tol.
 *
 * Cost: eleven RHS calls per attempt plus one per accepted step,
 * f(t + h, y_new).
 */
static int step_dop853(IntegratorAllData *integ) {

    const int dim   = integ->dim;
    double *state   = integ->y;
    double *temp    = integ->y_tmp;     // stage state, then y_new
    double *k       = integ->k;         // flat buffer: k[stage*dim + i] = h * f

    const double clock = integ->t;
    int steps    = 0;
    int rejected = 0;                   // a rejection happened in this call

    memcpy(integ->y_old, state, sizeof(double) * (size_t)dim);

    do {

        const double h = integ->h;

        /* Stage 1 is f(t, y): FSAL from the previous accepted step, or
         * from an earlier attempt of this one (the state has not moved). */
        if (!integ->fsal_valid) {
            int rc = integ->rhs(clock, state, integ->f_now, integ->user);
            integ->n_rhs++;
            if (rc != 0) return SPODY_INTEG_ERR_RHS;
            integ->fsal_valid = 1;
        }
        for (int i = 0; i < dim; i++) k[i] = integ->f_now[i] * h;

        for (int j = 1; j < DOP853_STAGES; j++) {
            for (int i = 0; i < dim; i++) {
                double acc = state[i];
                for (int x = 0; x < j; x++) {
                    if (dop853_a[j][x] != 0.0) acc += dop853_a[j][x] * k[x*dim + i];
                }
                temp[i] = acc;
            }
            double *kj = k + j*dim;
            int rc = integ->rhs(clock + h * dop853_c[j], temp, kj, integ->user);
            integ->n_rhs++;
            if (rc != 0) return SPODY_INTEG_ERR_RHS;
            for (int i = 0; i < dim; i++) kj[i] *= h;
        }

        /* y_new and the worst block of the stretched error estimate. */
        double err = 0.0;
        for (int i = 0; i + 2 < dim; i += 3) {
            double d2 = 0.0, e5_2 = 0.0, e3_2 = 0.0;
            for (int m = i; m < i + 3; m++) {
                double incr = 0.0, e5 = 0.0, bhh = 0.0;
                for (int s = 0; s < DOP853_STAGES; s++) {
                    const double ks = k[s*dim + m];
                    incr += dop853_b[s]   * ks;
                    e5   += dop853_e5[s]  * ks;
                    bhh  += dop853_bhh[s] * ks;
                }
                temp[m] = state[m] + incr;
                const double e3 = incr - bhh;
                d2   += incr * incr;
                e5_2 += e5 * e5;
                e3_2 += e3 * e3;
            }
            if (d2 > 0.1) { e5_2 /= d2; e3_2 /= d2; }
            const double deno = e5_2 + 0.01 * e3_2;
            const double err_i = (deno > 0.0) ? e5_2 / sqrt(deno) : 0.0;
            if (err_i > err) err = err_i;
        }
        /* a trailing component outside any 3-block (dim % 3 != 0) is
         * integrated but, as in step_rkdp45, not error-controlled */
        for (int m = dim - dim % 3; m < dim; m++) {
            double incr = 0.0;
            for (int s = 0; s < DOP853_STAGES; s++) incr += dop853_b[s] * k[s*dim + m];
            temp[m] = state[m] + incr;
        }

        steps++;

        #if DEBUG_INTEGRATORS == 1
        printf("[dop853_01] iter : %d | t : %.6f | h : %.6e | err : %.2e | tol : %.2e\n",
               steps, clock, h, err, integ->opt.rel_tol);
        #endif

        const double fac = integ->opt.safety * pow(integ->opt.rel_tol / err, 1.0 / 8.0);

        if (err < integ->opt.rel_tol) {

            /* FSAL: f(t + h, y_new) is the next step's stage 1. Evaluated
             * before committing so an RHS failure leaves the state as it
             * was. */
            int rc = integ->rhs(clock + h, temp, integ->f_new, integ->user);
            integ->n_rhs++;
            if (rc != 0) return SPODY_INTEG_ERR_RHS;

            memcpy(state, temp, sizeof(double) * (size_t)dim);
            {
                double *swap  = integ->f_now;
                integ->f_now  = integ->f_new;
                integ->f_new  = swap;
            }
            integ->fsal_valid = 1;
            integ->h_old = h;
            integ->t_old = clock;
            integ->t     = clock + h;

            /* no growth right after a rejection */
            double scale = fmin(DOP853_FACTOR_UPSCALE, fac);
            if (scale < DOP853_FACTOR_DOWNSCALE) scale = DOP853_FACTOR_DOWNSCALE;
            if (rejected && scale > 1.0) scale = 1.0;
            integ->h = h * scale;
            if (integ->h > integ->opt.h_max) integ->h = integ->opt.h_max;

            integ->n_accepted++;
            return SPODY_INTEG_OK;

        } else {

            integ->n_rejected++;
            rejected = 1;

            integ->h = h * fmax(DOP853_FACTOR_DOWNSCALE, fac);
            if (fabs(integ->h) < integ->opt.h_min) {
                return SPODY_INTEG_ERR_STEP_TOO_SMALL;
            }
        }

    } while (steps < DOP853_RETRIES_PER_STEP);

    spody_log_eprintf("[DOP853][MAX_ITER] !Iteration limit reached! ---> Bad constraints\n");
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

//----- public step / drive ----------------------------------------------

int spody_propagate_onestep(IntegratorAllData *integ) {
    if (!integ || !integ->rhs || !integ->y) return SPODY_INTEG_ERR_NULL;
    switch (integ->method) {
        case SPODY_INTEG_RK4:    return step_rk4(integ);
        case SPODY_INTEG_RK45:   return step_rkdp45(integ);
        case SPODY_INTEG_DOP853: return step_dop853(integ);
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
 * downstream localisation -- event checks run Brent on the dense
 * curve, so a non-integrator-consistent interpolant gives a wrong
 * t_trigger. (They now use spody_dense_state_rv6, the quintic on the
 * same endpoint data plus the accelerations; the reasoning holds for
 * both.)
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

int spody_dense_state_rv6(const IntegratorAllData *integ, double t,
                          double y_out[6]) {
    if (!integ || !y_out) return SPODY_INTEG_ERR_NULL;
    /* RK45 and DOP853: the quintic needs only the endpoint states and
     * accelerations. DOP853's own dense output (dop853.f, CONTD8) is
     * not implemented yet. */
    if ((integ->method != SPODY_INTEG_RK45 && integ->method != SPODY_INTEG_DOP853)
        || integ->dim != 6)
        return SPODY_INTEG_ERR_NULL;
    if (!integ->y_old || !integ->y || !integ->f_now || !integ->f_new)
        return SPODY_INTEG_ERR_NULL;

    /* The accepted step swapped the FSAL buffers: f_now now holds
     * f(t, y) at the END of the step and f_new f(t_old, y_old) at its
     * START -- the names describe the next step, not this one.
     * Components 3..5 of f are the accelerations. */
    const double *f_start = integ->f_new;
    const double *f_end   = integ->f_now;

    spody_hermite_quintic_rv6(t,
                              integ->t_old, integ->y_old, f_start + 3,
                              integ->t,     integ->y,     f_end + 3,
                              y_out);
    return SPODY_INTEG_OK;
}
