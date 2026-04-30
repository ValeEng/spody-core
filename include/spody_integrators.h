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
#ifndef SPODY_INTEGRATORS_H
#define SPODY_INTEGRATORS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

    //debug
#define DEBUG_INTEGRATORS 0 // 0 = no debug | 1 = debug |---> CODE TESTING

#define SPODY_INTEG_DEFAULT_DIM 6 // default state size: position(3) + velocity(3)

/* Return codes for the integrator step / drive functions. */
#define SPODY_INTEG_OK              0
#define SPODY_INTEG_ERR_NULL       -1
#define SPODY_INTEG_ERR_DIM        -2
#define SPODY_INTEG_ERR_RHS        -3
#define SPODY_INTEG_ERR_STEP_TOO_SMALL -4
#define SPODY_INTEG_ERR_MAX_STEPS  -5

/*
 * Right-hand side (RHS) callback.
 *
 *   t       : current independent variable (typically seconds from epoch).
 *   y       : current state vector. Size matches the `dim` passed to
 *             spody_setup_integrator and is implicit in the contract.
 *   dy      : output, dy/dt evaluated at (t, y). Same size as `y`.
 *   user    : opaque user pointer (forces context: gravity model, ephemeris
 *             handle, satellite parameters, and -- if a generic RHS needs it --
 *             the state size itself).
 *
 * Return 0 on success, non-zero to abort the integration step.
 *
 * The RHS must be reentrant: the integrator may call it multiple times per
 * step (Runge-Kutta stages, error estimate). It must not modify `y`.
 */
typedef int (*spody_rhs_fn)(double t, const double *y, double *dy, void *user);

/*
 * Available integrator methods. Fixed-step methods ignore tolerance fields,
 * adaptive methods ignore the requested step size after the first stage.
 */
typedef enum {
    SPODY_INTEG_RK4    = 0,  // classical 4th-order Runge-Kutta, fixed step
    SPODY_INTEG_RK45   = 1,  // Dormand-Prince 5(4), adaptive step
    SPODY_INTEG_RK78   = 2,  // Fehlberg 7(8), adaptive step (high order)
    SPODY_INTEG_VERLET = 3   // velocity-Verlet, fixed step (symplectic, 2nd order)
} spody_integrator_method;

/*
 * Per-method tuning. Not every field is read by every method.
 *
 *   h_init       : initial step size (s). For fixed-step methods, the step.
 *   h_min, h_max : bounds on the adaptive step size (s). Ignored if fixed.
 *   abs_tol      : absolute tolerance per state component (adaptive only).
 *   rel_tol      : relative tolerance per state component (adaptive only).
 *   max_steps    : safety cap on the number of steps per drive call. 0 = unlimited.
 *   safety       : safety factor on adaptive step update (typical 0.8 - 0.9).
 */
typedef struct {
    double h_init;
    double h_min;
    double h_max;
    double abs_tol;
    double rel_tol;
    size_t max_steps;
    double safety;
} IntegratorOptions;

/*
 * Per-thread integrator workspace.
 *
 * Holds the state vector, the RHS callback + user data, the chosen method,
 * its options, and scratch buffers for the Runge-Kutta stages. Sized for a
 * state of `dim` doubles.
 *
 * Threading model (mirrors MappedEphemeris / HarmonicGravity):
 *   - One spody_integrator per thread.
 *   - Multiple integrators can share the same RHS user-pointer payload as
 *     long as that payload is itself thread-safe (e.g. shared
 *     MappedEphemerisData + per-thread MappedEphemeris).
 */
typedef struct {
    spody_integrator_method method;
    IntegratorOptions opt;

    int dim;                      // state vector size
    double t;                     // current independent variable
    double *y;                    // current state, size `dim`
    double h;                     // next step size proposed by the controller
    double h_old;                 // step size actually used in the last accepted step

    spody_rhs_fn rhs;             // user-provided dynamics
    void *user;                   // opaque payload passed to rhs

    // scratch buffers, sized at setup time, reused across steps
    double *k;
    double *y_tmp;                // intermediate state for stage evaluation
    double *y_err; 
    /*
    double *k1, *k2, *k3, *k4;    // RK4 / RK45 base stages
    double *k5, *k6, *k7;         // additional stages (RK45 / RK78)
    double *k8, *k9, *k10, *k11, *k12, *k13; // RK78 high-order stages
               // error estimate (adaptive methods)
    */
} IntegratorAllData;

/*
 * Build a default options block for a method. Caller may then tweak fields
 * before passing it to spody_setup_integrator.
 */
void spody_default_integrator_options(spody_integrator_method method, IntegratorOptions *opt);

/*
 * Setup a per-thread integrator workspace.
 *
 *   integ   : workspace to initialize (caller-allocated).
 *   method  : integrator method.
 *   opt     : options (copied by value).
 *   dim     : state vector size (e.g. 6 for position+velocity).
 *   rhs     : dynamics callback.
 *   user    : opaque payload forwarded to rhs.
 *
 * Allocates internal scratch buffers sized to `dim`. Pair with
 * spody_free_integrator. Returns SPODY_INTEG_OK on success.
 */
int spody_setup_integrator(IntegratorAllData *integ,
                           spody_integrator_method method,
                           const IntegratorOptions *opt,
                           int dim,
                           spody_rhs_fn rhs,
                           void *user);

int spody_free_integrator(IntegratorAllData *integ);

/*
 * Set / reset the integration state without reallocating buffers.
 *
 *   t0  : initial independent variable.
 *   y0  : initial state vector, size matching the workspace `dim`. Copied in.
 */
int spody_set_integrator_state(IntegratorAllData *integ, double t0, const double *y0);

/*
 * Advance the state by exactly one step.
 *
 * For fixed-step methods the step size is integ->h. For adaptive methods the
 * step size may shrink (rejection) or grow (success); the final accepted h
 * is left in integ->h and integ->t/integ->y are updated accordingly.
 *
 * Returns SPODY_INTEG_OK on success.
 */
int spody_propagate_onestep(IntegratorAllData *integ);

/*
 * Drive the integration up to t_end (t_end may be in the past: negative h
 * is supported). The last step is clipped to land exactly on t_end.
 *
 * Returns SPODY_INTEG_OK if t_end is reached within opt.max_steps,
 * SPODY_INTEG_ERR_MAX_STEPS otherwise.
 */
int spody_propagate_untilend(IntegratorAllData *integ, double t_end);

#ifdef __cplusplus
}
#endif

#endif // SPODY_INTEGRATORS_H
