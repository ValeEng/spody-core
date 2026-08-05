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
 * Optional per-step callback, invoked ONCE before the stages of each
 * step -- and NOT again on the internal retries that a rejected trial
 * step triggers.
 *
 *   t    : independent variable at the START of the step.
 *   y    : state at the start of the step. Must not be modified.
 *   h    : step size about to be attempted. Retries only ever shrink
 *          it, so a decision taken for this h stays valid for them.
 *   user : the same opaque payload the RHS receives.
 *
 * Return 0 to proceed, non-zero to abort the step.
 *
 * It exists so a force model can retune ITSELF for the coming step --
 * typically to trade evaluation cost against an accuracy it can only
 * judge from the current state. Retuning inside the RHS instead would
 * be wrong: the stages would then sample different vector fields, the
 * embedded error estimate would read the model jump as truncation
 * error, and the controller would shrink h fighting a discontinuity
 * that is not in the dynamics. Deciding here, once, keeps every stage
 * of a step on ONE smooth field -- the property the Runge-Kutta order
 * derivation assumes.
 *
 * NULL (the value left by spody_setup_integrator) disables it, so
 * existing callers are unaffected bit for bit.
 */
typedef int (*spody_prestep_fn)(double t, const double *y, double h,
                                void *user);

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
    double t;                     // current independent variable (after the last step)
    double *y;                    // current state, size `dim`
    double h;                     // next step size proposed by the controller
    double h_old;                 // step size actually used in the last accepted step
    double t_old;                 // t at the START of the last accepted step (= t - h_old)
    double *y_old;                // state at the START of the last accepted step

    spody_rhs_fn rhs;             // user-provided dynamics
    void *user;                   // opaque payload passed to rhs
    spody_prestep_fn pre_step;    // optional per-step hook, NULL = none

    // scratch buffers, sized at setup time, reused across steps
    double *k;                    // RK stages, kept alive after step accept for dense output
    double *y_tmp;                // intermediate state for stage evaluation
    double *y_err;

    /* Cost counters. Zeroed by spody_setup_integrator, monotonically
     * increasing afterwards; read them whenever, typically once the run
     * is over. They measure the integrator's work in a way that does
     * not depend on machine, compiler or system load, which is what
     * makes a cost comparable across runs and across tools -- wall
     * clock alone is not.
     *
     * n_rhs counts every call to the RHS callback, including the ones
     * spent on trial steps that were later rejected: that work really
     * was done, and leaving it out would understate the cost. */
    size_t n_accepted;            // accepted steps
    size_t n_rejected;            // rejected trial steps
    size_t n_rhs;                 // RHS evaluations

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

/* ============================================================
 * Dense output (continuous interpolation inside the last accepted step)
 *
 * After spody_propagate_onestep returns SPODY_INTEG_OK, this function
 * can be called any number of times to evaluate the state at any point
 * within the just-completed interval [t_old, t_old + h_old]:
 *
 *     y(theta) = y(t_old + theta * h_old),    theta in [0, 1]
 *
 * Implementation (SPODY_INTEG_RK45):
 *   Cubic Hermite C^1 using the FSAL endpoint derivatives k_1 and k_7
 *   (which already sit in integ->k pre-multiplied by h_old). This is
 *   integrator-consistent for the 7S Butcher tableau we run, so a
 *   downstream root-finder (e.g. spody_event_check_refined) sees a
 *   curve that actually corresponds to the integrated trajectory.
 *   See the file-level comment in spody_integrators.c for why we use
 *   Hermite here instead of the classical DOPRI5 P-matrix.
 *
 *   Accuracy is 4th-order on the state (one order below the integrator
 *   itself); over a 30 s step at LRO this localises a surface crossing
 *   to well under one microsecond -- below any physically meaningful
 *   threshold for impact / altitude / eclipse events.
 *
 * Other methods (RK4, Verlet, RK78) currently return
 * SPODY_INTEG_ERR_NULL; their dense-output formulas can be added later.
 *
 * `theta` is clamped to [0, 1] internally. */
int spody_dense_eval(const IntegratorAllData *integ, double theta, double *y_out);

#ifdef __cplusplus
}
#endif

#endif // SPODY_INTEGRATORS_H
