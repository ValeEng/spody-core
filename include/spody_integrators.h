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
    SPODY_INTEG_DOP853 = 2   // Dormand-Prince 8(5,3), adaptive step
} spody_integrator_method;

/*
 * Per-method tuning. Not every field is read by every method.
 *
 *   h_init       : initial step size (s). For fixed-step methods, the step.
 *   h_min, h_max : bounds on the adaptive step size (s). Ignored if fixed.
 *   rel_tol      : step tolerance (adaptive only). RKDP45 splits the state
 *                  into 3-component blocks (position, velocity) and, per
 *                  block, takes the RSS of the embedded error estimate,
 *                  divided by the RSS of the step's own change of that
 *                  block when its square exceeds 0.1 (km^2 or km^2/s^2),
 *                  absolute otherwise; the step is accepted when the worst
 *                  block is below rel_tol (GMAT's RSS-step control). There
 *                  is no separate absolute tolerance. DOP853 applies the
 *                  same block norm to err5^2 / sqrt(err5^2 + 0.01 err3^2)
 *                  of its order-5 and order-3 estimates.
 *   max_steps    : safety cap on the number of steps per drive call. 0 = unlimited.
 *   safety       : safety factor on adaptive step update (typical 0.8 - 0.9).
 */
typedef struct {
    double h_init;
    double h_min;
    double h_max;
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

    // scratch buffers, sized at setup time, reused across steps
    double *k;                    // RK stages, kept alive after step accept for dense output
    double *y_tmp;                // intermediate state for stage evaluation
    double *y_err;

    /* First-same-as-last reuse, SPODY_INTEG_RK45 and SPODY_INTEG_DOP853:
     * the 7S tableau evaluates its last stage at (t+h, y_new), which is
     * exactly the first stage of the step after, so one RHS call per
     * accepted step is redundant; DOP853 evaluates f(t+h, y_new) once
     * the step is accepted and reuses it the same way. f_now holds
     * f(t, y) at the current state, unscaled by h; f_new is the scratch the
     * last stage writes into; the two swap on acceptance. fsal_valid
     * says whether f_now may seed the next step instead of a fresh RHS
     * call: set once the derivative at the current state has been
     * evaluated, cleared by spody_set_integrator_state and
     * spody_integrator_invalidate_fsal. NULL / 0 for every other method,
     * none of which reads them. */
    double *f_now;
    double *f_new;
    int     fsal_valid;

    /* Cost counters. Zeroed by spody_setup_integrator, monotonically
     * increasing afterwards; read them whenever, typically once the run
     * is over. They measure the integrator's work in a way that does
     * not depend on machine, compiler or system load, which is what
     * makes a cost comparable across runs and across tools -- wall
     * clock alone is not.
     *
     * n_rhs counts every call to the RHS callback, including the ones
     * spent on trial steps that were later rejected: that work really
     * was done, and leaving it out would understate the cost. With
     * RK45 that is six calls per attempted step plus one for the very
     * first stage of a run (or after an invalidation), not seven: the
     * seventh is the FSAL derivative carried over from the step
     * before. With DOP853 it is eleven calls per attempted step plus
     * one per accepted step (the FSAL derivative at the new state),
     * plus one for the very first stage. */
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
 *
 * Also drops the FSAL derivative (see spody_integrator_invalidate_fsal):
 * the next step re-evaluates the RHS at (t0, y0).
 */
int spody_set_integrator_state(IntegratorAllData *integ, double t0, const double *y0);

/*
 * Forget the derivative the RK45 step keeps from the previous step, so
 * the next step evaluates its first stage afresh.
 *
 * The FSAL reuse assumes the RHS is the same function of (t, y) at the
 * start of a step as it was at the end of the previous one. Call this
 * whenever that stops being true between two steps: the vector field
 * was retuned (an adaptive harmonics degree that moved, a spacecraft
 * parameter that changed, a force switched on or off) or the state was
 * edited in place rather than through spody_set_integrator_state.
 *
 * Costs one extra RHS evaluation on the next step and nothing else.
 * A no-op for methods without the FSAL property.
 */
void spody_integrator_invalidate_fsal(IntegratorAllData *integ);

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
 *   downstream root-finder sees a curve that actually corresponds to
 *   the integrated trajectory. For (r, v) states -- the event
 *   localisation and the fixed output grid -- spody_dense_state_rv6
 *   below is the one to use: quintic, it keeps the velocity at the
 *   integrator's accuracy.
 *   See the file-level comment in spody_integrators.c for why we use
 *   Hermite here instead of the classical DOPRI5 P-matrix.
 *
 *   The interpolation error is O(h^4), one order below the integrator
 *   itself, and on the velocity one order lower still: that is why the
 *   engine's own grid and events use spody_dense_state_rv6.
 *
 * RK4 and DOP853 return SPODY_INTEG_ERR_NULL (k_7 is RK45's FSAL
 * stage; DOP853 goes through spody_dense_state_rv6).
 *
 * `theta` is clamped to [0, 1] internally. */
int spody_dense_eval(const IntegratorAllData *integ, double theta, double *y_out);

/* Dense output for a second-order system in (r, v) layout -- dim 6,
 * y = [r(3), v(3)], f = (v, a) -- at time t inside the last accepted
 * step (clamped to [t_old, t]).
 *
 * Quintic Hermite (spody_hermite_quintic_rv6) on r, v and the
 * accelerations at both ends of the step. Those accelerations cost
 * nothing: they are the FSAL derivatives the step already holds,
 * unscaled by h -- after the step f(t_old, y_old) and f(t, y) --, so
 * no RHS evaluation is added and no k/h division rounds them. Position
 * and velocity keep the integrator's own accuracy, where the cubic of
 * spody_dense_eval loses an order on the velocity.
 *
 * Valid right after spody_propagate_onestep returns SPODY_INTEG_OK and
 * until the next step or state reset. RK45 and DOP853; returns
 * SPODY_INTEG_ERR_NULL otherwise. For DOP853 the quintic is of lower
 * order than the method; its own dense output is not implemented yet. */
int spody_dense_state_rv6(const IntegratorAllData *integ, double t,
                          double y_out[6]);

#ifdef __cplusplus
}
#endif

#endif // SPODY_INTEGRATORS_H
