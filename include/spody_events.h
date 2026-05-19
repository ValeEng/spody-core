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
#ifndef SPODY_EVENTS_H
#define SPODY_EVENTS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spody_forcemodels.h"   /* ForceModelContext */
#include "spody_integrators.h"   /* IntegratorAllData (for dense output)   */

/* ============================================================
 * Orbital events
 * ============================================================
 * An "event" is a predicate evaluated on the integrator state after
 * each accepted step. When the predicate fires, the runtime can:
 *   - log the trigger (action = LOG / LOG_AND_STOP);
 *   - stop the propagation (action = STOP / LOG_AND_STOP).
 *
 * Events are designed to be cheap (a handful of FLOPs + at most one
 * ephemeris query per check) and live OUTSIDE the RHS hot path: they
 * are evaluated by the drive loop (e.g. spody_run_propagate) on the
 * accepted state, never by the integrator's per-stage callback.
 *
 * Localization precision:
 *   The event fires at the first accepted step where the predicate
 *   becomes true; the timestamp of the trigger is therefore the
 *   integrator timestamp at that step (precision ~ accepted h, of
 *   order 30 s on LLO with rel_tol = 1e-9). Sub-step root-finding via
 *   dense output is a planned extension but not implemented today.
 *
 * Threading:
 *   SpodyEvent is small and copyable. Mission keeps a caller-owned
 *   array of events; each Mission (per-thread) has its own array.
 * ============================================================ */

/* Kinds of supported events. The enum is open: adding new kinds (e.g.
 * altitude, eclipse entry/exit, apsis passage) is non-breaking. */
typedef enum {
    SPODY_EVENT_KIND_IMPACT = 0       /* |r_sat - r_body| < radius_km */
    /* future: SPODY_EVENT_KIND_ALT_BELOW, _APSIS, _ECLIPSE_ENTRY, ... */
} spody_event_kind;

/* What the runtime should do when an event triggers. */
typedef enum {
    SPODY_EVENT_ACTION_LOG          = 0,  /* record the trigger, keep going */
    SPODY_EVENT_ACTION_STOP         = 1,  /* stop the propagation, no log   */
    SPODY_EVENT_ACTION_LOG_AND_STOP = 2   /* both                            */
} spody_event_action;

/* Generic event descriptor.
 *
 * IMPACT semantics:
 *   The check is geometric: at each step, compute the satellite-to-body
 *   distance and compare with `radius_km`. The body is identified by
 *   `naif_id`. If `naif_id == ctx->naif_central` the satellite position
 *   y[0..2] is already in the central frame; otherwise the body
 *   position is queried from the ephemeris in the central frame and
 *   subtracted from y[0..2]. Trigger condition:
 *
 *      |r_sat - r_body| < radius_km
 *
 *   At trigger the runtime fills the bookkeeping fields below.
 */
typedef struct {
    /* ---- caller-set ---- */
    spody_event_kind   kind;
    spody_event_action action;
    int    naif_id;                 /* IMPACT: body to check against */
    double radius_km;               /* IMPACT: trigger threshold     */

    /* ---- runtime-set (output) ---- */
    int    triggered;               /* 1 if the predicate has fired              */
    double t_trigger;               /* sim time of the trigger (seconds from t=0) */
    double y_trigger[6];            /* state at trigger: r(3), v(3)              */
    double distance_at_trigger;     /* IMPACT: |r_sat - r_body| at trigger        */

    /* ---- runtime-set (internal book-keeping for sign tracking) ----
     * Distance value at the previous step. Used to detect a sign change
     * of (distance - radius) between two accepted steps; root-finding is
     * launched only when the predicate flips. Sentinel value -1.0 means
     * "no previous value yet" (first call). */
    double prev_distance_signed;    /* (|r_sat - r_body| - radius) at last step  */
    int    prev_valid;              /* 0 = first call, 1 = prev_distance valid   */
} SpodyEvent;

/* Convenience constructor for an impact event. */
SpodyEvent spody_event_impact(int naif_id, double radius_km, spody_event_action action);

/* Evaluate the event predicate on the given state.
 *
 * Return value:
 *   1 -> NEW fire in this call. The trigger output fields
 *        (t_trigger, y_trigger, distance_at_trigger) have just been
 *        written; ev->triggered is now 1. The caller should react
 *        (log to file, request stop, etc.).
 *   0 -> no new fire. Either the predicate did not hold, or the event
 *        is a one-shot kind that already latched on a previous call
 *        (per-kind decision inside the dispatch -- IMPACT is one-shot,
 *        future recurring kinds like altitude / apsis / eclipse won't
 *        latch). The output fields are NOT updated in this case; if
 *        you need "has this event ever fired?", read ev->triggered.
 *
 * "Coarse" version: fires as soon as the predicate is true at the
 * accepted step (precision = step size). */
int spody_event_check(SpodyEvent *ev,
                      const ForceModelContext *ctx,
                      double t, const double *y);

/* "Refined" version: detects a sign change of the predicate between two
 * accepted steps and uses dense output + Brent root-finding to localise
 * the trigger inside the last step. Precision is at the level of the
 * dense interpolant (cubic Hermite C^1; sub-microsecond on a 30 s LRO
 * step). Same return semantics as the coarse version above.
 *
 * The integrator must be SPODY_INTEG_RK45 (other methods don't yet
 * provide dense_eval and silently fall back to the coarse path).
 * Pass `integ` as the integrator that just produced (t, y) via
 * spody_propagate_onestep -- the function reads integ->t_old, integ->y_old
 * and the stored RK stages to evaluate the interpolant. */
int spody_event_check_refined(SpodyEvent *ev,
                              const ForceModelContext *ctx,
                              const IntegratorAllData *integ);

/* On-disk record for the events log (magic SPDYEVNT). One record is
 * appended every time an event triggers; the runtime writes this struct
 * verbatim. 80 bytes total, no internal padding on standard ABIs. */
typedef struct {
    double t;                       /* sim time (s)                       */
    int    kind;                    /* spody_event_kind                   */
    int    naif_id;                 /* body involved                      */
    double radius_km;               /* threshold used                     */
    double distance_km;             /* distance at trigger                */
    double y[6];                    /* state at trigger                   */
} EventRecord;

#ifdef __cplusplus
}
#endif

#endif /* SPODY_EVENTS_H */
