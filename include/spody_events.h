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
 *   One-shot kinds (impact) fire at the first accepted step where
 *   the predicate becomes true; the timestamp is the integrator
 *   timestamp at that step (precision ~ accepted h, of order 30 s on
 *   LLO with rel_tol = 1e-9). Recurring kinds (eclipse, altitude
 *   crossing) refine the crossing inside the accepted step via Brent
 *   root-finding on the RK45 dense output (spody_event_check_refined,
 *   sub-microsecond localization); per-event opt-out via the
 *   `refined` flag falls back to step-boundary precision.
 *
 * Threading:
 *   SpodyEvent is small and copyable. Mission keeps a caller-owned
 *   array of events; each Mission (per-thread) has its own array.
 * ============================================================ */

/* Kinds of supported events. The enum is open: adding new kinds (e.g.
 * apsis passage, geodetic-region entry) is non-breaking -- the
 * dispatch switches in spody_event_check[_refined] add one case
 * each, the wire format stays as-is. */
typedef enum {
    SPODY_EVENT_KIND_IMPACT       = 0, /* one-shot: |r_sat - r_body| < radius_km   */
    SPODY_EVENT_KIND_ECLIPSE      = 1, /* recurring: eclipse fraction crosses threshold */
    SPODY_EVENT_KIND_ALT_CROSSING = 2  /* recurring: altitude wrt a body crosses a target value (asc + desc) */
    /* future: SPODY_EVENT_KIND_APSIS, _GEODETIC_REGION, ... */
} spody_event_kind;

/* What the runtime should do when an event triggers. */
typedef enum {
    SPODY_EVENT_ACTION_LOG          = 0,  /* record the trigger, keep going */
    SPODY_EVENT_ACTION_STOP         = 1,  /* stop the propagation, no log   */
    SPODY_EVENT_ACTION_LOG_AND_STOP = 2   /* both                            */
} spody_event_action;

/* Generic event descriptor.
 *
 * Field meanings depend on `kind` -- the same slots are reused across
 * predicates to keep the struct small and cheap to copy.
 *
 * IMPACT (one-shot):
 *   The check is geometric: compute the satellite-to-body distance and
 *   compare with `radius_km`. The body is identified by `naif_id`. If
 *   `naif_id == ctx->naif_central` the satellite position y[0..2] is
 *   already in the central frame; otherwise the body position is
 *   queried from the ephemeris in the central frame and subtracted
 *   from y[0..2]. Trigger condition:
 *
 *      |r_sat - r_body| < radius_km
 *
 *   At trigger the runtime writes |r_sat - r_body| into
 *   distance_at_trigger.
 *
 * ECLIPSE (recurring):
 *   `naif_id` is the occulting body, `radius_km` its physical radius
 *   (both must be set by the caller); `threshold_fraction` is the
 *   eclipse-fraction level whose crossing fires the event:
 *     1.0 -> any loss of sunlight (penumbra entry)
 *     0.5 -> middle of penumbra (default in spody CLI)
 *     0.0 -> full umbra entry
 *   The fraction is computed by spody_get_satlitfraction (Montenbruck
 *   & Gill) against THIS event's occulter alone -- one event per
 *   occulting body, unlike the SRP force, which combines every
 *   occulter into a single lit fraction. During a double eclipse the
 *   two therefore disagree on purpose. The Sun position is always
 *   queried from the ephemeris.
 *   At trigger the runtime writes the fraction at trigger into
 *   distance_at_trigger (semantic slot reuse).
 *
 *   ECLIPSE is recurring: every threshold crossing (entry AND exit)
 *   fires a fresh trigger. The caller distinguishes the direction by
 *   the value of prev_distance_signed at the previous fire, or by
 *   post-processing the events log.
 *
 * ALT_CROSSING (recurring):
 *   Trigger condition: |r_sat - r_body| - radius_km - altitude_km = 0
 *   where `radius_km` is the body's physical radius and `altitude_km`
 *   the user's target altitude above its surface. Fires on every
 *   sign change of the residual -> ascending AND descending crossings
 *   are both logged. `naif_id` selects the body (central or any third
 *   body whose position the ephemeris can resolve; the impact-style
 *   `has_ref_point` shortcut also works for CR3BP primaries which sit
 *   at fixed synodic positions).
 *   At trigger the runtime writes the satellite-to-body distance at
 *   trigger into distance_at_trigger; subtract `radius_km` to get the
 *   altitude actually attained (= altitude_km up to solver tolerance).
 *   Direction (ascending / descending) is recoverable from the radial
 *   component of y_trigger[3..5].
 */
typedef struct {
    /* ---- caller-set ---- */
    spody_event_kind   kind;
    spody_event_action action;
    int    naif_id;                 /* IMPACT:       body to check against
                                       ECLIPSE:      occulting body
                                       ALT_CROSSING: body to measure altitude from */
    double radius_km;               /* IMPACT:       trigger threshold
                                       ECLIPSE:      occulter physical radius
                                       ALT_CROSSING: body's physical radius        */
    double threshold_fraction;      /* ECLIPSE: fraction threshold in [0, 1] */
    double altitude_km;             /* ALT_CROSSING: target altitude above body
                                       surface (km); irrelevant to other kinds    */

    /* Per-event refinement opt-out. 1 (default) = Brent + dense output
     * to localise the trigger inside the accepted step (precision
     * sub-microsecond on a 30 s LRO step). 0 = trigger at the end of
     * the step where the predicate first holds (precision = step
     * size). The Brent path is essentially free except at the actual
     * crossing step, so the toggle exists for users who have hundreds
     * of altitude bands and want to skip the sign-change overhead. */
    int    refined;

    /* Optional explicit reference point for IMPACT and ALT_CROSSING,
     * in the integrator's working frame. When `has_ref_point` is
     * non-zero the distance check tests |y_xyz - ref_point| instead
     * of falling back to the central body (origin) or an ephemeris
     * query. Used by CR3BP where the two primaries sit at fixed
     * synodic-frame positions and there is no ephemeris. Default
     * 0 / {0,0,0} preserves the legacy HF behaviour. */
    int    has_ref_point;
    double ref_point[3];

    /* ---- runtime-set (output) ---- */
    int    triggered;               /* 1 if the predicate has ever fired        */
    double t_trigger;               /* sim time of the (last) trigger           */
    double y_trigger[6];            /* state at the (last) trigger              */
    double distance_at_trigger;     /* IMPACT: |r_sat - r_body| at trigger
                                       ECLIPSE: eclipse fraction at trigger      */

    /* ---- runtime-set (internal book-keeping for sign tracking) ----
     * Predicate scalar at the previous step. Used to detect a sign
     * change of the predicate between two accepted steps; the refined
     * root-finding is launched only when it flips. */
    double prev_distance_signed;    /* IMPACT: (|r_sat - r_body| - radius)
                                       ECLIPSE: (eclipse_fraction - threshold)   */
    int    prev_valid;              /* 0 = first call, 1 = prev_distance valid   */
} SpodyEvent;

/* Convenience constructor for an impact event. */
SpodyEvent spody_event_impact(int naif_id, double radius_km, spody_event_action action);

/* Constructor for an impact event with a fixed reference point in the
 * integrator's working frame. The distance check is
 *   |y_xyz - ref_point| < radius_km
 * regardless of `ctx->naif_central` or the ephemeris (both of which may
 * be absent, as in CR3BP). `naif_id` is carried purely for logging /
 * downstream identification. */
SpodyEvent spody_event_impact_at_point(int naif_id, const double ref_point[3],
                                        double radius_km, spody_event_action action);

/* Convenience constructor for an eclipse event.
 *
 *   occulter_naif_id    : NAIF id of the occulting body
 *   occulter_radius_km  : its physical radius
 *   threshold_fraction  : eclipse-fraction level whose crossing fires
 *                         the event (see ECLIPSE section above)
 *   action              : typically SPODY_EVENT_ACTION_LOG since
 *                         eclipse events are usually informational
 */
SpodyEvent spody_event_eclipse(int occulter_naif_id, double occulter_radius_km,
                               double threshold_fraction, spody_event_action action);

/* Convenience constructor for an altitude-crossing event (HF path).
 *
 *   naif_id         : body to measure altitude from. May be
 *                     ctx->naif_central (origin in HF) or any third
 *                     body whose position the ephemeris can resolve.
 *   body_radius_km  : physical radius of `naif_id`; the predicate is
 *                     `|r_sat - r_body| - body_radius_km - altitude_km`.
 *   altitude_km     : target altitude above the body's surface (km).
 *                     Must be > 0; use SPODY_EVENT_KIND_IMPACT for
 *                     surface impacts.
 *   action          : LOG (default for monitoring), STOP, or LOG_AND_STOP.
 *
 * Refinement defaults to ON; flip `ev.refined = 0` after construction
 * to opt out. */
SpodyEvent spody_event_altitude_crossing(int naif_id, double body_radius_km,
                                          double altitude_km,
                                          spody_event_action action);

/* Same as above but with a fixed reference point in the integrator's
 * working frame -- used by CR3BP where the body sits at a static
 * synodic position and the ephemeris is absent. `naif_id` is carried
 * for logging only. */
SpodyEvent spody_event_altitude_crossing_at_point(int naif_id,
                                                   const double ref_point[3],
                                                   double body_radius_km,
                                                   double altitude_km,
                                                   spody_event_action action);

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
