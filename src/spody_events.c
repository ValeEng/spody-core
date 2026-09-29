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
#include <string.h>
#include "spody_events.h"
#include "spody_solver.h"
#include "spody_eclipse.h"   /* spody_get_satlitfraction     */
#include "spody_const.h"     /* SUN_RADIUS, SUN_NAIF         */

SpodyEvent spody_event_impact(int naif_id, double radius_km, spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind       = SPODY_EVENT_KIND_IMPACT;
    ev.action     = action;
    ev.naif_id    = naif_id;
    ev.radius_km  = radius_km;
    ev.refined    = 1;
    ev.prev_valid = 0;
    return ev;
}

SpodyEvent spody_event_impact_at_point(int naif_id, const double ref_point[3],
                                        double radius_km, spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind          = SPODY_EVENT_KIND_IMPACT;
    ev.action        = action;
    ev.naif_id       = naif_id;
    ev.radius_km     = radius_km;
    ev.refined       = 1;
    ev.has_ref_point = 1;
    if (ref_point) {
        ev.ref_point[0] = ref_point[0];
        ev.ref_point[1] = ref_point[1];
        ev.ref_point[2] = ref_point[2];
    }
    ev.prev_valid = 0;
    return ev;
}

SpodyEvent spody_event_eclipse(int occulter_naif_id, double occulter_radius_km,
                               double threshold_fraction, spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind               = SPODY_EVENT_KIND_ECLIPSE;
    ev.action             = action;
    ev.naif_id            = occulter_naif_id;
    ev.radius_km          = occulter_radius_km;
    ev.threshold_fraction = threshold_fraction;
    ev.refined            = 1;
    ev.prev_valid         = 0;
    return ev;
}

SpodyEvent spody_event_altitude_crossing(int naif_id, double body_radius_km,
                                          double altitude_km,
                                          spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind         = SPODY_EVENT_KIND_ALT_CROSSING;
    ev.action       = action;
    ev.naif_id      = naif_id;
    ev.radius_km    = body_radius_km;
    ev.altitude_km  = altitude_km;
    ev.refined      = 1;
    ev.prev_valid   = 0;
    return ev;
}

SpodyEvent spody_event_altitude_crossing_at_point(int naif_id,
                                                   const double ref_point[3],
                                                   double body_radius_km,
                                                   double altitude_km,
                                                   spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind          = SPODY_EVENT_KIND_ALT_CROSSING;
    ev.action        = action;
    ev.naif_id       = naif_id;
    ev.radius_km     = body_radius_km;
    ev.altitude_km   = altitude_km;
    ev.refined       = 1;
    ev.has_ref_point = 1;
    if (ref_point) {
        ev.ref_point[0] = ref_point[0];
        ev.ref_point[1] = ref_point[1];
        ev.ref_point[2] = ref_point[2];
    }
    ev.prev_valid = 0;
    return ev;
}

/* Distance against the body's shape (see spody_event_body_distance).
 * With `bounded`, for the IMPACT / ALT_CROSSING predicates, a spheroid
 * is resolved only when it can matter: the geodetic altitude h of a
 * point at distance d lies between d - r_eq and d - r_pol, so when the
 * target altitude (0 for IMPACT) is outside that band the side is
 * already known and a bound with the right sign is returned instead
 * of the Bowring iteration. The predicate is evaluated at every step
 * against the always-on IMPACT events; this keeps it to a sqrt away
 * from the ground. The exact value is used wherever the number itself
 * is logged. */
static double body_distance(const SpodyEvent *ev,
                            const ForceModelContext *ctx,
                            double t, const double *y, int bounded)
{
    /* Satellite relative to the body in the integrator's working
     * frame. Three paths: (i) explicit fixed reference point set by
     * the caller (CR3BP), (ii) body is the central body (HF: origin),
     * (iii) body is some other body and its position comes from the
     * ephemeris (HF: third bodies). */
    double r_rel[3] = { y[0], y[1], y[2] };
    if (ev->has_ref_point) {
        for (int i = 0; i < 3; i++) r_rel[i] = y[i] - ev->ref_point[i];
    } else if (ev->naif_id != ctx->naif_central) {
        double body_pos[3] = {0.0, 0.0, 0.0};
        if (ctx->eph) {
            double et = ctx->et0 + t;
            spody_get_ephposition(ctx->eph, ctx->naif_central, ev->naif_id, et, body_pos);
        }
        for (int i = 0; i < 3; i++) r_rel[i] = y[i] - body_pos[i];
    }

    SpodyBodyShape shape = { ev->radius_km, ev->polar_radius_km,
                             { ev->pole[0], ev->pole[1], ev->pole[2] } };
    if (bounded && shape.r_pol > 0.0 && shape.r_pol != shape.r_eq) {
        double h_target = (ev->kind == SPODY_EVENT_KIND_ALT_CROSSING)
                          ? ev->altitude_km : 0.0;
        double d = sqrt(r_rel[0] * r_rel[0] + r_rel[1] * r_rel[1]
                        + r_rel[2] * r_rel[2]);
        if (d - shape.r_eq  > h_target) return d;                          /* surely above */
        if (d - shape.r_pol < h_target) return d - shape.r_pol + shape.r_eq; /* surely below */
    }
    return spody_body_shape_distance(&shape, r_rel);
}

double spody_event_body_distance(const SpodyEvent *ev,
                                 const ForceModelContext *ctx,
                                 double t, const double *y)
{
    return body_distance(ev, ctx, t, y, 0);
}

/* Satellite -> Sun and satellite -> occulter at state (t, y), in the
 * central frame (the occulter offset is zero when it is the central
 * body). An eclipse event is deliberately per-occulter ("eclipsed by
 * the Moon" and "eclipsed by the Earth" are two different events with
 * their own threshold), unlike the SRP force, which combines every
 * occulter into one lit fraction. Returns 0 without an ephemeris. */
static int eclipse_vectors(const SpodyEvent *ev,
                           const ForceModelContext *ctx,
                           double t, const double *y,
                           double sat2sun[3], double sat2occ[1][3])
{
    if (!ctx->eph) return 0;
    double et = ctx->et0 + t;

    double occ_pos_central[3] = { 0.0, 0.0, 0.0 };
    if (ev->naif_id != ctx->naif_central) {
        spody_get_ephposition(ctx->eph, ctx->naif_central, ev->naif_id,
                              et, occ_pos_central);
    }
    double sun_pos_central[3] = { 0.0, 0.0, 0.0 };
    spody_get_ephposition(ctx->eph, ctx->naif_central, SUN_NAIF,
                          et, sun_pos_central);

    for (int i = 0; i < 3; i++) {
        sat2occ[0][i] = occ_pos_central[i]  - y[i];
        sat2sun[i]    = sun_pos_central[i]  - y[i];
    }
    return 1;
}

/* Signed eclipse predicate at state (t, y): > 0 when the satellite is
 * more lit than the event threshold (spody_get_eclipse_residual).
 * Without an ephemeris the satellite counts as fully lit. */
static double eclipse_signed(const SpodyEvent *ev,
                             const ForceModelContext *ctx,
                             double t, const double *y)
{
    double sat2sun[3], sat2occ[1][3];
    if (!eclipse_vectors(ev, ctx, t, y, sat2sun, sat2occ)) return 1.0;
    SpodyBodyShape shape = { ev->radius_km, ev->polar_radius_km,
                             { ev->pole[0], ev->pole[1], ev->pole[2] } };
    return spody_get_eclipse_residual(sat2sun, SUN_RADIUS, sat2occ[0],
                                      &shape, ev->threshold_fraction);
}

int spody_event_check(SpodyEvent *ev,
                      const ForceModelContext *ctx,
                      double t, const double *y)
{
    if (!ev || !ctx || !y) return 0;

    switch (ev->kind) {
        case SPODY_EVENT_KIND_IMPACT: {
            /* IMPACT is one-shot: once latched, the predicate is not
             * re-evaluated. Output fields stay as they were on the
             * first fire (caller has already consumed them). */
            if (ev->triggered) break;

            double d  = body_distance(ev, ctx, t, y, 1);
            if (d >= ev->radius_km) break;

            ev->triggered = 1;
            ev->t_trigger = t;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = y[i];
            ev->distance_at_trigger = spody_event_body_distance(ev, ctx, t, y);
            return 1;
        }
        case SPODY_EVENT_KIND_ECLIPSE:
        case SPODY_EVENT_KIND_ALT_CROSSING: {
            /* Recurring kinds: correctly detecting transitions requires
             * sign tracking across two accepted steps, which only the
             * refined path implements. The coarse path cannot
             * disambiguate "currently inside the predicate" from "just
             * crossed it", so it does not fire. Use
             * spody_event_check_refined with a SPODY_INTEG_RK45
             * integrator instead. */
            break;
        }
        default:
            break;   /* unknown kind: never fires */
    }
    return 0;
}

/* ============================================================
 * Refined check: dense output + Brent root finding
 * ============================================================ */

/* Closure passed to the Brent solver. Generic over event kinds: each
 * kind plugs in its own residual function (impact_residual,
 * eclipse_residual, ...) that reads ev/ctx/integ through the closure
 * and writes the dense-evaluated state into y_buf as scratch. */
typedef struct {
    SpodyEvent              *ev;
    const ForceModelContext *ctx;
    const IntegratorAllData *integ;
    double y_buf[6];   /* scratch for the dense state (6-dim r, v) */
} EventClosure;

static double impact_residual(double theta, void *args) {
    EventClosure *c = (EventClosure*)args;
    /* state at theta on the just-completed step */
    /* time at theta, and the state there */
    double t_theta = c->integ->t_old + theta * c->integ->h_old;
    spody_dense_state_rv6(c->integ, t_theta, c->y_buf);

    /* distance to the body at that state */
    double d  = body_distance(c->ev, c->ctx, t_theta, c->y_buf, 1);
    return d - c->ev->radius_km;
}

static double eclipse_residual(double theta, void *args) {
    EventClosure *c = (EventClosure*)args;
    double t_theta = c->integ->t_old + theta * c->integ->h_old;
    spody_dense_state_rv6(c->integ, t_theta, c->y_buf);
    return eclipse_signed(c->ev, c->ctx, t_theta, c->y_buf);
}

static double alt_crossing_residual(double theta, void *args) {
    EventClosure *c = (EventClosure*)args;
    double t_theta = c->integ->t_old + theta * c->integ->h_old;
    spody_dense_state_rv6(c->integ, t_theta, c->y_buf);
    double d  = body_distance(c->ev, c->ctx, t_theta, c->y_buf, 1);
    return d - c->ev->radius_km - c->ev->altitude_km;
}

int spody_event_check_refined(SpodyEvent *ev,
                              const ForceModelContext *ctx,
                              const IntegratorAllData *integ)
{
    if (!ev || !ctx || !integ) return 0;

    /* dense output is only available for RK45 today -- everything else
     * falls back to the coarse check (which also handles the per-kind
     * latch). */
    if (integ->method != SPODY_INTEG_RK45) {
        return spody_event_check(ev, ctx, integ->t, integ->y);
    }

    switch (ev->kind) {
        case SPODY_EVENT_KIND_IMPACT: {
            /* IMPACT is one-shot: once latched, no further refinement. */
            if (ev->triggered) break;

            /* signed distance at the two ends of the just-completed
             * step. On the very first call (no cache yet) f_start is
             * computed from integ->y_old; on subsequent calls it comes
             * from the value cached on the previous call. */
            double f_end   = body_distance(ev, ctx, integ->t, integ->y, 1) - ev->radius_km;
            double f_start = ev->prev_valid ? ev->prev_distance_signed
                : (body_distance(ev, ctx, integ->t_old, integ->y_old, 1) - ev->radius_km);
            ev->prev_distance_signed = f_end;
            ev->prev_valid = 1;

            /* No sign change -> surface not crossed in [t_old, t]. */
            if ((f_start > 0.0) == (f_end > 0.0)) break;

            /* Bracket Brent on theta in [0, 1] using the closure that
             * evaluates the dense state + body distance at each probe. */
            EventClosure cl;
            cl.ev    = ev;
            cl.ctx   = ctx;
            cl.integ = integ;

            double theta_root = 0.0;
            int rc = spody_solver_brent(impact_residual, &cl,
                                        /*x_lo=*/0.0, /*x_hi=*/1.0,
                                        /*f_lo=*/f_start, /*f_hi=*/f_end,
                                        /*use_provided=*/1,
                                        /*tol=*/1e-12, /*max_iter=*/60,
                                        &theta_root);
            if (rc != SPODY_SOLVER_OK) {
                /* solver failed (shouldn't happen with a valid
                 * bracket): fall back to declaring the trigger at
                 * the end of the step. */
                theta_root = 1.0;
            }

            /* Explicit evaluation at theta_root: Brent's last probe is
             * not guaranteed to be exactly at the converged root, so
             * we re-evaluate the dense state here to make sure y_trigger
             * matches t_trigger. */
            double t_trigger = integ->t_old + theta_root * integ->h_old;
            spody_dense_state_rv6(integ, t_trigger, cl.y_buf);
            double d_trig    = spody_event_body_distance(ev, ctx, t_trigger, cl.y_buf);

            ev->triggered = 1;
            ev->t_trigger = t_trigger;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = cl.y_buf[i];
            ev->distance_at_trigger = d_trig;
            return 1;
        }
        case SPODY_EVENT_KIND_ECLIPSE: {
            /* ECLIPSE is recurring: NO latch on ev->triggered. Every
             * threshold crossing fires a fresh trigger; the caller is
             * expected to consume the output fields (via emit_event)
             * before the next call overwrites them. */

            /* signed predicate at the two ends of the just-completed step */
            double f_end   = eclipse_signed(ev, ctx, integ->t, integ->y);
            double f_start = ev->prev_valid ? ev->prev_distance_signed
                : eclipse_signed(ev, ctx, integ->t_old, integ->y_old);
            ev->prev_distance_signed = f_end;
            ev->prev_valid = 1;

            /* No sign change -> no threshold crossing in [t_old, t]. */
            if ((f_start > 0.0) == (f_end > 0.0)) break;

            EventClosure cl;
            cl.ev    = ev;
            cl.ctx   = ctx;
            cl.integ = integ;

            double theta_root = 0.0;
            int rc = spody_solver_brent(eclipse_residual, &cl,
                                        /*x_lo=*/0.0, /*x_hi=*/1.0,
                                        /*f_lo=*/f_start, /*f_hi=*/f_end,
                                        /*use_provided=*/1,
                                        /*tol=*/1e-12, /*max_iter=*/60,
                                        &theta_root);
            if (rc != SPODY_SOLVER_OK) {
                theta_root = 1.0;
            }

            double t_trigger = integ->t_old + theta_root * integ->h_old;
            spody_dense_state_rv6(integ, t_trigger, cl.y_buf);
            double sat2sun[3], sat2occ[1][3];
            SpodyBodyShape shape = { ev->radius_km, ev->polar_radius_km,
                                     { ev->pole[0], ev->pole[1], ev->pole[2] } };
            double frac_trig = eclipse_vectors(ev, ctx, t_trigger, cl.y_buf, sat2sun, sat2occ)
                ? spody_get_satlitfraction(sat2sun, SUN_RADIUS, sat2occ, &shape, 1)
                : 1.0;

            ev->triggered = 1;
            ev->t_trigger = t_trigger;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = cl.y_buf[i];
            ev->distance_at_trigger = frac_trig;   /* semantic reuse: fraction in [0,1] */
            return 1;
        }
        case SPODY_EVENT_KIND_ALT_CROSSING: {
            /* ALT_CROSSING is recurring: NO latch on ev->triggered.
             * Every sign change of (|r_sat - r_body| - body_radius -
             * altitude) fires, so both ascending and descending
             * crossings are logged. The caller is expected to consume
             * the output fields before the next call overwrites them.
             *
             * `ev->refined == 0` opts out of the Brent root finding:
             * we still do sign tracking but skip the localisation,
             * firing at the END of the accepted step (precision =
             * step size). The toggle is essentially free in steady
             * state -- Brent runs only at the actual crossing step --
             * but it's exposed for users with many altitude bands. */
            double threshold = ev->radius_km + ev->altitude_km;
            double f_end   = body_distance(ev, ctx, integ->t, integ->y, 1)
                             - threshold;
            double f_start = ev->prev_valid ? ev->prev_distance_signed
                : (body_distance(ev, ctx, integ->t_old, integ->y_old, 1)
                   - threshold);
            ev->prev_distance_signed = f_end;
            ev->prev_valid = 1;

            /* No sign change -> altitude not crossed in [t_old, t]. */
            if ((f_start > 0.0) == (f_end > 0.0)) break;

            double theta_root = 1.0;
            EventClosure cl;
            cl.ev    = ev;
            cl.ctx   = ctx;
            cl.integ = integ;
            if (ev->refined) {
                int rc = spody_solver_brent(alt_crossing_residual, &cl,
                                             /*x_lo=*/0.0, /*x_hi=*/1.0,
                                             /*f_lo=*/f_start, /*f_hi=*/f_end,
                                             /*use_provided=*/1,
                                             /*tol=*/1e-12, /*max_iter=*/60,
                                             &theta_root);
                if (rc != SPODY_SOLVER_OK) theta_root = 1.0;
            }

            double t_trigger = integ->t_old + theta_root * integ->h_old;
            spody_dense_state_rv6(integ, t_trigger, cl.y_buf);
            double d_trig    = spody_event_body_distance(ev, ctx, t_trigger, cl.y_buf);

            ev->triggered = 1;
            ev->t_trigger = t_trigger;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = cl.y_buf[i];
            ev->distance_at_trigger = d_trig;
            return 1;
        }
        default:
            break;   /* unknown kind: never fires */
    }
    return 0;
}
