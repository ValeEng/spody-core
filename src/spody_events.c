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

SpodyEvent spody_event_impact(int naif_id, double radius_km, spody_event_action action) {
    SpodyEvent ev;
    memset(&ev, 0, sizeof(ev));
    ev.kind       = SPODY_EVENT_KIND_IMPACT;
    ev.action     = action;
    ev.naif_id    = naif_id;
    ev.radius_km  = radius_km;
    ev.prev_valid = 0;
    return ev;
}

/* Compute |r_sat - r_body|^2 in the central-body inertial frame. */
static double impact_distance2(const SpodyEvent *ev,
                               const ForceModelContext *ctx,
                               double t, const double *y)
{
    if (ev->naif_id == ctx->naif_central) {
        /* satellite position is already in the central frame */
        return y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
    }

    /* body other than the central: query its position in the central frame */
    double body_pos[3] = {0.0, 0.0, 0.0};
    if (ctx->eph) {
        double et = ctx->et0 + t;
        spody_get_ephposition(ctx->eph, ctx->naif_central, ev->naif_id, et, body_pos);
    }
    double dx = y[0] - body_pos[0];
    double dy = y[1] - body_pos[1];
    double dz = y[2] - body_pos[2];
    return dx*dx + dy*dy + dz*dz;
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

            double d2 = impact_distance2(ev, ctx, t, y);
            if (d2 >= ev->radius_km * ev->radius_km) break;

            ev->triggered = 1;
            ev->t_trigger = t;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = y[i];
            ev->distance_at_trigger = sqrt(d2);
            return 1;
        }
        default:
            break;   /* unknown kind: never fires */
    }
    return 0;
}

/* ============================================================
 * Refined check: dense output + Brent root finding
 * ============================================================ */

/* Closure passed to the Brent solver: evaluate the predicate
 * (distance - radius) at a given theta in [0, 1]. */
typedef struct {
    SpodyEvent              *ev;
    const ForceModelContext *ctx;
    const IntegratorAllData *integ;
    double y_buf[6];   /* scratch for dense_eval (caller's frame is 6-dim) */
} ImpactClosure;

static double impact_residual(double theta, void *args) {
    ImpactClosure *c = (ImpactClosure*)args;
    /* state at theta on the just-completed step */
    spody_dense_eval(c->integ, theta, c->y_buf);

    /* time at theta */
    double t_theta = c->integ->t_old + theta * c->integ->h_old;

    /* distance to the body at that state */
    double d2 = impact_distance2(c->ev, c->ctx, t_theta, c->y_buf);
    return sqrt(d2) - c->ev->radius_km;
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
            double f_end   = sqrt(impact_distance2(ev, ctx, integ->t, integ->y)) - ev->radius_km;
            double f_start = ev->prev_valid ? ev->prev_distance_signed
                : (sqrt(impact_distance2(ev, ctx, integ->t_old, integ->y_old)) - ev->radius_km);
            ev->prev_distance_signed = f_end;
            ev->prev_valid = 1;

            /* No sign change -> surface not crossed in [t_old, t]. */
            if ((f_start > 0.0) == (f_end > 0.0)) break;

            /* Bracket Brent on theta in [0, 1] using the closure that
             * evaluates Hermite + impact_distance2 at each probe. */
            ImpactClosure cl;
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
             * we re-evaluate Hermite here to make sure y_trigger
             * matches t_trigger. */
            spody_dense_eval(integ, theta_root, cl.y_buf);
            double t_trigger = integ->t_old + theta_root * integ->h_old;
            double d2_trig   = impact_distance2(ev, ctx, t_trigger, cl.y_buf);

            ev->triggered = 1;
            ev->t_trigger = t_trigger;
            for (int i = 0; i < 6; i++) ev->y_trigger[i] = cl.y_buf[i];
            ev->distance_at_trigger = sqrt(d2_trig);
            return 1;
        }
        default:
            break;   /* unknown kind: never fires */
    }
    return 0;
}
