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
#include <stdlib.h>
#include <string.h>
#include "spody_forcemodels.h"

/* ============================================================
 * Spacecraft init
 * ============================================================ */
void spody_init_Spacecraft(Spacecraft *sc) {
    if (!sc) return;
    if (sc->mass > 0.0) {
        sc->am_drag = sc->area_drag / sc->mass;
        sc->am_srp  = sc->area_srp  / sc->mass;
    } else {
        sc->am_drag = 0.0;
        sc->am_srp  = 0.0;
    }
}

/* ============================================================
 * Atomic force: two-body
 * ============================================================ */
void spody_force_twobody(double mu, const double r[3], double acc[3]) {
    double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double inv_r3 = 1.0 / (r2 * sqrt(r2));
    double k = -mu * inv_r3;
    acc[0] = k * r[0];
    acc[1] = k * r[1];
    acc[2] = k * r[2];
}

/* ============================================================
 * Atomic force: spherical harmonics (disturbing part)
 *
 * Wraps the rotation pipeline ICRF -> body-fixed -> harm eval ->
 * ICRF. The central-body orientation provider (`get_R`) is passed
 * in by the caller so this kernel stays body-agnostic: any central
 * body whose orientation can be expressed as ICRF<->body-fixed
 * rotation matrices is supported. The application registers the
 * right provider for the configured central body
 * (e.g. spody_bf_rotation_moon) when it builds the
 * ForceModelContext.
 *
 * Selector for the body-fixed kernel: defaults to the HPC variant. Set the
 * env var SPODY_HG_NONHPC=1 (any non-zero int) to switch the whole process
 * to the reference (non-HPC) kernel for audit / regression comparison.
 * The env is sampled lazily on the first call and cached in a static
 * function pointer, so the per-call cost in the RHS hot loop is one
 * indirect call (same as previous direct call after inlining is disabled).
 * ============================================================ */
typedef void (*spody_hg_kernel_fn)(HarmonicGravity *hg, double pos[3], double acc_out[3]);

static spody_hg_kernel_fn spody_select_hg_kernel(void) {
    static spody_hg_kernel_fn cached = NULL;
    if (!cached) {
        const char *env = getenv("SPODY_HG_NONHPC");
        cached = (env && atoi(env)) ? spody_get_hgaccbodyfixed
                                    : spody_get_hgaccbodyfixed_hpc;
    }
    return cached;
}

void spody_bf_rotation_moon(const ForceModelContext *ctx, double et,
                             double R_icrf_to_bf[3][3],
                             double R_bf_to_icrf[3][3]) {
    double angles[3];
    spody_get_lunarlibrationangles(ctx->eph, et, angles);
    spody_getrotmatrix_icrf2moonpa(angles[0], angles[1], angles[2], R_icrf_to_bf);
    spody_getrotmatrix_moonpa2icrf(angles[0], angles[1], angles[2], R_bf_to_icrf);
}

void spody_force_sphericalharmonics(const ForceModelContext *ctx,
                                    double et, const double r[3],
                                    double acc[3]) {
    double R_i2bf[3][3];
    double R_bf2i[3][3];
    ctx->get_bf_rotation(ctx, et, R_i2bf, R_bf2i);

    /* r in body-fixed frame */
    double r_bf[3];
    r_bf[0] = R_i2bf[0][0]*r[0] + R_i2bf[0][1]*r[1] + R_i2bf[0][2]*r[2];
    r_bf[1] = R_i2bf[1][0]*r[0] + R_i2bf[1][1]*r[1] + R_i2bf[1][2]*r[2];
    r_bf[2] = R_i2bf[2][0]*r[0] + R_i2bf[2][1]*r[1] + R_i2bf[2][2]*r[2];

    /* harmonic disturbing acc in body-fixed frame */
    double acc_bf[3];
    spody_select_hg_kernel()(ctx->hg, r_bf, acc_bf);

    /* back to ICRF */
    acc[0] = R_bf2i[0][0]*acc_bf[0] + R_bf2i[0][1]*acc_bf[1] + R_bf2i[0][2]*acc_bf[2];
    acc[1] = R_bf2i[1][0]*acc_bf[0] + R_bf2i[1][1]*acc_bf[1] + R_bf2i[1][2]*acc_bf[2];
    acc[2] = R_bf2i[2][0]*acc_bf[0] + R_bf2i[2][1]*acc_bf[1] + R_bf2i[2][2]*acc_bf[2];
}

/* ============================================================
 * Atomic force: third body (Cowell)
 * ============================================================ */
void spody_force_thirdbody_cowell(double mu_3, const double r_3[3],
                                  const double r_sat[3], double acc[3]) {
    /* direct term: -(r_sat - r_3) / |r_sat - r_3|^3 */
    double dx = r_sat[0] - r_3[0];
    double dy = r_sat[1] - r_3[1];
    double dz = r_sat[2] - r_3[2];
    double d2 = dx*dx + dy*dy + dz*dz;
    double inv_d3 = 1.0 / (d2 * sqrt(d2));

    /* indirect term: -r_3 / |r_3|^3 */
    double s2 = r_3[0]*r_3[0] + r_3[1]*r_3[1] + r_3[2]*r_3[2];
    double inv_s3 = 1.0 / (s2 * sqrt(s2));

    acc[0] = -mu_3 * (dx * inv_d3 + r_3[0] * inv_s3);
    acc[1] = -mu_3 * (dy * inv_d3 + r_3[1] * inv_s3);
    acc[2] = -mu_3 * (dz * inv_d3 + r_3[2] * inv_s3);
}

/* ============================================================
 * Atomic force: SRP (cannonball)
 * ============================================================ */
void spody_force_srp(const Spacecraft *sat, double fraction_sunlight,
                     const double r_sat_to_sun[3], double acc[3]) {
    double s2 = r_sat_to_sun[0]*r_sat_to_sun[0]
              + r_sat_to_sun[1]*r_sat_to_sun[1]
              + r_sat_to_sun[2]*r_sat_to_sun[2];
    double inv_s3 = 1.0 / (s2 * sqrt(s2));

    /* The acceleration points from the Sun towards the satellite, which is
     * the reverse direction of r_sat_to_sun. We absorb the minus sign in f. */
    double f = -SOLAR_LUMINOSITY_4PIC * sat->Cr * sat->am_srp * fraction_sunlight * inv_s3;

    acc[0] = f * r_sat_to_sun[0];
    acc[1] = f * r_sat_to_sun[1];
    acc[2] = f * r_sat_to_sun[2];
}

/* ============================================================
 * Atomic force: atmospheric drag - PLACEHOLDER
 * ============================================================ */
void spody_force_drag(const Spacecraft *sat,
                      const double r_sat[3], const double v_sat[3],
                      double acc[3]) {
    (void)sat; (void)r_sat; (void)v_sat;
    /* TBD: needs atmospheric density model + central-body rotation rate
     * to compute air-relative velocity. Returns zero for now. */
    acc[0] = 0.0;
    acc[1] = 0.0;
    acc[2] = 0.0;
}

/* ============================================================
 * Composite default RHS
 *
 * Sum order (from smallest to largest typical magnitude on LLO),
 * to keep the round-off error of the final summation small:
 *   1. SRP                    ~ 5e-11 km/s^2
 *   2. drag (placeholder, 0)  ~ 1e-9  km/s^2 (LEO only; lunar atmos negligible)
 *   3. third bodies           ~ 3e-9  km/s^2 (Sun, Earth wrt Moon)
 *   4. spherical harmonics    ~ 1e-5  km/s^2
 *   5. central two-body       ~ 1.5e-3 km/s^2  (added last)
 * ============================================================ */
int spody_force_rhs_default(double t, const double *y, double *dy, void *user) {
    ForceModelContext *ctx = (ForceModelContext*)user;
    const double *r = y;
    const double *v = y + 3;

    double et = ctx->et0 + t;

    /* perturbation accumulator (not the 2-body) */
    double acc_pert[3] = { 0.0, 0.0, 0.0 };
    double acc_tmp[3];

    /* ---- 1. SRP ------------------------------------------------ */
    if (ctx->enable_srp && ctx->eph) {
        /* Sun position in central-body frame */
        double sun_pos[3];
        spody_get_ephposition(ctx->eph, ctx->naif_central, 10, et, sun_pos);

        double r_sat_to_sun[3];
        r_sat_to_sun[0] = sun_pos[0] - r[0];
        r_sat_to_sun[1] = sun_pos[1] - r[1];
        r_sat_to_sun[2] = sun_pos[2] - r[2];

        /* eclipse fraction (1 in full sunlight) */
        double fraction = 1.0;
        if (ctx->srp_occulter_radius > 0.0 && ctx->sun_radius > 0.0) {
            double r_occ_to_sat[3];
            double r_occ_to_sun[3];
            if (ctx->srp_occulter_naif == ctx->naif_central) {
                r_occ_to_sat[0] = r[0]; r_occ_to_sat[1] = r[1]; r_occ_to_sat[2] = r[2];
                r_occ_to_sun[0] = sun_pos[0]; r_occ_to_sun[1] = sun_pos[1]; r_occ_to_sun[2] = sun_pos[2];
            } else {
                double occ_pos[3];
                spody_get_ephposition(ctx->eph, ctx->naif_central, ctx->srp_occulter_naif, et, occ_pos);
                r_occ_to_sat[0] = r[0] - occ_pos[0];
                r_occ_to_sat[1] = r[1] - occ_pos[1];
                r_occ_to_sat[2] = r[2] - occ_pos[2];
                r_occ_to_sun[0] = sun_pos[0] - occ_pos[0];
                r_occ_to_sun[1] = sun_pos[1] - occ_pos[1];
                r_occ_to_sun[2] = sun_pos[2] - occ_pos[2];
            }
            fraction = spody_get_sateclipsestatus(r_occ_to_sat, r_occ_to_sun, r_sat_to_sun,
                                                  ctx->sun_radius, ctx->srp_occulter_radius);
        }

        spody_force_srp(ctx->sat, fraction, r_sat_to_sun, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 2. drag (placeholder) -------------------------------- */
    if (ctx->enable_drag) {
        spody_force_drag(ctx->sat, r, v, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 3. third bodies -------------------------------------- */
    if (ctx->n_third > 0 && ctx->eph) {
        for (int i = 0; i < ctx->n_third; i++) {
            double r_3[3];
            spody_get_ephposition(ctx->eph, ctx->naif_central, ctx->third_naif[i], et, r_3);
            spody_force_thirdbody_cowell(ctx->third_mu[i], r_3, r, acc_tmp);
            acc_pert[0] += acc_tmp[0];
            acc_pert[1] += acc_tmp[1];
            acc_pert[2] += acc_tmp[2];
        }
    }

    /* ---- 4. spherical harmonics (disturbing) ------------------ */
    if (ctx->hg && ctx->eph) {
        spody_force_sphericalharmonics(ctx, et, r, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 5. central two-body (largest term, summed last) ------ */
    double acc_2body[3];
    spody_force_twobody(ctx->mu_central, r, acc_2body);

    /* state derivative */
    dy[0] = v[0];
    dy[1] = v[1];
    dy[2] = v[2];
    dy[3] = acc_pert[0] + acc_2body[0];
    dy[4] = acc_pert[1] + acc_2body[1];
    dy[5] = acc_pert[2] + acc_2body[2];

    return 0;
}

/* ============================================================
 * CR3BP context init + RHS
 *
 * Circular Restricted 3-Body Problem in the synodic rotating frame,
 * dimensional units (km, km/s, rad/s). The two primaries are fixed
 * on the x-axis at x1 = -(mu2/(mu1+mu2)) * L (bigger primary) and
 * x2 = +(mu1/(mu1+mu2)) * L (smaller primary). The frame rotates
 * with omega = sqrt((mu1+mu2)/L^3) about +z. Equations of motion:
 *
 *   ax = -mu1 (x-x1)/r1^3 - mu2 (x-x2)/r2^3 + omega^2 * x + 2*omega*vy
 *   ay = -mu1  y    /r1^3 - mu2  y    /r2^3 + omega^2 * y - 2*omega*vx
 *   az = -mu1  z    /r1^3 - mu2  z    /r2^3
 *
 * Time is autonomous (t and et0 are not consulted). The RHS reads only
 * the cr3bp_* fields; HF fields may be NULL/zero.
 * ============================================================ */
void spody_init_CR3BPContext(ForceModelContext *ctx) {
    if (!ctx) return;
    double mu1 = ctx->cr3bp_mu1;
    double mu2 = ctx->cr3bp_mu2;
    double L   = ctx->cr3bp_L;
    if (mu1 <= 0.0 || mu2 <= 0.0 || L <= 0.0) {
        ctx->cr3bp_omega = 0.0;
        ctx->cr3bp_x1    = 0.0;
        ctx->cr3bp_x2    = 0.0;
        return;
    }
    double mu_tot = mu1 + mu2;
    ctx->cr3bp_omega = sqrt(mu_tot / (L * L * L));
    ctx->cr3bp_x1    = -(mu2 / mu_tot) * L;
    ctx->cr3bp_x2    = +(mu1 / mu_tot) * L;
}

int spody_force_rhs_cr3bp(double t, const double *y, double *dy, void *user) {
    (void)t;
    ForceModelContext *ctx = (ForceModelContext*)user;

    double mu1   = ctx->cr3bp_mu1;
    double mu2   = ctx->cr3bp_mu2;
    double omega = ctx->cr3bp_omega;
    double x1    = ctx->cr3bp_x1;
    double x2    = ctx->cr3bp_x2;

    double rx = y[0], ry = y[1], rz = y[2];
    double vx = y[3], vy = y[4];

    /* relative to bigger primary at (x1, 0, 0) */
    double dx1 = rx - x1;
    double r1_sq = dx1*dx1 + ry*ry + rz*rz;
    double inv_r1_3 = 1.0 / (r1_sq * sqrt(r1_sq));

    /* relative to smaller primary at (x2, 0, 0) */
    double dx2 = rx - x2;
    double r2_sq = dx2*dx2 + ry*ry + rz*rz;
    double inv_r2_3 = 1.0 / (r2_sq * sqrt(r2_sq));

    double g1 = -mu1 * inv_r1_3;
    double g2 = -mu2 * inv_r2_3;
    double omega2 = omega * omega;

    dy[0] = y[3];
    dy[1] = y[4];
    dy[2] = y[5];
    dy[3] = g1 * dx1 + g2 * dx2 + omega2 * rx + 2.0 * omega * vy;
    dy[4] = g1 * ry  + g2 * ry  + omega2 * ry - 2.0 * omega * vx;
    dy[5] = g1 * rz  + g2 * rz;

    return 0;
}

/* ============================================================
 * Diagnostic: force breakdown (post-step)
 *
 * Calls each spody_force_* exactly as spody_force_rhs_default does,
 * but stores the per-force contribution into the ForceBreakdown
 * struct. acc_total uses the same sum order as the default RHS so
 * the two are bit-equivalent at the same (t, y).
 * ============================================================ */
void spody_force_breakdown(const ForceModelContext *ctx,
                           double t, const double *y,
                           ForceBreakdown *bd) {
    if (!ctx || !y || !bd) return;
    const double *r = y;
    const double *v = y + 3;

    double et = ctx->et0 + t;

    /* zero everything (covers also the unused part of acc_thirdbody[]) */
    memset(bd, 0, sizeof(*bd));
    bd->t  = t;
    bd->eclipse_fraction = 1.0;

    /* central two-body */
    spody_force_twobody(ctx->mu_central, r, bd->acc_2body);

    /* spherical harmonics */
    if (ctx->hg && ctx->eph) {
        spody_force_sphericalharmonics(ctx, et, r,
                                       bd->acc_sphericalharmonics);
    }

    /* third bodies (per-body + total) */
    if (ctx->n_third > 0 && ctx->eph) {
        int n = ctx->n_third;
        if (n > SPODY_FM_MAX_THIRD) n = SPODY_FM_MAX_THIRD;
        bd->n_third = n;
        for (int i = 0; i < n; i++) {
            double r_3[3];
            spody_get_ephposition(ctx->eph, ctx->naif_central,
                                  ctx->third_naif[i], et, r_3);
            spody_force_thirdbody_cowell(ctx->third_mu[i], r_3, r,
                                         bd->acc_thirdbody[i]);
            bd->acc_thirdbody_total[0] += bd->acc_thirdbody[i][0];
            bd->acc_thirdbody_total[1] += bd->acc_thirdbody[i][1];
            bd->acc_thirdbody_total[2] += bd->acc_thirdbody[i][2];
        }
    }

    /* SRP (with eclipse fraction) */
    if (ctx->enable_srp && ctx->eph) {
        double sun_pos[3];
        spody_get_ephposition(ctx->eph, ctx->naif_central, 10, et, sun_pos);

        double r_sat_to_sun[3];
        r_sat_to_sun[0] = sun_pos[0] - r[0];
        r_sat_to_sun[1] = sun_pos[1] - r[1];
        r_sat_to_sun[2] = sun_pos[2] - r[2];

        double fraction = 1.0;
        if (ctx->srp_occulter_radius > 0.0 && ctx->sun_radius > 0.0) {
            double r_occ_to_sat[3];
            double r_occ_to_sun[3];
            if (ctx->srp_occulter_naif == ctx->naif_central) {
                r_occ_to_sat[0] = r[0]; r_occ_to_sat[1] = r[1]; r_occ_to_sat[2] = r[2];
                r_occ_to_sun[0] = sun_pos[0]; r_occ_to_sun[1] = sun_pos[1]; r_occ_to_sun[2] = sun_pos[2];
            } else {
                double occ_pos[3];
                spody_get_ephposition(ctx->eph, ctx->naif_central,
                                      ctx->srp_occulter_naif, et, occ_pos);
                r_occ_to_sat[0] = r[0] - occ_pos[0];
                r_occ_to_sat[1] = r[1] - occ_pos[1];
                r_occ_to_sat[2] = r[2] - occ_pos[2];
                r_occ_to_sun[0] = sun_pos[0] - occ_pos[0];
                r_occ_to_sun[1] = sun_pos[1] - occ_pos[1];
                r_occ_to_sun[2] = sun_pos[2] - occ_pos[2];
            }
            fraction = spody_get_sateclipsestatus(r_occ_to_sat, r_occ_to_sun,
                                                  r_sat_to_sun,
                                                  ctx->sun_radius,
                                                  ctx->srp_occulter_radius);
        }
        bd->eclipse_fraction = fraction;
        spody_force_srp(ctx->sat, fraction, r_sat_to_sun, bd->acc_srp);
    }

    /* drag (placeholder; writes zero today) */
    if (ctx->enable_drag) {
        spody_force_drag(ctx->sat, r, v, bd->acc_drag);
    }

    /* total: same sum order as rhs_default
     *   acc_pert = SRP + drag + sum(third) + harmonics
     *   acc      = acc_pert + 2body
     */
    bd->acc_total[0] = bd->acc_srp[0] + bd->acc_drag[0]
                     + bd->acc_thirdbody_total[0]
                     + bd->acc_sphericalharmonics[0]
                     + bd->acc_2body[0];
    bd->acc_total[1] = bd->acc_srp[1] + bd->acc_drag[1]
                     + bd->acc_thirdbody_total[1]
                     + bd->acc_sphericalharmonics[1]
                     + bd->acc_2body[1];
    bd->acc_total[2] = bd->acc_srp[2] + bd->acc_drag[2]
                     + bd->acc_thirdbody_total[2]
                     + bd->acc_sphericalharmonics[2]
                     + bd->acc_2body[2];
}
