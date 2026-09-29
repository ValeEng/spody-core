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
 * Kernel for the body-fixed evaluation: the HPC variant by default, the
 * reference (non-HPC) one when SPODY_HG_NONHPC=1 (any non-zero int) for
 * audit / regression comparison. The choice is made once, when the field
 * is loaded (HarmonicGravityData.use_reference_kernel), so the RHS hot
 * loop only reads a flag shared read-only by every worker thread.
 * ============================================================ */

void spody_bf_rotation_moon(const ForceModelContext *ctx, double et,
                             double R_icrf_to_bf[3][3],
                             double R_bf_to_icrf[3][3]) {
    double angles[3];
    spody_get_lunarlibrationangles(ctx->eph, et, angles);
    spody_getrotmatrix_icrf2moonpa(angles[0], angles[1], angles[2], R_icrf_to_bf);
    spody_getrotmatrix_moonpa2icrf(angles[0], angles[1], angles[2], R_bf_to_icrf);
}

void spody_bf_angular_velocity_icrf(const ForceModelContext *ctx, double et,
                                    double omega_icrf[3]) {
    double R_i2bf[3][3], R_bf2i[3][3];
    if (ctx->naif_central == EARTH_NAIF) {
        ctx->get_bf_rotation(ctx, et, R_i2bf, R_bf2i);
        for (int i = 0; i < 3; i++)
            omega_icrf[i] = EARTH_ROT_RATE_RADPS * R_bf2i[i][2];
        return;
    }
    /* R = R_bf2i(t). Its derivative satisfies dR/dt = [omega]x R, so
     * [omega]x = dR/dt R^T; the skew part is averaged over its two
     * mirrored entries. */
    const double h = SPODY_BF_OMEGA_FD_STEP_S;
    double Rp[3][3], Rm[3][3];
    ctx->get_bf_rotation(ctx, et + h, R_i2bf, Rp);
    ctx->get_bf_rotation(ctx, et - h, R_i2bf, Rm);
    ctx->get_bf_rotation(ctx, et,     R_i2bf, R_bf2i);
    double W[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++)
                s += (Rp[i][k] - Rm[i][k]) / (2.0 * h) * R_bf2i[j][k];
            W[i][j] = s;
        }
    omega_icrf[0] = 0.5 * (W[2][1] - W[1][2]);
    omega_icrf[1] = 0.5 * (W[0][2] - W[2][0]);
    omega_icrf[2] = 0.5 * (W[1][0] - W[0][1]);
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
    if (ctx->hg->hgd->use_reference_kernel)
        spody_get_hgaccbodyfixed(ctx->hg, r_bf, acc_bf);
    else
        spody_get_hgaccbodyfixed_hpc(ctx->hg, r_bf, acc_bf);

    /* back to ICRF */
    acc[0] = R_bf2i[0][0]*acc_bf[0] + R_bf2i[0][1]*acc_bf[1] + R_bf2i[0][2]*acc_bf[2];
    acc[1] = R_bf2i[1][0]*acc_bf[0] + R_bf2i[1][1]*acc_bf[1] + R_bf2i[1][2]*acc_bf[2];
    acc[2] = R_bf2i[2][0]*acc_bf[0] + R_bf2i[2][1]*acc_bf[1] + R_bf2i[2][2]*acc_bf[2];
}

/* ============================================================
 * Atomic force: solid-body tide (IERS 2010 sec. 6.2.1, step 1)
 * ============================================================ */

/* Highest degree the tide corrections reach: 3 directly, 4 through
 * the k^(+) coupling. The V/W recursion below needs one more. */
#define TIDE_NMAX 4

/* Acceleration of a field made only of the (unnormalized) coefficients
 * C, S for 2 <= n <= TIDE_NMAX, at the body-fixed position r_bf:
 * the V/W recursion of Montenbruck & Gill (2000) sec. 3.2.4-3.2.5.
 * The static harmonics kernel is not reused: it is built once for a
 * file's coefficient table, while these change at every call. */
static void tide_field_accel(double gm, double r_ref,
                             const double C[TIDE_NMAX + 1][TIDE_NMAX + 1],
                             const double S[TIDE_NMAX + 1][TIDE_NMAX + 1],
                             const double r_bf[3], double acc_bf[3]) {
    double V[TIDE_NMAX + 3][TIDE_NMAX + 3] = {{0.0}};
    double W[TIDE_NMAX + 3][TIDE_NMAX + 3] = {{0.0}};
    double r2  = r_bf[0]*r_bf[0] + r_bf[1]*r_bf[1] + r_bf[2]*r_bf[2];
    double rho = r_ref * r_ref / r2;
    double x0  = r_ref * r_bf[0] / r2;
    double y0  = r_ref * r_bf[1] / r2;
    double z0  = r_ref * r_bf[2] / r2;

    V[0][0] = r_ref / sqrt(r2);
    W[0][0] = 0.0;
    for (int m = 0; m <= TIDE_NMAX + 1; ++m) {
        if (m > 0) {
            V[m][m] = (2*m - 1) * (x0 * V[m-1][m-1] - y0 * W[m-1][m-1]);
            W[m][m] = (2*m - 1) * (x0 * W[m-1][m-1] + y0 * V[m-1][m-1]);
        }
        if (m <= TIDE_NMAX) {
            V[m+1][m] = (2*m + 1) * z0 * V[m][m];
            W[m+1][m] = (2*m + 1) * z0 * W[m][m];
        }
        for (int n = m + 2; n <= TIDE_NMAX + 1; ++n) {
            V[n][m] = ((2*n - 1) * z0 * V[n-1][m] - (n + m - 1) * rho * V[n-2][m]) / (n - m);
            W[n][m] = ((2*n - 1) * z0 * W[n-1][m] - (n + m - 1) * rho * W[n-2][m]) / (n - m);
        }
    }

    double ax = 0.0, ay = 0.0, az = 0.0;
    for (int n = 2; n <= TIDE_NMAX; ++n) {
        ax -= C[n][0] * V[n+1][1];
        ay -= C[n][0] * W[n+1][1];
        az += (n + 1) * (-C[n][0] * V[n+1][0]);
        for (int m = 1; m <= n; ++m) {
            double fac = 0.5 * (n - m + 1) * (n - m + 2);
            ax += 0.5 * (-C[n][m] * V[n+1][m+1] - S[n][m] * W[n+1][m+1])
                + fac * ( C[n][m] * V[n+1][m-1] + S[n][m] * W[n+1][m-1]);
            ay += 0.5 * (-C[n][m] * W[n+1][m+1] + S[n][m] * V[n+1][m+1])
                + fac * (-C[n][m] * W[n+1][m-1] + S[n][m] * V[n+1][m-1]);
            az += (n - m + 1) * (-C[n][m] * V[n+1][m] - S[n][m] * W[n+1][m]);
        }
    }
    double k = gm / (r_ref * r_ref);
    acc_bf[0] = k * ax;
    acc_bf[1] = k * ay;
    acc_bf[2] = k * az;
}

void spody_force_solidtides(const ForceModelContext *ctx, double et,
                            const double r[3], double acc[3]) {
    const SpodySolidTides *td = ctx->tides;
    double R_i2bf[3][3], R_bf2i[3][3];
    ctx->get_bf_rotation(ctx, et, R_i2bf, R_bf2i);

    /* normalized corrections, eq. 6.6 (degree n) and 6.7 (k^(+)) */
    double dC[TIDE_NMAX + 1][TIDE_NMAX + 1] = {{0.0}};
    double dS[TIDE_NMAX + 1][TIDE_NMAX + 1] = {{0.0}};
    for (int j = 0; j < td->n_raisers; ++j) {
        double p[3], pb[3];
        spody_get_ephposition(ctx->eph, ctx->naif_central,
                              td->raiser_naif[j], et, p);
        for (int i = 0; i < 3; ++i)
            pb[i] = R_i2bf[i][0]*p[0] + R_i2bf[i][1]*p[1] + R_i2bf[i][2]*p[2];
        double rj  = sqrt(pb[0]*pb[0] + pb[1]*pb[1] + pb[2]*pb[2]);
        double s   = pb[2] / rj;                    /* sin(latitude) */
        double c   = sqrt(pb[0]*pb[0] + pb[1]*pb[1]) / rj;
        double lam = atan2(pb[1], pb[0]);
        double q   = td->r_ref / rj;
        double f2  = (td->raiser_mu[j] / td->gm) * q * q * q;
        double f3  = f2 * q;

        /* fully normalized Legendre functions, no Condon-Shortley phase */
        double P[4][4] = {{0.0}};
        P[2][0] = sqrt(5.0) * 0.5 * (3.0*s*s - 1.0);
        P[2][1] = sqrt(15.0) * s * c;
        P[2][2] = sqrt(15.0) * 0.5 * c * c;
        P[3][0] = sqrt(7.0) * 0.5 * s * (5.0*s*s - 3.0);
        P[3][1] = sqrt(42.0) * 0.25 * c * (5.0*s*s - 1.0);
        P[3][2] = sqrt(105.0) * 0.5 * s * c * c;
        P[3][3] = sqrt(70.0) * 0.25 * c * c * c;

        for (int n = 2; n <= td->max_degree; ++n) {
            double f = (n == 2 ? f2 : f3) / (2*n + 1);
            for (int m = 0; m <= n; ++m) {
                double cm = cos(m * lam), sm = sin(m * lam);
                /* (kr + i ki)(cos - i sin) = dC - i dS */
                dC[n][m] += f * P[n][m] * (td->k_re[n][m]*cm + td->k_im[n][m]*sm);
                dS[n][m] += f * P[n][m] * (td->k_re[n][m]*sm - td->k_im[n][m]*cm);
            }
        }
        for (int m = 0; m <= 2; ++m) {
            double f = f2 / 5.0 * td->kplus[m] * P[2][m];
            dC[4][m] += f * cos(m * lam);
            dS[4][m] += f * sin(m * lam);
        }
    }
    dC[2][0] -= td->dc20_perm;

    /* normalized -> unnormalized: N_nm = sqrt((2-d0m)(2n+1)(n-m)!/(n+m)!) */
    static const double fact[2 * TIDE_NMAX + 1] =
        { 1, 1, 2, 6, 24, 120, 720, 5040, 40320 };
    for (int n = 2; n <= TIDE_NMAX; ++n)
        for (int m = 0; m <= n; ++m) {
            double N = sqrt((m == 0 ? 1.0 : 2.0) * (2*n + 1) * fact[n-m] / fact[n+m]);
            dC[n][m] *= N;
            dS[n][m] *= N;
        }

    double r_bf[3], a_bf[3];
    for (int i = 0; i < 3; ++i)
        r_bf[i] = R_i2bf[i][0]*r[0] + R_i2bf[i][1]*r[1] + R_i2bf[i][2]*r[2];
    tide_field_accel(td->gm, td->r_ref, dC, dS, r_bf, a_bf);
    for (int i = 0; i < 3; ++i)
        acc[i] = R_bf2i[i][0]*a_bf[0] + R_bf2i[i][1]*a_bf[1] + R_bf2i[i][2]*a_bf[2];
}

/* ============================================================
 * Atomic force: Earth radiation pressure (Knocke et al. 1988)
 * ============================================================ */
void spody_force_earthradiation(const ForceModelContext *ctx, double et,
                                const double r[3], double acc[3]) {
    static const double gx[EARTHRAD_N_NADIR] = { -GL6_X3, -GL6_X2, -GL6_X1,
                                                  GL6_X1,  GL6_X2,  GL6_X3 };
    static const double gw[EARTHRAD_N_NADIR] = {  GL6_W3,  GL6_W2,  GL6_W1,
                                                  GL6_W1,  GL6_W2,  GL6_W3 };
    double R  = ctx->R_central;
    double rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    acc[0] = acc[1] = acc[2] = 0.0;
    if (!(rn > R)) return;

    /* Sun direction, and the solar flux at the Earth as F/c in the
     * units spody_force_srp uses (km/s^2 per m^2/kg of A/m) */
    double sun[3];
    spody_get_ephposition(ctx->eph, ctx->naif_central, SUN_NAIF, et, sun);
    double ds = sqrt(sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2]);
    double p_sun = SOLAR_LUMINOSITY_4PIC / (ds * ds);
    for (int i = 0; i < 3; ++i) sun[i] /= ds;

    /* body-fixed z axis in ICRF: sin(latitude) = z_bf . n */
    double R_i2bf[3][3], R_bf2i[3][3];
    ctx->get_bf_rotation(ctx, et, R_i2bf, R_bf2i);
    const double *zbf = R_i2bf[2];

    /* season */
    double wt = 2.0 * PI * (et - KNOCKE_T0_ET) / JULIAN_YEAR_S;
    double a1 = KNOCKE_C0 + KNOCKE_C1 * cos(wt) + KNOCKE_C2 * sin(wt);
    double e1 = KNOCKE_K0 + KNOCKE_K1 * cos(wt) + KNOCKE_K2 * sin(wt);

    /* ray basis: z toward nadir, x and y across it */
    double zh[3] = { -r[0]/rn, -r[1]/rn, -r[2]/rn };
    double hp[3] = { 1.0, 0.0, 0.0 };
    if (fabs(zh[0]) > 0.9) { hp[0] = 0.0; hp[1] = 1.0; }
    double xh[3] = { zh[1]*hp[2] - zh[2]*hp[1],
                     zh[2]*hp[0] - zh[0]*hp[2],
                     zh[0]*hp[1] - zh[1]*hp[0] };
    double xn = sqrt(xh[0]*xh[0] + xh[1]*xh[1] + xh[2]*xh[2]);
    for (int i = 0; i < 3; ++i) xh[i] /= xn;
    double yh[3] = { zh[1]*xh[2] - zh[2]*xh[1],
                     zh[2]*xh[0] - zh[0]*xh[2],
                     zh[0]*xh[1] - zh[1]*xh[0] };

    /* cos(nadir angle) from cos(eta_max) = sqrt(1 - (R/r)^2) to 1 */
    double ce_min = sqrt(1.0 - (R / rn) * (R / rn));
    double half   = 0.5 * (1.0 - ce_min);
    double dlam   = 2.0 * PI / EARTHRAD_N_AZIMUTH;
    double E[3]   = { 0.0, 0.0, 0.0 };   /* sum of (M/F) dOmega (-ray) */
    /* azimuth directions, (j + 1/2) dlam, by rotating the first one:
     * the same for every ring, so computed once per call */
    double caz[EARTHRAD_N_AZIMUTH], saz[EARTHRAD_N_AZIMUTH];
    double cd = cos(dlam), sd = sin(dlam);
    caz[0] = cos(0.5 * dlam);
    saz[0] = sin(0.5 * dlam);
    for (int j = 1; j < EARTHRAD_N_AZIMUTH; ++j) {
        caz[j] = caz[j-1] * cd - saz[j-1] * sd;
        saz[j] = saz[j-1] * cd + caz[j-1] * sd;
    }
    for (int k = 0; k < EARTHRAD_N_NADIR; ++k) {
        double ce = ce_min + half * (gx[k] + 1.0);
        double se = sqrt(1.0 - ce * ce);
        double w  = gw[k] * half * dlam;
        for (int j = 0; j < EARTHRAD_N_AZIMUTH; ++j) {
            double cl = caz[j], sl = saz[j], ray[3], n[3];
            for (int i = 0; i < 3; ++i)
                ray[i] = ce * zh[i] + se * (cl * xh[i] + sl * yh[i]);
            /* nearest intersection of r + t ray with the sphere |p| = R */
            double b = r[0]*ray[0] + r[1]*ray[1] + r[2]*ray[2];
            double disc = b * b - (rn * rn - R * R);
            double t = -b - sqrt(disc > 0.0 ? disc : 0.0);
            for (int i = 0; i < 3; ++i) n[i] = (r[i] + t * ray[i]) / R;
            double sphi = zbf[0]*n[0] + zbf[1]*n[1] + zbf[2]*n[2];
            double p2   = 1.5 * sphi * sphi - 0.5;
            double M    = 0.25 * (KNOCKE_E0 + e1 * sphi + KNOCKE_E2 * p2);
            double cz   = n[0]*sun[0] + n[1]*sun[1] + n[2]*sun[2];
            if (cz > 0.0)
                M += (KNOCKE_A0 + a1 * sphi + KNOCKE_A2 * p2) * cz;
            for (int i = 0; i < 3; ++i) E[i] -= M * w * ray[i];
        }
    }
    double f = ctx->sat->Cr * ctx->sat->am_srp * p_sun / PI;
    for (int i = 0; i < 3; ++i) acc[i] = f * E[i];
}

/* ============================================================
 * Atomic force: general relativity (IERS 2010 eq. 10.12, line 1)
 * ============================================================ */
void spody_force_relativity(double mu, const double r[3], const double v[3],
                            double acc[3]) {
    double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double rn = sqrt(r2);
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double rv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
    double k  = mu / (SPEED_OF_LIGHT_KMS * SPEED_OF_LIGHT_KMS * r2 * rn);
    double cr = 2.0 * (PPN_BETA + PPN_GAMMA) * mu / rn - PPN_GAMMA * v2;
    double cv = 2.0 * (1.0 + PPN_GAMMA) * rv;
    for (int i = 0; i < 3; ++i)
        acc[i] = k * (cr * r[i] + cv * v[i]);
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
 * Atomic force: atmospheric drag
 * ============================================================
 *   a_drag = -0.5 * rho * |v_rel|^2 * Cd * (A/m) * v_rel_hat
 *
 * Density comes from the central body's atmosphere callback
 * (`ctx->atmosphere->density`); air-relative velocity uses the
 * body's BF +Z axis rotated into ICRF as the rotation-axis
 * direction, scaled by `ctx->body_spin_rad_s`. Wrappers convert
 * spacecraft am_drag from m^2/kg + density from kg/m^3 + velocity
 * km/s into the km/s^2 the integrator expects (the trailing /KM2M
 * is the m^2 -> km^2 contraction of A/m times the km^2 -> m^2
 * expansion of |v|^2 -- net 1/KM2M factor).
 */
void spody_force_drag(const ForceModelContext *ctx, double et,
                      const double r_sat[3], const double v_sat[3],
                      double acc[3]) {
    acc[0] = 0.0; acc[1] = 0.0; acc[2] = 0.0;
    if (!ctx || !ctx->sat || !ctx->atmosphere || !ctx->atmosphere->density
            || !ctx->get_bf_rotation || ctx->body_spin_rad_s <= 0.0) {
        return;
    }
    double R_icrf_to_bf[3][3], R_bf_to_icrf[3][3];
    ctx->get_bf_rotation(ctx, et, R_icrf_to_bf, R_bf_to_icrf);

    double r_bf[3];
    r_bf[0] = R_icrf_to_bf[0][0]*r_sat[0] + R_icrf_to_bf[0][1]*r_sat[1] + R_icrf_to_bf[0][2]*r_sat[2];
    r_bf[1] = R_icrf_to_bf[1][0]*r_sat[0] + R_icrf_to_bf[1][1]*r_sat[1] + R_icrf_to_bf[1][2]*r_sat[2];
    r_bf[2] = R_icrf_to_bf[2][0]*r_sat[0] + R_icrf_to_bf[2][1]*r_sat[1] + R_icrf_to_bf[2][2]*r_sat[2];

    double rho_kg_m3 = 0.0;
    if (ctx->atmosphere->density(ctx, et, r_bf, &rho_kg_m3) != 0) return;
    if (!(rho_kg_m3 > 0.0)) return;
    rho_kg_m3 *= spody_interpolate_density_scale(ctx->density_scale, et);

    /* omega vector in ICRF = body's BF +Z rotated into ICRF, scaled. */
    double omega[3];
    omega[0] = ctx->body_spin_rad_s * R_bf_to_icrf[0][2];
    omega[1] = ctx->body_spin_rad_s * R_bf_to_icrf[1][2];
    omega[2] = ctx->body_spin_rad_s * R_bf_to_icrf[2][2];

    double v_rel[3];
    v_rel[0] = v_sat[0] - (omega[1]*r_sat[2] - omega[2]*r_sat[1]);
    v_rel[1] = v_sat[1] - (omega[2]*r_sat[0] - omega[0]*r_sat[2]);
    v_rel[2] = v_sat[2] - (omega[0]*r_sat[1] - omega[1]*r_sat[0]);

    double v_rel_mag2 = v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2];
    double v_rel_mag  = sqrt(v_rel_mag2);
    if (!(v_rel_mag > 0.0)) return;

    /* Mixed units: rho [kg/m^3], am_drag [m^2/kg], v [km/s].
     * 0.5 * rho * v_m^2 * (A/m) gives a [m/s^2]; convert to km/s^2
     * via M2KM. Equivalently 0.5 * rho * v_km^2 * KM2M^2 * (A/m) *
     * M2KM = 0.5 * rho * v_km^2 * KM2M * (A/m). */
    double f = -0.5 * rho_kg_m3 * v_rel_mag * ctx->sat->Cd * ctx->sat->am_drag * KM2M;
    acc[0] = f * v_rel[0];
    acc[1] = f * v_rel[1];
    acc[2] = f * v_rel[2];
}

/* ============================================================
 * Lit fraction seen by the satellite, combining every occulter in
 * ctx->srp_occulter_naif[].
 *
 * Pure plumbing: it pulls the occulter positions out of the
 * ephemeris and hands the geometry to spody_get_satlitfraction. Both
 * the RHS and the diagnostic breakdown go through it so the two
 * cannot drift apart.
 *
 * Costs no extra Chebyshev evaluation: the occulters are a subset of
 * the third bodies, which the same step queries at the same epoch for
 * their gravity, and spody_get_ephposition caches per (target, epoch).
 * The central body needs no query at all -- it sits at the origin of
 * the frame the states are integrated in.
 * ============================================================ */
static double srp_lit_fraction(const ForceModelContext *ctx, double et,
                               const double r[3], const double sat2sun[3]) {
    if (ctx->srp_n_occulters <= 0 || ctx->sun_radius <= 0.0) return 1.0;

    double sat2occ[SPODY_ECL_MAX_OCCULTERS][3];
    SpodyBodyShape shape[SPODY_ECL_MAX_OCCULTERS];
    int    n = 0;

    for (int i = 0; i < ctx->srp_n_occulters; ++i) {
        if (ctx->srp_occulter_shape[i].r_eq <= 0.0) continue;
        double occ_pos[3] = { 0.0, 0.0, 0.0 };
        if (ctx->srp_occulter_naif[i] != ctx->naif_central) {
            spody_get_ephposition(ctx->eph, ctx->naif_central,
                                  ctx->srp_occulter_naif[i], et, occ_pos);
        }
        sat2occ[n][0] = occ_pos[0] - r[0];
        sat2occ[n][1] = occ_pos[1] - r[1];
        sat2occ[n][2] = occ_pos[2] - r[2];
        shape[n]      = ctx->srp_occulter_shape[i];
        n++;
    }
    if (n == 0) return 1.0;

    return spody_get_satlitfraction(sat2sun, ctx->sun_radius,
                                    sat2occ, shape, n);
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
/* The rule itself, kept free of every struct in this file so it can be
 * read, checked and reasoned about as plain arithmetic: the smallest
 * degree whose truncation error stays below a relative accuracy of
 * exp(-ln_inv_eps) at radius r_km, clamped to [2, N_max].
 *
 * At or below the reference sphere the ratio is <= 1 and there is no
 * geometric decay to exploit, so nothing can be dropped. */
static int hgdegree_for_radius(double r_km, double R_ref_km,
                               double ln_inv_eps, int N_max) {
    if (!(r_km > R_ref_km) || !(R_ref_km > 0.0)) return N_max;

    double n = ln_inv_eps / log(r_km / R_ref_km);
    if (!(n > 2.0)) return 2;                  /* also catches NaN */
    if (n >= (double)N_max) return N_max;
    return (int)ceil(n);
}

int spody_adapt_hgdegree(double t, const double *y, double h, void *user) {
    (void)t;
    ForceModelContext *ctx = (ForceModelContext*)user;
    if (!ctx || !ctx->hg || !ctx->hg->hgd) return 0;   /* harmonics off */
    const HarmonicGravityData *hgd = ctx->hg->hgd;

    double r = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    double v = sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);

    /* Lowest radius the step can reach, bounded by the straight-line
     * excursion. Under gravity the path bends inward of that line, so
     * the accuracy margin on the degree absorbs the difference. */
    double r_bound = r - SPODY_HG_ADAPTIVE_STEP_MARGIN * v * fabs(h);

    int n_before = ctx->hg->N_eval;
    ctx->hg->N_eval = hgdegree_for_radius(r_bound, hgd->R_ref,
                                          SPODY_HG_ADAPTIVE_LN_INV_EPS,
                                          hgd->N);
    return ctx->hg->N_eval != n_before;
}

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
        spody_get_ephposition(ctx->eph, ctx->naif_central, SUN_NAIF, et, sun_pos);

        double r_sat_to_sun[3];
        r_sat_to_sun[0] = sun_pos[0] - r[0];
        r_sat_to_sun[1] = sun_pos[1] - r[1];
        r_sat_to_sun[2] = sun_pos[2] - r[2];

        /* eclipse fraction over every occulter (1 in full sunlight) */
        double fraction = srp_lit_fraction(ctx, et, r, r_sat_to_sun);

        spody_force_srp(ctx->sat, fraction, r_sat_to_sun, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 1b. Earth radiation (albedo + infrared) ------------- */
    if (ctx->enable_earthradiation && ctx->eph) {
        spody_force_earthradiation(ctx, et, r, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 2. drag ---------------------------------------------- */
    if (ctx->enable_drag) {
        spody_force_drag(ctx, et, r, v, acc_tmp);
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

    /* ---- 4b. solid-body tide ---------------------------------- */
    if (ctx->tides && ctx->eph) {
        spody_force_solidtides(ctx, et, r, acc_tmp);
        acc_pert[0] += acc_tmp[0];
        acc_pert[1] += acc_tmp[1];
        acc_pert[2] += acc_tmp[2];
    }

    /* ---- 4c. general relativity (Schwarzschild) --------------- */
    if (ctx->enable_relativity) {
        spody_force_relativity(ctx->mu_central, r, v, acc_tmp);
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

void spody_inertial_to_cr3bp_synodic(
        const double r_primary_inertial[3],
        const double v_primary_inertial[3],
        double mu1_km3_s2, double mu2_km3_s2, double L_km,
        int    primary_index,
        double r_synodic[3], double v_synodic[3]) {
    double mu_tot = mu1_km3_s2 + mu2_km3_s2;
    double omega  = sqrt(mu_tot / (L_km * L_km * L_km));
    /* Primary positions on the synodic x-axis (matches the convention
     * used by spody_init_CR3BPContext above): primary_1 at
     * -mu2/mu_tot * L, primary_2 at +mu1/mu_tot * L. At t = 0 these
     * also coincide with the underlying inertial x-axis. */
    double x_primary = (primary_index == 2)
                       ?  (mu1_km3_s2 / mu_tot) * L_km
                       : -(mu2_km3_s2 / mu_tot) * L_km;

    /* Step 1: primary-centered inertial -> barycenter-centered inertial.
     * At t = 0 the primary sits at (x_primary, 0, 0) and moves
     * tangentially with v = omega * x_primary along +y. */
    double r_bary[3] = {
        r_primary_inertial[0] + x_primary,
        r_primary_inertial[1],
        r_primary_inertial[2],
    };
    double v_bary[3] = {
        v_primary_inertial[0],
        v_primary_inertial[1] + omega * x_primary,
        v_primary_inertial[2],
    };

    /* Step 2: barycenter-inertial -> synodic-rotating at t = 0.
     * Rotation is identity; velocity loses omega x r (omega = (0,0,omega)).
     * The sign of the omega-x-r terms below comes from
     * v_rot = v_inert - omega x r, with omega x r = (-w*ry, +w*rx, 0). */
    r_synodic[0] = r_bary[0];
    r_synodic[1] = r_bary[1];
    r_synodic[2] = r_bary[2];
    v_synodic[0] = v_bary[0] + omega * r_synodic[1];
    v_synodic[1] = v_bary[1] - omega * r_synodic[0];
    v_synodic[2] = v_bary[2];
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

    /* solid-body tide */
    if (ctx->tides && ctx->eph) {
        spody_force_solidtides(ctx, et, r, bd->acc_solidtides);
    }

    /* general relativity (Schwarzschild) */
    if (ctx->enable_relativity) {
        spody_force_relativity(ctx->mu_central, r, v, bd->acc_relativity);
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
        spody_get_ephposition(ctx->eph, ctx->naif_central, SUN_NAIF, et, sun_pos);

        double r_sat_to_sun[3];
        r_sat_to_sun[0] = sun_pos[0] - r[0];
        r_sat_to_sun[1] = sun_pos[1] - r[1];
        r_sat_to_sun[2] = sun_pos[2] - r[2];

        double fraction = srp_lit_fraction(ctx, et, r, r_sat_to_sun);
        bd->eclipse_fraction = fraction;
        spody_force_srp(ctx->sat, fraction, r_sat_to_sun, bd->acc_srp);
    }

    /* Earth radiation (albedo + infrared) */
    if (ctx->enable_earthradiation && ctx->eph) {
        spody_force_earthradiation(ctx, et, r, bd->acc_earthradiation);
    }

    /* drag */
    if (ctx->enable_drag) {
        spody_force_drag(ctx, et, r, v, bd->acc_drag);
    }

    /* total: accumulated exactly as rhs_default does -- a perturbation
     * sum starting from zero that takes SRP, drag, each third body in
     * turn, the harmonics, the solid tide and relativity, then the
     * two-body term last; Earth radiation comes right after SRP. Adding the
     * pre-summed acc_thirdbody_total instead regroups the third bodies
     * and, with two or more of them next to a non-zero SRP, rounds
     * differently from the RHS. A disabled force is exactly zero here,
     * and adding zero leaves the sum unchanged. */
    for (int k = 0; k < 3; k++) {
        double acc_pert = 0.0;
        acc_pert += bd->acc_srp[k];
        if (ctx->enable_earthradiation && ctx->eph)
            acc_pert += bd->acc_earthradiation[k];
        acc_pert += bd->acc_drag[k];
        for (int i = 0; i < bd->n_third; i++)
            acc_pert += bd->acc_thirdbody[i][k];
        acc_pert += bd->acc_sphericalharmonics[k];
        if (ctx->tides && ctx->eph)
            acc_pert += bd->acc_solidtides[k];
        if (ctx->enable_relativity)
            acc_pert += bd->acc_relativity[k];
        bd->acc_total[k] = acc_pert + bd->acc_2body[k];
    }
}
