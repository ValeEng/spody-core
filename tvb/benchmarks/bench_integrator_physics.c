/*
 * Benchmark: time per integration step as we add physics to the RHS.
 *
 * Scenarios (low lunar orbit, RKDP45 adaptive):
 *   1. Two-body around the Moon (baseline)
 *   2. + Earth as third body
 *   3. + Sun as third body
 *   4. + SRP (solar radiation pressure)
 *   5. Spherical harmonics N=10
 *   6. Spherical harmonics N=50
 *   7. Spherical harmonics N=100
 *   8. Spherical harmonics N=200
 *   9. Full: harmonics N=100 + Earth + Sun + SRP
 *
 * For each scenario reports:
 *   accepted steps, RHS calls, total wall time, time per step, time per RHS.
 *
 * Usage:
 *   bench_integrator_physics <de440.spody> <grgm1200a.tab>
 *
 * If a data file is missing the dependent scenarios are skipped.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "spody_core.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

// ---- RHS payload ------------------------------------------------------

typedef struct {
    double mu_moon;

    HarmonicGravity *hg;        // NULL = no harmonics
    MappedEphemeris *eph;       // NULL = no ephemeris-based perturbations

    int use_earth_3rd;
    int use_sun_3rd;
    double srp_AoM_Cr;          // 0.0 = SRP disabled, else = Cr * A/m  (m^2/kg)

    double t0_jd;               // JD at integration t = 0
    long n_rhs_calls;
} OrbitParams;

// ---- helpers ----------------------------------------------------------

static double now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

static void mat_mul_vec(const double R[3][3], const double v[3], double out[3]) {
    out[0] = R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2];
    out[1] = R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2];
    out[2] = R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2];
}

// Cowell third-body acceleration (Battin form)
//   a = -mu_3 * ( (r - s)/|r-s|^3 + s/|s|^3 )
// r is the satellite position in moon-centered inertial,
// s is the third body position in moon-centered inertial.
static void cowell_3rd(double mu_3, const double r[3], const double s[3], double acc_out[3]) {
    double dx = r[0] - s[0];
    double dy = r[1] - s[1];
    double dz = r[2] - s[2];
    double d2 = dx*dx + dy*dy + dz*dz;
    double inv_d3 = 1.0 / (d2 * sqrt(d2));
    double s2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
    double inv_s3 = 1.0 / (s2 * sqrt(s2));

    acc_out[0] = -mu_3 * (dx*inv_d3 + s[0]*inv_s3);
    acc_out[1] = -mu_3 * (dy*inv_d3 + s[1]*inv_s3);
    acc_out[2] = -mu_3 * (dz*inv_d3 + s[2]*inv_s3);
}

// ---- the RHS ----------------------------------------------------------

static int orbit_rhs(double t, const double *y, double *dy, void *user) {
    OrbitParams *p = (OrbitParams*)user;
    p->n_rhs_calls++;

    const double *r = y;
    double a[3] = {0.0, 0.0, 0.0};

    // 1. Two-body around the Moon (always present)
    double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double inv_r3 = 1.0 / (r2 * sqrt(r2));
    double k2 = -p->mu_moon * inv_r3;
    a[0] = k2 * r[0];
    a[1] = k2 * r[1];
    a[2] = k2 * r[2];

    double jd = p->t0_jd + t / SECONDSxDAY;

    // 2. Spherical harmonics (disturbing part only -- 2-body already added)
    if (p->hg) {
        double angles[3], R_i2pa[3][3], R_pa2i[3][3];
        spody_get_lunarlibrationangles(p->eph, jd, angles);
        spody_getrotmatrix_icrf2moonpa(angles[0], angles[1], angles[2], R_i2pa);
        spody_getrotmatrix_moonpa2icrf(angles[0], angles[1], angles[2], R_pa2i);

        double r_pa[3], acc_pa[3], acc_i[3];
        mat_mul_vec(R_i2pa, r, r_pa);
        spody_get_hgaccbodyfixed(p->hg, r_pa, acc_pa);
        mat_mul_vec(R_pa2i, acc_pa, acc_i);

        a[0] += acc_i[0];
        a[1] += acc_i[1];
        a[2] += acc_i[2];
    }

    // 3. Third bodies (Earth and Sun)
    if (p->use_earth_3rd && p->eph) {
        double earth_pos[3], acc[3];
        spody_get_ephposition(p->eph, 301, 399, jd, earth_pos);
        cowell_3rd(EARTH_MU, r, earth_pos, acc);
        a[0] += acc[0]; a[1] += acc[1]; a[2] += acc[2];
    }
    if (p->use_sun_3rd && p->eph) {
        double sun_pos[3], acc[3];
        spody_get_ephposition(p->eph, 301, 10, jd, sun_pos);
        cowell_3rd(SUN_MU, r, sun_pos, acc);
        a[0] += acc[0]; a[1] += acc[1]; a[2] += acc[2];
    }

    // 4. SRP (cannonball model). Acceleration pushes away from Sun.
    if (p->srp_AoM_Cr > 0.0 && p->eph) {
        double sun_pos[3];
        spody_get_ephposition(p->eph, 301, 10, jd, sun_pos);
        double sx = r[0] - sun_pos[0];
        double sy = r[1] - sun_pos[1];
        double sz = r[2] - sun_pos[2];
        double s2 = sx*sx + sy*sy + sz*sz;
        double inv_s3 = 1.0 / (s2 * sqrt(s2));
        // f * (sat - sun) points away from the Sun
        double f = SOLAR_LUMINOSITY_4PIC * p->srp_AoM_Cr * inv_s3;
        a[0] += f * sx;
        a[1] += f * sy;
        a[2] += f * sz;
    }

    dy[0] = y[3]; dy[1] = y[4]; dy[2] = y[5];
    dy[3] = a[0]; dy[4] = a[1]; dy[5] = a[2];
    return 0;
}

// ---- scenario runner --------------------------------------------------

static void run_scenario(const char *label,
                         OrbitParams *p,
                         const double *y0,
                         double t_end) {

    p->n_rhs_calls = 0;

    IntegratorOptions opt;
    spody_default_integrator_options(SPODY_INTEG_RK45, &opt);
    opt.rel_tol = 1e-9;
    opt.h_init  = 10.0;
    opt.h_max   = 200.0;

    IntegratorAllData integ;
    if (spody_setup_integrator(&integ, SPODY_INTEG_RK45, &opt, 6, orbit_rhs, p) != SPODY_INTEG_OK) {
        printf("  setup failed\n");
        return;
    }
    spody_set_integrator_state(&integ, 0.0, y0);

    long n_steps = 0;
    double t0 = now_ms();
    while (integ.t < t_end) {
        if (spody_propagate_onestep(&integ) != SPODY_INTEG_OK) {
            printf("  step failed at t=%.3f\n", integ.t);
            spody_free_integrator(&integ);
            return;
        }
        if (integ.t + integ.h > t_end) integ.h = t_end - integ.t;
        n_steps++;
    }
    double dt = now_ms() - t0;

    double us_per_step = (n_steps > 0) ? (dt * 1000.0 / (double)n_steps) : 0.0;
    double us_per_rhs  = (p->n_rhs_calls > 0) ? (dt * 1000.0 / (double)p->n_rhs_calls) : 0.0;

    printf("  %-40s | steps=%6ld | rhs=%7ld | %8.2f ms | %8.2f us/step | %7.2f us/rhs\n",
           label, n_steps, p->n_rhs_calls, dt, us_per_step, us_per_rhs);

    spody_free_integrator(&integ);
}

// ---- main -------------------------------------------------------------

int main(int argc, char **argv) {

    const char *eph_path  = (argc >= 2) ? argv[1] : NULL;
    const char *grav_path = (argc >= 3) ? argv[2] : NULL;

    // ---- load shared data (each per-thread handle bound to it) ----
    MappedEphemerisData med = {0};
    MappedEphemeris    eph_handle = {0};
    int eph_ok = 0;
    if (eph_path) {
        if (spody_setup_MappedEphemerisData(&med, eph_path) == 0 &&
            spody_setup_MappedEphemeris(&eph_handle, &med) == 0) {
            eph_ok = 1;
        } else {
            printf("WARN: ephemeris setup failed for '%s' - perturbation scenarios will be skipped\n", eph_path);
        }
    } else {
        printf("WARN: no ephemeris path provided - perturbation scenarios will be skipped\n");
    }

    // ---- low lunar orbit setup ----
    double mu = MOON_MU;
    double R_moon = 1737.4;             // km
    double altitude = 100.0;            // km
    double r0 = R_moon + altitude;
    double v0 = sqrt(mu / r0);
    double y0[6] = { r0, 0.0, 0.0,  0.0, v0, 0.0 };

    double T = 2.0 * PI * sqrt(r0*r0*r0 / mu);
    int n_orbits = 5;
    double t_end = n_orbits * T;

    double t0_jd = JD_JAN_2000;         // anchor for ephemeris time

    printf("\n=== Integrator physics breakdown benchmark ===\n");
    printf("Low lunar orbit: r=%.1f km (alt %.0f km), v=%.4f km/s, T=%.1f s\n",
           r0, altitude, v0, T);
    printf("Propagation: %d orbits = %.1f s\n", n_orbits, t_end);
    printf("Integrator: RKDP45 rel_tol=1e-9\n\n");

    OrbitParams base = {
        .mu_moon       = mu,
        .hg            = NULL,
        .eph           = NULL,
        .use_earth_3rd = 0,
        .use_sun_3rd   = 0,
        .srp_AoM_Cr    = 0.0,
        .t0_jd         = t0_jd,
        .n_rhs_calls   = 0
    };

    // 1. Two-body only (baseline)
    run_scenario("2-body only", &base, y0, t_end);

    if (eph_ok) {
        // 2. + Earth 3rd body
        OrbitParams p = base;
        p.eph = &eph_handle;
        p.use_earth_3rd = 1;
        run_scenario("2-body + Earth 3rd body", &p, y0, t_end);

        // 3. + Sun 3rd body
        p = base;
        p.eph = &eph_handle;
        p.use_sun_3rd = 1;
        run_scenario("2-body + Sun 3rd body", &p, y0, t_end);

        // 4. + SRP
        p = base;
        p.eph = &eph_handle;
        p.srp_AoM_Cr = CR_ALUMINUM_COATED_MYLAR * 0.02;   // Cr * (A/m) m^2/kg
        run_scenario("2-body + SRP", &p, y0, t_end);
    } else {
        printf("  (3rd body / SRP scenarios skipped: ephemeris not available)\n");
    }

    // ---- harmonics scenarios (need both ephemeris and gravity file) ----
    if (eph_ok && grav_path) {
        int N_list[] = {10, 50, 100, 200};
        for (size_t i = 0; i < sizeof(N_list)/sizeof(N_list[0]); i++) {
            int N = N_list[i];

            HarmonicGravityData hgd = {0};
            HarmonicGravity     hg  = {0};
            if (spody_load_HarmonicGravityData(&hgd, grav_path, N) != 0 ||
                spody_setup_HarmonicGravity(&hg, &hgd) != 0) {
                printf("  WARN: harmonics N=%d setup failed - skipping\n", N);
                continue;
            }

            OrbitParams p = base;
            p.eph = &eph_handle;
            p.hg  = &hg;
            char label[64];
            snprintf(label, sizeof(label), "Harmonics N=%d", N);
            run_scenario(label, &p, y0, t_end);

            spody_free_HarmonicGravity(&hg);
            spody_free_HarmonicGravityData(&hgd);
        }

        // Full: harmonics N=100 + 3rd bodies + SRP
        {
            int N_full = 100;
            HarmonicGravityData hgd = {0};
            HarmonicGravity     hg  = {0};
            if (spody_load_HarmonicGravityData(&hgd, grav_path, N_full) == 0 &&
                spody_setup_HarmonicGravity(&hg, &hgd) == 0) {

                OrbitParams p = base;
                p.eph = &eph_handle;
                p.hg  = &hg;
                p.use_earth_3rd = 1;
                p.use_sun_3rd   = 1;
                p.srp_AoM_Cr    = CR_ALUMINUM_COATED_MYLAR * 0.02;
                char label[64];
                snprintf(label, sizeof(label), "Full: N=%d + Earth + Sun + SRP", N_full);
                run_scenario(label, &p, y0, t_end);

                spody_free_HarmonicGravity(&hg);
                spody_free_HarmonicGravityData(&hgd);
            }
        }
    } else if (!grav_path) {
        printf("  (harmonics scenarios skipped: gravity file not provided as 2nd argument)\n");
    }

    printf("\n");

    if (eph_ok) {
        spody_free_MappedEphemeris(&eph_handle);
        spody_free_MappedEphemerisData(&med);
    }

    return 0;
}
