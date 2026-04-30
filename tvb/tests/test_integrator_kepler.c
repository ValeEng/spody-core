/*
 * Regression test for the integrators on the 2-body Kepler problem.
 *
 * Setup:
 *   - circular LEO around the Earth (r0 = 7000 km, v0 = sqrt(mu/r0))
 *   - propagate for one orbital period
 *   - check that the satellite returns near the initial state
 *   - check that the specific orbital energy is conserved
 *
 * No external data required.
 *
 * Exit code 0 on success, non-zero on any failure.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spody_core.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef struct {
    double mu;
} TwoBodyParams;

static int two_body_rhs(double t, const double *y, double *dy, void *user) {
    (void)t;
    TwoBodyParams *p = (TwoBodyParams*)user;

    double r2 = y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
    double r3 = r2 * sqrt(r2);
    double k  = -p->mu / r3;

    dy[0] = y[3];
    dy[1] = y[4];
    dy[2] = y[5];
    dy[3] = k * y[0];
    dy[4] = k * y[1];
    dy[5] = k * y[2];

    return 0;
}

static double specific_energy(const double *y, double mu) {
    double r  = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    double v2 = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    return 0.5 * v2 - mu / r;
}

static int run_one_method(spody_integrator_method method,
                          const char *name,
                          const IntegratorOptions *opt,
                          double tol_pos_km,
                          double tol_energy_rel) {
    TwoBodyParams params = { .mu = EARTH_MU };

    // initial state: circular LEO in the equatorial plane
    double r0 = 7000.0;                        // km
    double v0 = sqrt(params.mu / r0);          // km/s, circular
    double y0[6] = { r0, 0.0, 0.0,  0.0, v0, 0.0 };
    double E0 = specific_energy(y0, params.mu);

    // one orbital period
    double T = 2.0 * PI * sqrt(r0*r0*r0 / params.mu);

    IntegratorAllData integ;
    int rc = spody_setup_integrator(&integ, method, opt, 6, two_body_rhs, &params);
    if (rc != SPODY_INTEG_OK) {
        printf("FAIL [%s]: setup returned %d\n", name, rc);
        return 1;
    }

    spody_set_integrator_state(&integ, 0.0, y0);
    rc = spody_propagate_untilend(&integ, T);
    if (rc != SPODY_INTEG_OK) {
        printf("FAIL [%s]: propagate returned %d (t=%.3f, expected %.3f)\n",
               name, rc, integ.t, T);
        spody_free_integrator(&integ);
        return 1;
    }

    // check 1: position returns close to initial
    double dx = integ.y[0] - y0[0];
    double dy = integ.y[1] - y0[1];
    double dz = integ.y[2] - y0[2];
    double pos_err = sqrt(dx*dx + dy*dy + dz*dz);

    // check 2: specific orbital energy is conserved
    double E1 = specific_energy(integ.y, params.mu);
    double energy_rel_err = fabs((E1 - E0) / E0);

    int fail = 0;
    printf("[%s] T=%.3f s | pos_err=%.6e km | dE/E=%.6e\n",
           name, T, pos_err, energy_rel_err);

    if (pos_err > tol_pos_km) {
        printf("  FAIL: pos_err %.6e km > tol %.6e km\n", pos_err, tol_pos_km);
        fail = 1;
    }
    if (energy_rel_err > tol_energy_rel) {
        printf("  FAIL: dE/E %.6e > tol %.6e\n", energy_rel_err, tol_energy_rel);
        fail = 1;
    }

    spody_free_integrator(&integ);
    return fail;
}

int main(void) {

    int fail = 0;

    // RK4 fixed step at h = 5 s -> ~ 1180 steps per orbit
    {
        IntegratorOptions opt;
        spody_default_integrator_options(SPODY_INTEG_RK4, &opt);
        opt.h_init = 5.0;
        fail |= run_one_method(SPODY_INTEG_RK4, "RK4 (h=5s)", &opt,
                               /*pos_err <*/ 5.0,
                               /*dE/E   <*/ 1e-7);
    }

    // RK4 with tighter step
    {
        IntegratorOptions opt;
        spody_default_integrator_options(SPODY_INTEG_RK4, &opt);
        opt.h_init = 1.0;
        fail |= run_one_method(SPODY_INTEG_RK4, "RK4 (h=1s)", &opt,
                               /*pos_err <*/ 1.0e-2,
                               /*dE/E   <*/ 1e-10);
    }

    // RKDP45 adaptive at default tolerance
    {
        IntegratorOptions opt;
        spody_default_integrator_options(SPODY_INTEG_RK45, &opt);
        opt.rel_tol = 1e-10;
        opt.h_init  = 30.0;
        opt.h_max   = 200.0;
        fail |= run_one_method(SPODY_INTEG_RK45, "RKDP45 (rel_tol=1e-10)", &opt,
                               /*pos_err <*/ 1.0e-2,
                               /*dE/E   <*/ 1e-10);
    }

    if (!fail) printf("\nPASS\n");
    return fail;
}
