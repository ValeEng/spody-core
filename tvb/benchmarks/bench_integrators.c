/*
 * Benchmark: compare RK4 (fixed step, varying h) and RKDP45 (adaptive,
 * varying tolerance) on the 2-body Kepler problem.
 *
 * Reports for each scenario: number of steps, wall time, and the final
 * relative energy drift over a fixed propagation horizon (10 orbits).
 *
 * No external data required.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "spody_core.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef struct {
    double mu;
    long   n_rhs_calls;   // counter
} TwoBodyParams;

static int two_body_rhs(double t, const double *y, double *dy, void *user) {
    (void)t;
    TwoBodyParams *p = (TwoBodyParams*)user;
    p->n_rhs_calls++;

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

static double now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

static void run_scenario(const char *label,
                         spody_integrator_method method,
                         IntegratorOptions *opt,
                         double t_end,
                         const double *y0,
                         double mu) {

    TwoBodyParams params = { .mu = mu, .n_rhs_calls = 0 };

    IntegratorAllData integ;
    if (spody_setup_integrator(&integ, method, opt, 6, two_body_rhs, &params) != SPODY_INTEG_OK) {
        printf("  setup failed\n");
        return;
    }
    spody_set_integrator_state(&integ, 0.0, y0);

    double E0 = specific_energy(y0, mu);

    long n_steps = 0;
    double t0 = now_ms();
    while (integ.t < t_end) {
        if (spody_propagate_onestep(&integ) != SPODY_INTEG_OK) {
            printf("  step failed at t=%.3f\n", integ.t);
            spody_free_integrator(&integ);
            return;
        }
        // clip last step
        if (integ.t + integ.h > t_end) integ.h = t_end - integ.t;
        n_steps++;
    }
    double dt = now_ms() - t0;

    double E1 = specific_energy(integ.y, mu);
    double dE = fabs((E1 - E0) / E0);

    printf("  %-28s | steps=%7ld | rhs_calls=%8ld | time=%8.2f ms | dE/E=%.3e\n",
           label, n_steps, params.n_rhs_calls, dt, dE);

    spody_free_integrator(&integ);
}

int main(void) {

    double r0 = 7000.0;
    double mu = EARTH_MU;
    double v0 = sqrt(mu / r0);
    double y0[6] = { r0, 0.0, 0.0,  0.0, v0, 0.0 };

    double T = 2.0 * PI * sqrt(r0*r0*r0 / mu);
    int n_orbits = 10;
    double t_end = n_orbits * T;

    printf("\n=== Integrator benchmark: %d Kepler orbits, T=%.3f s each ===\n", n_orbits, T);
    printf("Total horizon: %.1f s\n\n", t_end);

    // ---- RK4 sweep on fixed step ----
    printf("RK4 (fixed step):\n");
    double h_list[] = {30.0, 10.0, 5.0, 1.0};
    for (int i = 0; i < (int)(sizeof(h_list)/sizeof(h_list[0])); i++) {
        IntegratorOptions opt;
        spody_default_integrator_options(SPODY_INTEG_RK4, &opt);
        opt.h_init = h_list[i];
        char label[64];
        snprintf(label, sizeof(label), "h=%.1f s", h_list[i]);
        run_scenario(label, SPODY_INTEG_RK4, &opt, t_end, y0, mu);
    }
    printf("\n");

    // ---- RKDP45 sweep on tolerance ----
    printf("RKDP45 (adaptive step):\n");
    double tol_list[] = {1e-6, 1e-9, 1e-12};
    for (int i = 0; i < (int)(sizeof(tol_list)/sizeof(tol_list[0])); i++) {
        IntegratorOptions opt;
        spody_default_integrator_options(SPODY_INTEG_RK45, &opt);
        opt.rel_tol = tol_list[i];
        opt.h_init  = 30.0;
        opt.h_max   = 600.0;
        char label[64];
        snprintf(label, sizeof(label), "rel_tol=%.0e", tol_list[i]);
        run_scenario(label, SPODY_INTEG_RK45, &opt, t_end, y0, mu);
    }
    printf("\n");

    return 0;
}
