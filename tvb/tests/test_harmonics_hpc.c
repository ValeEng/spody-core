/*
 * Numerical regression test: high-performance harmonic gravity vs. reference.
 *
 * Verifies that spody_get_hgaccbodyfixed_hpc produces results bit-identical
 * (or within strict epsilon-machine tolerance) to spody_get_hgaccbodyfixed
 * over a sampled set of positions and degrees.
 *
 * Usage:
 *   test_harmonics_hpc <path_to_grgm1200a_sha.tab>
 *
 * Exit code 0 on PASS, non-zero on any failure.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "spody_core.h"

static int compare_against_ref(HarmonicGravity *hg, double pos[3],
                               int n_idx, int p_idx,
                               double abs_tol, double rel_tol) {
    double acc_ref[3], acc_hpc[3];
    spody_get_hgaccbodyfixed    (hg, pos, acc_ref);
    spody_get_hgaccbodyfixed_hpc(hg, pos, acc_hpc);

    int fail = 0;
    for (int i = 0; i < 3; i++) {
        double diff = fabs(acc_hpc[i] - acc_ref[i]);
        double scale = fmax(fabs(acc_ref[i]), fabs(acc_hpc[i]));
        double rel  = (scale > 0.0) ? diff / scale : 0.0;

        if (diff > abs_tol && rel > rel_tol) {
            printf("  FAIL n_idx=%d p_idx=%d comp=%d: ref=%.21e hpc=%.21e diff=%.6e rel=%.6e\n",
                   n_idx, p_idx, i, acc_ref[i], acc_hpc[i], diff, rel);
            fail = 1;
        }
    }
    return fail;
}

int main(int argc, char **argv) {

    if (argc < 2) {
        printf("usage: %s <path_to_grgm1200a_sha.tab>\n", argv[0]);
        return 1;
    }

    int N_list[] = {2, 10, 50, 100, 200};
    int n_N = (int)(sizeof(N_list)/sizeof(N_list[0]));

    // sample positions covering different geometries (km, lunar PA frame).
    // R_moon ~ 1738 km. Test from low orbit out to high orbit.
    double positions[][3] = {
        { 1837.0,    0.0,    0.0 },     // 100 km LLO equator, x-axis
        {    0.0, 1837.0,    0.0 },     // y-axis
        {    0.0,    0.0, 1837.0 },     // polar overhead (z-axis: triggers near-pole behavior)
        { 1500.0,  900.0,  600.0 },     // generic position
        { 5000.0,    0.0,    0.0 },     // medium altitude
        { -2500.0, 3000.0, 1000.0 },    // negative x, mixed
        { 100000.0, 100.0, 100.0 },     // far field (Earth-Moon distance)
    };
    int n_pos = (int)(sizeof(positions)/sizeof(positions[0]));

    int total_fail = 0;
    double abs_tol = 1e-15;
    double rel_tol = 1e-12;

    for (int j = 0; j < n_N; j++) {
        int N = N_list[j];

        HarmonicGravityData hgd = {0};
        HarmonicGravity     hg  = {0};
        if (spody_load_HarmonicGravityData(&hgd, argv[1], N) != 0) {
            printf("FAIL: cannot load gravity data N=%d\n", N);
            spody_free_HarmonicGravityData(&hgd);
            return 1;
        }
        if (spody_setup_HarmonicGravity(&hg, &hgd) != 0) {
            printf("FAIL: cannot setup gravity N=%d\n", N);
            spody_free_HarmonicGravityData(&hgd);
            return 1;
        }

        int N_fail = 0;
        for (int p = 0; p < n_pos; p++) {
            N_fail += compare_against_ref(&hg, positions[p], j, p, abs_tol, rel_tol);
        }

        printf("[N=%4d] hpc: %s (over %d positions, rel < %.0e)\n",
               N, (N_fail == 0) ? "PASS" : "FAIL", n_pos, rel_tol);

        total_fail += N_fail;

        spody_free_HarmonicGravity(&hg);
        spody_free_HarmonicGravityData(&hgd);
    }

    if (total_fail == 0) {
        printf("\nPASS\n");
        return 0;
    } else {
        printf("\nFAIL: %d total mismatches\n", total_fail);
        return 1;
    }
}
