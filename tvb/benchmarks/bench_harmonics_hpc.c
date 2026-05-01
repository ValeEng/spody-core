/*
 * Benchmark: spody_get_hgaccbodyfixed (reference) vs.
 *            spody_get_hgaccbodyfixed_hpc (production high-performance variant).
 *
 * For each degree N, both functions are called many times on a fixed
 * position; total wall time and ns per call are reported.
 *
 * The hpc variant carries SIMD pragmas that take effect only when the
 * project is compiled with SPODY_ENABLE_OMP_SIMD=ON. Without the flag
 * hpc still benefits from branch-free + peeled layout, but no explicit
 * vectorization.
 *
 * Usage:
 *   bench_harmonics_hpc <path_to_grgm1200a_sha.tab>
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "spody_core.h"

static double now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

static void run_one(const char *grav_path, int N, long iters,
                    const double pos[3]) {

    HarmonicGravityData hgd = {0};
    HarmonicGravity     hg  = {0};
    if (spody_load_HarmonicGravityData(&hgd, grav_path, N) != 0 ||
        spody_setup_HarmonicGravity(&hg, &hgd) != 0) {
        printf("  N=%4d: setup failed\n", N);
        spody_free_HarmonicGravity(&hg);
        spody_free_HarmonicGravityData(&hgd);
        return;
    }

    double sink = 0.0;       // prevent dead-code elimination
    double acc[3];

    // ---- reference ----
    double t0 = now_ms();
    for (long i = 0; i < iters; i++) {
        spody_get_hgaccbodyfixed(&hg, (double*)pos, acc);
        sink += acc[0] + acc[1] + acc[2];
    }
    double t_ref = now_ms() - t0;

    // ---- high-performance variant ----
    t0 = now_ms();
    for (long i = 0; i < iters; i++) {
        spody_get_hgaccbodyfixed_hpc(&hg, (double*)pos, acc);
        sink += acc[0] + acc[1] + acc[2];
    }
    double t_hpc = now_ms() - t0;

    double ns_ref = (iters > 0) ? (t_ref * 1e6 / (double)iters) : 0.0;
    double ns_hpc = (iters > 0) ? (t_hpc * 1e6 / (double)iters) : 0.0;
    double speedup = (t_hpc > 0.0) ? t_ref / t_hpc : 0.0;

    printf("  N=%4d (%6ld iters) | ref %8.1f ms (%7.1f ns/call) | hpc %8.1f ms (%7.1f ns/call) | speedup %5.2fx | sink=%.6e\n",
           N, iters, t_ref, ns_ref, t_hpc, ns_hpc, speedup, sink);

    spody_free_HarmonicGravity(&hg);
    spody_free_HarmonicGravityData(&hgd);
}

int main(int argc, char **argv) {

    if (argc < 2) {
        printf("usage: %s <path_to_grgm1200a_sha.tab>\n", argv[0]);
        return 1;
    }

    double pos[3] = { 1500.0, 900.0, 600.0 }; // generic position in PA frame

    printf("\n=== Harmonic gravity hpc benchmark ===\n");
    printf("Position: (%.1f, %.1f, %.1f) km, |r| = %.3f km\n\n",
           pos[0], pos[1], pos[2],
           sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]));

    struct { int N; long iters; } cases[] = {
        {  10, 200000 },
        {  50,  50000 },
        { 100,  10000 },
        { 200,   3000 },
    };
    int n_cases = (int)(sizeof(cases)/sizeof(cases[0]));

    for (int i = 0; i < n_cases; i++) {
        run_one(argv[1], cases[i].N, cases[i].iters, pos);
    }
    printf("\n");
    return 0;
}
