/*
 * Benchmark: compare single-call vs batch vs single-call-no-cache
 * for spody_get_ephposition on multiple targets sharing the same JD.
 *
 * Usage:
 *   bench_ephemeris_batch <path_to_de440.spody>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "spody_ephemeris.h"

#define N_DATES 100000

static double now_ms(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC * 1000.0;
}

static void invalidate_cache(MappedEphemeris *map) {
    for (int i = 0; i < EPH_CACHE_SLOTS; i++) map->cache_valid[i] = 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("usage: %s <path_to_de440.spody>\n", argv[0]);
        return 1;
    }

    MappedEphemeris map = {0};
    if (spody_setup_MappedEphemeris(&map, argv[1]) != 0) {
        printf("failed to setup ephemeris\n");
        return 1;
    }

    // pick a safe JD window well inside the file
    double jd_start = map.header->start_epoch + 100.0;
    double jd_step  = 1.0;

    // from Moon: Sun, Earth, Jupiter, Saturn
    int targets[] = {10, 399, 5, 6};
    const int n_targets = sizeof(targets) / sizeof(targets[0]);

    double sink = 0.0; // prevent dead-code elimination
    double t0, tA, tB, tC;

    /* -------- Scenario A: single calls, cache forcibly invalidated
                between each call (simulates pre-cache behavior) -------- */
    invalidate_cache(&map);
    t0 = now_ms();
    for (int k = 0; k < N_DATES; k++) {
        double jd = jd_start + k * jd_step;
        double r[3];
        for (int t = 0; t < n_targets; t++) {
            invalidate_cache(&map);
            spody_get_ephposition(&map, 301, targets[t], jd, r);
            sink += r[0] + r[1] + r[2];
        }
    }
    tA = now_ms() - t0;

    /* -------- Scenario B: single calls, cache on (current default) -------- */
    invalidate_cache(&map);
    t0 = now_ms();
    for (int k = 0; k < N_DATES; k++) {
        double jd = jd_start + k * jd_step;
        double r[3];
        for (int t = 0; t < n_targets; t++) {
            spody_get_ephposition(&map, 301, targets[t], jd, r);
            sink += r[0] + r[1] + r[2];
        }
    }
    tB = now_ms() - t0;

    /* -------- Scenario C: batch API, cache on -------- */
    invalidate_cache(&map);
    double rb[12]; // flat buffer: 4 targets * 3 components
    t0 = now_ms();
    for (int k = 0; k < N_DATES; k++) {
        double jd = jd_start + k * jd_step;
        spody_get_ephposition_batch(&map, 301, targets, n_targets, jd, rb);
        for (int i = 0; i < 3 * n_targets; i++) {
            sink += rb[i];
        }
    }
    tC = now_ms() - t0;

    printf("\n=== Ephemeris query benchmark ===\n");
    printf("Iterations      : %d dates * %d bodies = %d queries\n",
           N_DATES, n_targets, N_DATES * n_targets);
    printf("\n");
    printf("A) single calls, no cache   : %8.2f ms  (%6.0f ns / query)\n",
           tA, tA * 1e6 / (N_DATES * n_targets));
    printf("B) single calls, with cache : %8.2f ms  (%6.0f ns / query)\n",
           tB, tB * 1e6 / (N_DATES * n_targets));
    printf("C) batch API, with cache    : %8.2f ms  (%6.0f ns / query)\n",
           tC, tC * 1e6 / (N_DATES * n_targets));
    printf("\n");
    printf("speedup B vs A (cache)      : %.2fx\n", tA / tB);
    printf("speedup C vs A (cache+batch): %.2fx\n", tA / tC);
    printf("speedup C vs B (batch only) : %.2fx\n", tB / tC);
    printf("\nsink (ignore) = %.6e\n", sink);

    return 0;
}
