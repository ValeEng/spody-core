/*
 * Regression test for the Moon->Sun computation and related invariants.
 *
 * Usage:
 *   test_ephemeris_moon_sun <path_to_de440.spody>
 *
 * Exit code 0 on success, non-zero on any failure.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "spody_ephemeris.h"

static double vec_norm(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("usage: %s <path_to_de440.spody>\n", argv[0]);
        return 1;
    }

    MappedEphemeris map = {0};
    if (spody_setup_MappedEphemeris(&map, argv[1]) != 0) {
        printf("FAIL: setup\n");
        return 1;
    }

    // some JD well inside the file
    double jd = map.header->start_epoch + 100.0;

    double earth_sun[3], moon_sun[3], earth_moon[3], moon_earth[3];
    spody_get_ephposition(&map, 399, 10,  jd, earth_sun);
    spody_get_ephposition(&map, 301, 10,  jd, moon_sun);
    spody_get_ephposition(&map, 399, 301, jd, earth_moon);
    spody_get_ephposition(&map, 301, 399, jd, moon_earth);

    double d_earth_sun  = vec_norm(earth_sun);
    double d_moon_sun   = vec_norm(moon_sun);
    double d_earth_moon = vec_norm(earth_moon);

    printf("|earth->sun|   = %12.3f km (expected ~1.496e8, perihelion/aphelion 1.47-1.52e8)\n", d_earth_sun);
    printf("|moon->sun|    = %12.3f km (expected ~1.496e8 +/- ~4e5)\n",               d_moon_sun);
    printf("|earth->moon|  = %12.3f km (expected ~384400, perigee/apogee 356k-407k)\n",d_earth_moon);

    int fail = 0;

    // 1. earth-sun in reasonable AU range
    if (d_earth_sun < 1.45e8 || d_earth_sun > 1.55e8) {
        printf("FAIL: earth-sun distance out of expected range\n");
        fail = 1;
    }

    // 2. earth-moon in reasonable range
    if (d_earth_moon < 350000.0 || d_earth_moon > 410000.0) {
        printf("FAIL: earth-moon distance out of expected range\n");
        fail = 1;
    }

    // 3. moon-sun and earth-sun must differ. They can differ at most by |earth->moon|.
    double diff = fabs(d_moon_sun - d_earth_sun);
    if (diff < 1.0) {
        printf("FAIL: moon-sun == earth-sun. The old bug is back!\n");
        fail = 1;
    }
    if (diff > 500000.0) {
        printf("FAIL: |moon-sun| and |earth-sun| differ by %.0f km, max expected ~4e5\n", diff);
        fail = 1;
    }

    // 4. symmetry: r(A->B) = -r(B->A)
    for (int i = 0; i < 3; i++) {
        if (fabs(earth_moon[i] + moon_earth[i]) > 1e-6) {
            printf("FAIL: earth->moon != -(moon->earth) on component %d\n", i);
            fail = 1;
        }
    }

    // 5. triangular closure: r(Moon->Sun) = r(Earth->Sun) + r(Moon->Earth)
    for (int i = 0; i < 3; i++) {
        double expected = earth_sun[i] + moon_earth[i];
        double err = fabs(expected - moon_sun[i]);
        if (err > 1e-3) {
            printf("FAIL: triangular closure off by %.6e km on component %d\n", err, i);
            fail = 1;
        }
    }

    // 6. batch API must agree with single-call API (bit-identical)
    int targets[] = {10, 399, 5, 6};
    double batch_result[12]; // flat: 4 targets * 3 components
    spody_get_ephposition_batch(&map, 301, targets, 4, jd, batch_result);
    double ref[4][3];
    spody_get_ephposition(&map, 301, 10,  jd, ref[0]);
    spody_get_ephposition(&map, 301, 399, jd, ref[1]);
    spody_get_ephposition(&map, 301, 5,   jd, ref[2]);
    spody_get_ephposition(&map, 301, 6,   jd, ref[3]);
    for (int t = 0; t < 4; t++) {
        for (int i = 0; i < 3; i++) {
            double err = fabs(batch_result[3*t + i] - ref[t][i]);
            if (err > 1e-6) {
                printf("FAIL: batch vs single differ on target[%d] comp %d by %.6e km\n",
                       t, i, err);
                fail = 1;
            }
        }
    }

    if (!fail) printf("PASS\n");
    return fail;
}
