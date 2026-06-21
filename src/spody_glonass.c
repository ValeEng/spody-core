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
/*
 * Implementation of spody_convert_glonass_to_state_icrf. See
 * spody_glonass.h for the format references (RINEX 3.05 sect. 6 for
 * GLONASS nav messages; PZ-90 frame; IERS nominal omega_earth) and
 * the ECEF -> ICRF state-vector transformation.
 *
 * RINEX 3.x GLONASS record layout (5 lines as written by BKG / IGS
 * since RINEX 3.05). Each numeric field is in Fortran D-format which
 * sscanf %lf reads natively after a D->E swap:
 *
 *   Line 1 (PRN row):     R<nn> YYYY MM DD HH MM SS  clk_bias clk_drift TOM
 *   Line 2 (X):           <indent> pos_x[km] vel_x[km/s] acc_x[km/s^2] health
 *   Line 3 (Y):           <indent> pos_y[km] vel_y[km/s] acc_y[km/s^2] frqnum
 *   Line 4 (Z):           <indent> pos_z[km] vel_z[km/s] acc_z[km/s^2] age_op
 *   Line 5 (status):      <indent> status_flags delta_tau_n urai health_flags
 *
 * The 5th line was added in RINEX 3.05; older 3.x files (3.0 - 3.04)
 * omit it. We default to reading 5 lines per record (current IGS
 * practice). A future revision can sniff "RINEX VERSION / TYPE" in
 * the header and adapt.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_glonass.h"
#include "spody_eop.h"
#include "spody_earth_orientation.h"
#include "spody_forcemodels.h"
#include "spody_const.h"

#define SPODY_GLONASS_OUT_MAGIC      "SPDYOUT_"
#define SPODY_GLONASS_OUT_VERSION    1u
#define SPODY_GLONASS_OUT_STATE_DIM  6u

/* RINEX-NAV nav lines fit in 80 cols + newline + sentinel; 256 covers
 * any RINEX 3.x version with comfortable slack. */
#define SPODY_RINEX_LINE             256

/* TAI - UTC, constant since the 2017-01-01 leap second. For dates
 * before 2017 a proper leap-second table (spody_eop.c has one) would
 * be needed; this converter is documented to target post-2017 data. */
#define SPODY_TAI_MINUS_UTC_SEC      37.0

/* TT - TAI is exact by IAU definition. TDB - TT is <2 ms, far below
 * broadcast precision (~few m / cm-per-s), so we use TDB ~ TT. */
#define SPODY_TT_MINUS_TAI_SEC       32.184

/* Nominal sidereal Earth rotation rate, IERS Conventions 2010 sec.
 * 1.4. Polar-motion corrections to the rotation axis are <1e-7 of
 * this magnitude and ignored here (broadcast (r, v) is itself only
 * decimetre / cm-per-s, well above any omega correction). */
#define SPODY_OMEGA_EARTH_RADPS      7.2921151467e-5

#define SPODY_JD_J2000_TT            2451545.0

/* Fortran-D exponent normaliser for sscanf. RINEX 3 still emits
 * 1.23D-04 -- swap to E in place (private buffer copy). */
static double _parse_rinex_double(const char *token) {
    char buf[40];
    size_t i = 0;
    for (; token[i] && i < sizeof(buf) - 1; ++i) {
        char c = token[i];
        if (c == 'D' || c == 'd') c = 'E';
        buf[i] = c;
    }
    buf[i] = '\0';
    return strtod(buf, NULL);
}

/* Meeus Gregorian -> JD. Returns JD at midnight + day fraction. The
 * y/m/d/hh/mn/ss values are taken verbatim from the RINEX TOC line:
 * for GLONASS records that is UTC per RINEX 3 spec sect. 6.10.5. */
static double _greg_to_jd(int y, int m, int d, int hh, int mn, double ss) {
    if (m <= 2) { y -= 1; m += 12; }
    int A = y / 100;
    int B = 2 - A + (A / 4);
    double jd_midnight = floor(365.25 * (double)(y + 4716))
                        + floor(30.6001 * (double)(m + 1))
                        + (double)d + (double)B - 1524.5;
    double day_frac = ((double)hh * 3600.0 + (double)mn * 60.0 + ss) / 86400.0;
    return jd_midnight + day_frac;
}

/* 3x3 matrix-vector multiply (out = R . v). */
static void _rot3_apply(const double R[3][3], const double v[3], double out[3]) {
    out[0] = R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2];
    out[1] = R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2];
    out[2] = R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2];
}

/* SPDYOUT_ 24-byte preamble identical to sim_run.c's writer. */
static int _write_glonass_out_header(FILE *fp) {
    if (fwrite(SPODY_GLONASS_OUT_MAGIC, 1, 8, fp) != 8) return -1;
    uint32_t hdr[4] = {
        SPODY_GLONASS_OUT_VERSION,
        SPODY_GLONASS_OUT_STATE_DIM,
        0u, 0u
    };
    if (fwrite(hdr, sizeof(uint32_t), 4, fp) != 4) return -1;
    return 0;
}

int spody_convert_glonass_to_state_icrf(const char *input_rnx,
                                        const char *output_bin,
                                        const char *sat_id,
                                        const char *eop_file,
                                        const char *iau2006_dir) {
    if (!input_rnx || !output_bin || !sat_id || !eop_file || !iau2006_dir) {
        fprintf(stderr, "glonass: NULL argument\n");
        return 1;
    }
    if (strlen(sat_id) != 3 || (sat_id[0] != 'R')) {
        fprintf(stderr,
            "glonass: sat_id must be 3 chars starting with 'R' "
            "(got '%s')\n", sat_id);
        return 1;
    }

    /* --- Bring up EOP + IAU 2006 -------------------------------- */
    MappedEOPData eop_data = {0};
    if (spody_setup_MappedEOPData(&eop_data, eop_file) != 0) {
        fprintf(stderr, "glonass: cannot load EOP from '%s'\n", eop_file);
        return 1;
    }
    MappedIAU2006Data iau_data = {0};
    if (spody_setup_MappedIAU2006Data(&iau_data, iau2006_dir) != 0) {
        fprintf(stderr, "glonass: cannot load IAU 2006 from '%s'\n",
                iau2006_dir);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    MappedEOP     eop_map = {0};
    MappedIAU2006 iau_map = {0};
    spody_setup_MappedEOP(&eop_map, &eop_data);
    spody_setup_MappedIAU2006(&iau_map, &iau_data);
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    /* --- IO ----------------------------------------------------- */
    FILE *fin = fopen(input_rnx, "r");
    if (!fin) {
        fprintf(stderr, "glonass: cannot open input '%s'\n", input_rnx);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        fprintf(stderr, "glonass: cannot open output '%s'\n", output_bin);
        fclose(fin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    if (_write_glonass_out_header(fout) != 0) {
        fprintf(stderr, "glonass: cannot write output header\n");
        fclose(fout); fclose(fin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    /* --- Skip header until "END OF HEADER" ---------------------- */
    char line[SPODY_RINEX_LINE];
    int in_body = 0;
    while (fgets(line, sizeof line, fin)) {
        /* The "END OF HEADER" marker is right-justified in the
         * RINEX label column (col 61..). We do a substring search
         * because the leading bytes are version-dependent padding. */
        if (strstr(line, "END OF HEADER")) {
            in_body = 1;
            break;
        }
    }
    if (!in_body) {
        fprintf(stderr,
            "glonass: input '%s' has no 'END OF HEADER' marker\n",
            input_rnx);
        fclose(fout); fclose(fin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    /* --- Body: scan for sat_id records, write SPDYOUT_ rows ----- */
    size_t n_total      = 0;
    size_t n_written    = 0;
    double et_first     = 0.0, et_last = 0.0;

    while (fgets(line, sizeof line, fin)) {
        /* PRN row begins with the 3-char sat id at column 0. Other
         * rows start with whitespace. */
        if (line[0] != 'R') continue;

        char rec_id[4] = {0};
        memcpy(rec_id, line, 3);
        ++n_total;
        int sat_match = (strcmp(rec_id, sat_id) == 0);

        /* Parse TOC from line 1 regardless: epoch matters for the
         * frame rotation even on non-target records (cheap, keeps
         * code linear without conditional skips). For non-target
         * records we still need to consume the 4 continuation lines
         * so the body sweep doesn't run away. */
        int    y, mo, d, h, mi;
        double sec;
        double clk_bias, clk_drift, tom;
        /* Line 1 layout after the PRN: " YYYY MM DD HH MM SS  D D D" */
        if (sscanf(line + 3, "%d %d %d %d %d %lf %lf %lf %lf",
                   &y, &mo, &d, &h, &mi, &sec,
                   &clk_bias, &clk_drift, &tom) < 6) {
            continue;
        }
        (void)clk_bias; (void)clk_drift; (void)tom;

        /* Read the 4 continuation lines (X / Y / Z / status). The
         * 5th line is the RINEX 3.05 "status_flags / delta_tau" row
         * and is silently consumed too. */
        char xl[SPODY_RINEX_LINE], yl[SPODY_RINEX_LINE];
        char zl[SPODY_RINEX_LINE], sl[SPODY_RINEX_LINE];
        if (!fgets(xl, sizeof xl, fin)) break;
        if (!fgets(yl, sizeof yl, fin)) break;
        if (!fgets(zl, sizeof zl, fin)) break;
        if (!fgets(sl, sizeof sl, fin)) break;
        (void)sl;

        if (!sat_match) continue;

        /* Each continuation line is "  posα  velα  accα  flagα" --
         * 4 numeric fields after the indent. */
        double px, vx, ax, hx;
        double py, vy, ay, hy;
        double pz, vz, az, hz;
        if (sscanf(xl, "%lf %lf %lf %lf", &px, &vx, &ax, &hx) != 4 ||
            sscanf(yl, "%lf %lf %lf %lf", &py, &vy, &ay, &hy) != 4 ||
            sscanf(zl, "%lf %lf %lf %lf", &pz, &vz, &az, &hz) != 4) {
            continue;
        }
        (void)ax; (void)ay; (void)az;
        (void)hx; (void)hy; (void)hz;

        /* Some legacy GLONASS files use Fortran D exponents in the
         * data rows. sscanf %lf on most C runtimes rejects D, so we
         * re-parse the field via _parse_rinex_double when sscanf
         * looks suspicious. In practice modern IGS files use E and
         * sscanf works; we keep the explicit reparse as a safety net. */
        if (px == 0.0 && (strchr(xl, 'D') || strchr(xl, 'd'))) {
            char buf[SPODY_RINEX_LINE]; size_t i;
            for (i = 0; i + 1 < sizeof buf && xl[i]; ++i)
                buf[i] = (xl[i] == 'D' || xl[i] == 'd') ? 'E' : xl[i];
            buf[i] = '\0';
            sscanf(buf, "%lf %lf %lf %lf", &px, &vx, &ax, &hx);
        }
        if (py == 0.0 && (strchr(yl, 'D') || strchr(yl, 'd'))) {
            char buf[SPODY_RINEX_LINE]; size_t i;
            for (i = 0; i + 1 < sizeof buf && yl[i]; ++i)
                buf[i] = (yl[i] == 'D' || yl[i] == 'd') ? 'E' : yl[i];
            buf[i] = '\0';
            sscanf(buf, "%lf %lf %lf %lf", &py, &vy, &ay, &hy);
        }
        if (pz == 0.0 && (strchr(zl, 'D') || strchr(zl, 'd'))) {
            char buf[SPODY_RINEX_LINE]; size_t i;
            for (i = 0; i + 1 < sizeof buf && zl[i]; ++i)
                buf[i] = (zl[i] == 'D' || zl[i] == 'd') ? 'E' : zl[i];
            buf[i] = '\0';
            sscanf(buf, "%lf %lf %lf %lf", &pz, &vz, &az, &hz);
        }
        (void)_parse_rinex_double;  /* kept for documentation symmetry */

        /* --- Time UTC -> ET TDB --------------------------------- */
        double jd_utc = _greg_to_jd(y, mo, d, h, mi, sec);
        double jd_tt  = jd_utc +
                        (SPODY_TAI_MINUS_UTC_SEC + SPODY_TT_MINUS_TAI_SEC)
                        / 86400.0;
        double et     = (jd_tt - SPODY_JD_J2000_TT) * 86400.0;

        /* --- Rotation ECEF -> ICRF ----------------------------- */
        double R_i2bf[3][3], R_bf2i[3][3];
        spody_bf_rotation_earth(&ctx, et, R_i2bf, R_bf2i);

        double r_ecef[3] = { px, py, pz };
        double v_ecef[3] = { vx, vy, vz };
        double r_icrf[3], v_icrf_rot[3];
        _rot3_apply(R_bf2i, r_ecef, r_icrf);
        _rot3_apply(R_bf2i, v_ecef, v_icrf_rot);

        /* omega_earth x r_ICRF with the FULL rotation axis. At J2024 the
         * ITRS z-axis (= Earth's mean rotation axis) is tilted by
         * X ~ 480 arcsec from the ICRF z-axis due to precession; using
         * a nominal omega = (0, 0, omega_earth) would miss ~omega*X*r ~
         * 4 m/s at GLONASS altitude. The true axis in ICRF is the third
         * column of R_bf_to_icrf (image of ITRS z-hat under the rotation).
         * Polar motion of the actual instantaneous axis from the mean
         * ITRS axis is sub-arcsec and ignored. */
        double omega_icrf[3] = {
            SPODY_OMEGA_EARTH_RADPS * R_bf2i[0][2],
            SPODY_OMEGA_EARTH_RADPS * R_bf2i[1][2],
            SPODY_OMEGA_EARTH_RADPS * R_bf2i[2][2],
        };
        double v_icrf[3];
        v_icrf[0] = v_icrf_rot[0]
                  + omega_icrf[1] * r_icrf[2] - omega_icrf[2] * r_icrf[1];
        v_icrf[1] = v_icrf_rot[1]
                  + omega_icrf[2] * r_icrf[0] - omega_icrf[0] * r_icrf[2];
        v_icrf[2] = v_icrf_rot[2]
                  + omega_icrf[0] * r_icrf[1] - omega_icrf[1] * r_icrf[0];

        double rec[7] = {
            et,
            r_icrf[0], r_icrf[1], r_icrf[2],
            v_icrf[0], v_icrf[1], v_icrf[2]
        };
        if (fwrite(rec, sizeof(double), 7, fout) != 7) {
            fprintf(stderr,
                "glonass: short write at record %zu (et=%.6f)\n",
                n_written, et);
            fclose(fout); fclose(fin);
            spody_free_MappedIAU2006(&iau_map);
            spody_free_MappedEOP(&eop_map);
            spody_free_MappedIAU2006Data(&iau_data);
            spody_free_MappedEOPData(&eop_data);
            return 1;
        }
        if (n_written == 0) et_first = et;
        et_last = et;
        ++n_written;
    }

    fclose(fout);
    fclose(fin);

    if (n_written == 0) {
        fprintf(stderr,
            "glonass: WARNING -- no records written for sat_id '%s' "
            "(scanned %zu total RINEX records in '%s')\n",
            sat_id, n_total, input_rnx);
    } else {
        double duration_h = (et_last - et_first) / 3600.0;
        fprintf(stderr,
            "glonass: '%s' -> %zu records (sat=%s, et=%.6f..%.6f, "
            "%.3f h, scanned %zu nav messages)\n",
            input_rnx, n_written, sat_id, et_first, et_last,
            duration_h, n_total);
    }

    spody_free_MappedIAU2006(&iau_map);
    spody_free_MappedEOP(&eop_map);
    spody_free_MappedIAU2006Data(&iau_data);
    spody_free_MappedEOPData(&eop_data);
    return 0;
}
