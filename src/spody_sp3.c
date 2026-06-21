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
 * Implementation of spody_convert_sp3_to_state_icrf. See spody_sp3.h
 * for the format references (Hilla 2010 / 2016 SP3-c / SP3-d) and the
 * contract.
 *
 * Wire format we EMIT (must match spody/src/sim_run.c):
 *   8 bytes : ASCII "SPDYOUT_" (no NUL)
 *   4 bytes : uint32 LE version = 1
 *   4 bytes : uint32 LE state_dim = 6
 *   8 bytes : reserved (two uint32 zero)
 *   then per record: 7 little-endian doubles
 *     (t [s past J2000 TDB], x, y, z [km, ICRF], vx, vy, vz [km/s])
 *   We always write vx=vy=vz=0 because the SP3 'P' record carries
 *   position only; downstream diff routines compare positions only
 *   for SP3-derived references.
 *
 * SP3 wire format we READ (subset sufficient for IGS orbit files):
 *   - Header rows up to and including any line starting with "/* ";
 *     we read until the body keywords ("*  YYYY MM DD ..." or "P<id>")
 *     start appearing -- the converter does not validate the header.
 *   - An epoch row is exactly "*  YYYY MM DD hh mm ss.ssssssss" with
 *     fixed columns per the spec; we use sscanf with 6 fields.
 *   - Each satellite position row is "P<id> x y z clock_bias  [flags]"
 *     where <id> is a 3-char identifier (G11, E03, R23, ...). x, y, z
 *     are in km (ITRF); clock_bias in microseconds (ignored here).
 *   - We stop at "EOF" or end-of-file (whichever comes first).
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_sp3.h"
#include "spody_eop.h"
#include "spody_earth_orientation.h"
#include "spody_forcemodels.h"   /* ForceModelContext */

/* Static offsets carried by the SpOdy SPDYOUT_ format. Kept in sync
 * with spody/src/sim_run.c -- see comment block at the top of the
 * file. The magic is the ASCII literal, no NUL terminator. */
#define SPODY_SP3_OUT_MAGIC       "SPDYOUT_"
#define SPODY_SP3_OUT_VERSION     1u
#define SPODY_SP3_OUT_STATE_DIM   6u

/* GPS time -> TT offset in seconds. GPS time is locked to TAI - 19 s
 * since the GPS epoch (1980-01-06 UTC = TAI 1980-01-06 + 19 s); TT
 * runs TAI + 32.184 s. So TT = GPS + 19 + 32.184 = GPS + 51.184 s.
 * No leap-second table required: this offset is constant by
 * definition of the GPS time scale. */
#define SPODY_GPS_TO_TT_SEC       51.184

/* TT epoch for the J2000.0 reference: JD 2451545.0 TT (= 2000-01-01
 * 12:00 TT). spody-core's `et` is TDB seconds past this epoch; TDB-TT
 * is sub-millisecond and negligible for orbital validation. */
#define SPODY_JD_J2000_TT         2451545.0

/* SP3 line buffer. The spec caps records at 80 columns but real-world
 * IGS files sometimes carry trailing whitespace and the trailing
 * accuracy/exclude flag block, so 256 is comfortably above. */
#define SPODY_SP3_LINE            256

/* Gregorian (UTC-like) date -> Julian Day at midnight + fractional
 * day. Standard Meeus algorithm (Astronomical Algorithms, ch. 7);
 * valid for any date in the Gregorian calendar (post 1582-10-15)
 * which more than covers any GPS epoch (post 1980). The y/m/d are
 * the wall-clock fields verbatim from the SP3 epoch row. */
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

/* Matrix-vector multiply for the body-fixed to ICRF rotation. The
 * SP3 position is in ITRF (Earth body-fixed for our purposes); we
 * rotate it into ICRF for the output binary. Inlined here rather
 * than pulling spody_math just for one mat*vec. */
static void _rot3_apply(const double R[3][3], const double v[3], double out[3]) {
    out[0] = R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2];
    out[1] = R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2];
    out[2] = R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2];
}

/* Write the 24-byte SPDYOUT_ preamble. Returns 0 on success, non-zero
 * on a short write (which on a fresh fopen("wb") essentially never
 * happens, but we still check for diagnostic clarity). */
static int _write_sp3_out_header(FILE *fp) {
    if (fwrite(SPODY_SP3_OUT_MAGIC, 1, 8, fp) != 8) return -1;
    uint32_t hdr[4] = {
        SPODY_SP3_OUT_VERSION,
        SPODY_SP3_OUT_STATE_DIM,
        0u,
        0u
    };
    if (fwrite(hdr, sizeof(uint32_t), 4, fp) != 4) return -1;
    return 0;
}

int spody_convert_sp3_to_state_icrf(const char *input_sp3,
                                    const char *output_bin,
                                    const char *sat_id,
                                    const char *eop_file,
                                    const char *iau2006_dir) {
    if (!input_sp3 || !output_bin || !sat_id || !eop_file || !iau2006_dir) {
        fprintf(stderr, "sp3: NULL argument\n");
        return 1;
    }
    if (strlen(sat_id) != 3) {
        fprintf(stderr,
            "sp3: sat_id must be 3 chars (got '%s', length %zu)\n",
            sat_id, strlen(sat_id));
        return 1;
    }

    /* --- Bring up the EOP + IAU 2006 machinery one-shot --------- */
    MappedEOPData eop_data = {0};
    if (spody_setup_MappedEOPData(&eop_data, eop_file) != 0) {
        fprintf(stderr, "sp3: cannot load EOP from '%s'\n", eop_file);
        return 1;
    }
    MappedIAU2006Data iau_data = {0};
    if (spody_setup_MappedIAU2006Data(&iau_data, iau2006_dir) != 0) {
        fprintf(stderr, "sp3: cannot load IAU 2006 tables from '%s'\n",
                iau2006_dir);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    MappedEOP     eop_map = {0};
    MappedIAU2006 iau_map = {0};
    spody_setup_MappedEOP(&eop_map, &eop_data);
    spody_setup_MappedIAU2006(&iau_map, &iau_data);

    /* Minimal ForceModelContext for spody_bf_rotation_earth. The
     * earth-orientation routine reads ONLY ctx->eop and ctx->iau2006;
     * everything else (hg, eph, third bodies, sat, ...) is irrelevant.
     * Stack-allocated, zero-initialised; freed at function exit. */
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    /* --- Open IO ------------------------------------------------- */
    FILE *fin = fopen(input_sp3, "r");
    if (!fin) {
        fprintf(stderr, "sp3: cannot open input '%s'\n", input_sp3);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        fprintf(stderr, "sp3: cannot open output '%s'\n", output_bin);
        fclose(fin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    if (_write_sp3_out_header(fout) != 0) {
        fprintf(stderr, "sp3: cannot write output header\n");
        fclose(fout); fclose(fin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    /* --- Sweep the SP3 body ------------------------------------- */
    /* The parser is intentionally permissive: any line not matching
     * the epoch or position pattern is silently ignored. Real-world
     * IGS files have header rows, comment rows ("/* "), correlation
     * blocks, the trailing "EOF" sentinel; none of those need any
     * special handling because they don't match the regex above. */
    char line[SPODY_SP3_LINE];
    double cur_et = 0.0;
    int    have_epoch = 0;
    size_t n_records = 0;
    double et_first = 0.0, et_last = 0.0;

    /* The position row pattern is "P<id> x y z clock ..." where <id>
     * is exactly 3 chars and starts in column 1 (0-indexed). We do a
     * 4-byte prefix compare ("P<id>") rather than scan with %s so a
     * different satellite's record never accidentally matches. */
    char pos_prefix[5];
    snprintf(pos_prefix, sizeof pos_prefix, "P%s", sat_id);
    size_t pos_prefix_len = strlen(pos_prefix);   /* always 4 */

    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '*' && line[1] == ' ') {
            int    yy, mm, dd, hh, mn;
            double ss;
            /* SP3 epoch line layout (SP3-d sec. 3.2.4):
             *   "*  YYYY MM DD HH MM SS.SSSSSSSS"
             * Column widths are exact, but sscanf handles
             * whitespace-tolerant input. We accept any spacing. */
            if (sscanf(line, "* %d %d %d %d %d %lf",
                       &yy, &mm, &dd, &hh, &mn, &ss) != 6) {
                continue;
            }
            double jd_gps = _greg_to_jd(yy, mm, dd, hh, mn, ss);
            double jd_tt  = jd_gps + SPODY_GPS_TO_TT_SEC / 86400.0;
            cur_et        = (jd_tt - SPODY_JD_J2000_TT) * 86400.0;
            have_epoch    = 1;
            continue;
        }
        if (!have_epoch) continue;
        if (strncmp(line, pos_prefix, pos_prefix_len) != 0) continue;

        /* Position record. Per SP3-d sec. 3.2.5 the four numeric
         * fields are x [km], y [km], z [km], clock_bias [us]. We
         * read with %lf; the 'V' velocity record (rare) would have
         * leading 'V' and is filtered out by the prefix match. */
        double x_itrf, y_itrf, z_itrf, clk;
        if (sscanf(line + pos_prefix_len, "%lf %lf %lf %lf",
                   &x_itrf, &y_itrf, &z_itrf, &clk) < 3) {
            continue;
        }
        (void)clk;   /* clock bias not used in the position binary */

        /* Bad-position sentinel: SP3 uses 0.000000 as the "no data"
         * marker per spec sec. 3.2.5. Skip cleanly to avoid emitting
         * an origin-of-frame artefact that diff tools would flag. */
        if (x_itrf == 0.0 && y_itrf == 0.0 && z_itrf == 0.0) continue;

        double R_i2bf[3][3], R_bf2i[3][3];
        spody_bf_rotation_earth(&ctx, cur_et, R_i2bf, R_bf2i);
        double pos_itrf[3] = { x_itrf, y_itrf, z_itrf };
        double pos_icrf[3];
        _rot3_apply(R_bf2i, pos_itrf, pos_icrf);

        /* Time column = integrator's 0-based t, matching sim_run.c
         * emit_trajectory. See the same fix on the GLONASS path; the
         * GUI analysis diff aligns by this column. */
        if (n_records == 0) et_first = cur_et;
        double rec[7] = {
            cur_et - et_first,
            pos_icrf[0], pos_icrf[1], pos_icrf[2],
            0.0, 0.0, 0.0
        };
        if (fwrite(rec, sizeof(double), 7, fout) != 7) {
            fprintf(stderr,
                "sp3: short write at record %zu (et=%.6f)\n",
                n_records, cur_et);
            fclose(fout); fclose(fin);
            spody_free_MappedIAU2006(&iau_map);
            spody_free_MappedEOP(&eop_map);
            spody_free_MappedIAU2006Data(&iau_data);
            spody_free_MappedEOPData(&eop_data);
            return 1;
        }
        et_last = cur_et;
        ++n_records;
    }
    fclose(fout);
    fclose(fin);

    /* Operator-visible summary. Stderr so stdout stays clean for any
     * downstream piping (consistent with the other 'spody convert *'
     * subcommands). */
    if (n_records == 0) {
        fprintf(stderr,
            "sp3: WARNING -- no records written for sat_id '%s' "
            "(no matching P%s row found in '%s')\n",
            sat_id, sat_id, input_sp3);
    } else {
        double duration_h = (et_last - et_first) / 3600.0;
        fprintf(stderr,
            "sp3: '%s' -> %zu records (sat=%s, et=%.6f..%.6f, %.3f h)\n",
            input_sp3, n_records, sat_id, et_first, et_last, duration_h);
    }

    spody_free_MappedIAU2006(&iau_map);
    spody_free_MappedEOP(&eop_map);
    spody_free_MappedIAU2006Data(&iau_data);
    spody_free_MappedEOPData(&eop_data);
    return 0;
}
