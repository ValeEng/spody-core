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
 *   - Header: the version letter (line 1, "#a".."#d") and the time
 *     system of the first "%c" line (columns 10-12: GPS, GAL, QZS,
 *     IRN, BDT, TAI, UTC, GLO). SP3-a/-b predate that field and are
 *     GPS time by definition; any other value is refused, never
 *     guessed (a UTC file read as GPS lands 18 s off: ~100 km along a
 *     LAGEOS track). The rest of the header is not validated.
 *   - An epoch row is exactly "*  YYYY MM DD hh mm ss.ssssssss" with
 *     fixed columns per the spec; we use sscanf with 6 fields.
 *   - Each satellite position row is "P<id> x y z clock_bias  [flags]"
 *     where <id> is a 3-char identifier (G11, E03, R23, ...). x, y, z
 *     are in km (ITRF); clock_bias in microseconds (ignored here).
 *   - We stop at "EOF" or end-of-file (whichever comes first).
 */
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_sp3.h"
#include "spody_eop.h"
#include "spody_earth_orientation.h"
#include "spody_forcemodels.h"   /* ForceModelContext */
#include "spody_const.h"
#include "spody_math.h"
#include "spody_time.h"
#include "spody_io.h"

/* Static offsets carried by the SpOdy SPDYOUT_ format. Kept in sync
 * with spody/src/sim_run.c -- see comment block at the top of the
 * file. The magic is the ASCII literal, no NUL terminator. */
#define SPODY_SP3_OUT_MAGIC       "SPDYOUT_"
#define SPODY_SP3_OUT_VERSION     1u
#define SPODY_SP3_OUT_STATE_DIM   6u

/* SP3 line buffer. The spec caps records at 80 columns but real-world
 * IGS files sometimes carry trailing whitespace and the trailing
 * accuracy/exclude flag block, so 256 is comfortably above. */
#define SPODY_SP3_LINE            256

/* Time scales an SP3 epoch can be written in, after folding the ones
 * that share GPS time's epoch and rate (GAL, QZS, IRN) into GPS. */
typedef enum {
    SP3_TS_GPS = 0,
    SP3_TS_BDT,
    SP3_TS_TAI,
    SP3_TS_UTC,
    SP3_TS_GLO
} Sp3TimeScale;

/* Resolve the "%c" time-system field; -1 when it is not one we can
 * convert exactly. */
static int sp3_time_scale(char version, const char ts[4]) {
    if (!strcmp(ts, "GPS") || !strcmp(ts, "GAL") ||
        !strcmp(ts, "QZS") || !strcmp(ts, "IRN")) return SP3_TS_GPS;
    if (!strcmp(ts, "BDT")) return SP3_TS_BDT;
    if (!strcmp(ts, "TAI")) return SP3_TS_TAI;
    if (!strcmp(ts, "UTC")) return SP3_TS_UTC;
    if (!strcmp(ts, "GLO")) return SP3_TS_GLO;
    /* SP3-a and SP3-b have no time-system field (a "ccc" placeholder
     * or no %c line at all): their epochs are GPS time. */
    if ((version == 'a' || version == 'b') &&
        (ts[0] == '\0' || !strcmp(ts, "ccc"))) return SP3_TS_GPS;
    return -1;
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

/* Scan one open SP3 file and append SPDYOUT_ position records for
 * *sat_id* to *fout*. Cross-file accumulators (et_first, et_last,
 * n_records) live in the caller so the time column stays 0-anchored
 * to the very first record across all inputs and the multi-file
 * summary is correct. Returns 0 on success, non-zero on parse /
 * write failure. */
static int _sp3_scan_file(FILE *fin,
                          FILE *fout,
                          const char *input_sp3,
                          const char *sat_id,
                          const ForceModelContext *ctx,
                          size_t *n_records_inout,
                          double *et_first_inout,
                          double *et_last_inout) {
    char line[SPODY_SP3_LINE];
    double cur_et = 0.0;
    int    have_epoch = 0;
    char   version = 0;          /* 'a'..'d', from line 1 */
    char   ts_field[4] = "";     /* first %c line, columns 10-12 */
    int    scale = -1;           /* resolved at the first epoch row */
    size_t n_records_this = 0;
    double et_first_this  = 0.0;
    double et_last_this   = 0.0;

    /* "P<id>" 4-byte prefix (always 4 chars; sat_id is enforced to 3
     * by the public entry point). */
    char pos_prefix[5];
    snprintf(pos_prefix, sizeof pos_prefix, "P%s", sat_id);
    size_t pos_prefix_len = strlen(pos_prefix);

    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '#' && line[1] != '#' && version == 0) {
            version = line[1];
            continue;
        }
        if (line[0] == '%' && line[1] == 'c' && ts_field[0] == '\0') {
            if (strlen(line) >= 12) {
                memcpy(ts_field, line + 9, 3);
                ts_field[3] = '\0';
            }
            continue;
        }
        if (line[0] == '*' && line[1] == ' ') {
            int    yy, mm, dd, hh, mn;
            double ss;
            if (sscanf(line, "* %d %d %d %d %d %lf",
                       &yy, &mm, &dd, &hh, &mn, &ss) != 6) {
                continue;
            }
            if (scale < 0) {
                scale = sp3_time_scale(version, ts_field);
                if (scale < 0) {
                    spody_log_eprintf(
                           "sp3: '%s': time system '%s' (SP3-%c) is not "
                           "supported (GPS, GAL, QZS, IRN, BDT, TAI, UTC, "
                           "GLO); file refused\n",
                           input_sp3, ts_field, version ? version : '?');
                    return 1;
                }
            }
            /* Epoch -> TT, exactly, per scale: GPS-aligned scales are
             * TT - 51.184 s, BDT a further 14 s behind, TAI TT - 32.184
             * s; UTC needs the leap offset of its calendar day (right
             * for a 23:59:60 epoch too), and GLONASS time is UTC + 3 h,
             * so its day is the UTC one three hours earlier. Then
             * TT -> TDB adds the deltet periodic term (+/-1.657 ms). */
            double base = spody_greg_to_sec_j2000(yy, mm, dd, hh, mn, ss);
            double tt_sec;
            if (scale == SP3_TS_GPS) {
                tt_sec = base + GPST2TT_SEC;
            } else if (scale == SP3_TS_BDT) {
                tt_sec = base + BDT2GPST_SEC + GPST2TT_SEC;
            } else if (scale == SP3_TS_TAI) {
                tt_sec = base - TT2TAI_SEC;
            } else {
                double jd_day = spody_greg_to_jd(yy, mm, dd, 0, 0, 0.0);
                if (scale == SP3_TS_GLO) {
                    base -= GLOT_MINUS_UTC_SEC;
                    if (hh < 3) jd_day -= 1.0;
                }
                tt_sec = base + spody_tai_minus_utc(jd_day - JD_MJD_EPOCH)
                       - TT2TAI_SEC;
            }
            cur_et        = tt_sec + spody_tdb_minus_tt(tt_sec);
            have_epoch    = 1;
            continue;
        }
        if (!have_epoch) continue;
        if (strncmp(line, pos_prefix, pos_prefix_len) != 0) continue;

        double x_itrf, y_itrf, z_itrf, clk;
        if (sscanf(line + pos_prefix_len, "%lf %lf %lf %lf",
                   &x_itrf, &y_itrf, &z_itrf, &clk) < 3) {
            continue;
        }
        (void)clk;

        /* SP3 0.000000 = "no data" sentinel (SP3-d sect. 3.2.5). */
        if (x_itrf == 0.0 && y_itrf == 0.0 && z_itrf == 0.0) continue;

        /* Outside the EOP table the rotation falls back to the
         * identity: refuse instead of writing ITRF labelled ICRF. */
        const MappedEOPData *eop = ctx->eop->med;
        double mjd_utc = spody_et_to_mjd_utc(cur_et);
        if (!spody_eop_covers_mjd(eop, mjd_utc)) {
            spody_log_eprintf(
                   "sp3: epoch UTC MJD %.5f in '%s' is outside the EOP table "
                   "(UTC MJD %.2f .. %.2f); update finals2000A.all\n",
                   mjd_utc, input_sp3, eop->mjd_first, eop->mjd_last_predicted);
            return 1;
        }

        double R_i2bf[3][3], R_bf2i[3][3];
        spody_bf_rotation_earth(ctx, cur_et, R_i2bf, R_bf2i);
        double pos_itrf[3] = { x_itrf, y_itrf, z_itrf };
        double pos_icrf[3];
        spody_rotate_vector(R_bf2i, pos_itrf, pos_icrf);

        /* Time column = integrator's 0-based t. Multi-file mode:
         * et_first anchored on the FIRST written record across ALL
         * input files via the caller-owned *et_first_inout. */
        if (*n_records_inout == 0) *et_first_inout = cur_et;
        if (n_records_this == 0)   et_first_this   = cur_et;
        double rec[7] = {
            cur_et - *et_first_inout,
            pos_icrf[0], pos_icrf[1], pos_icrf[2],
            0.0, 0.0, 0.0
        };
        if (fwrite(rec, sizeof(double), 7, fout) != 7) {
            spody_log_eprintf(
                   "sp3: short write at record %zu (et=%.6f)\n",
                   *n_records_inout, cur_et);
            return 1;
        }
        *et_last_inout = cur_et;
        et_last_this   = cur_et;
        ++(*n_records_inout);
        ++n_records_this;
    }

    if (n_records_this == 0) {
        spody_log_eprintf(
               "sp3: WARNING -- no records written for sat_id '%s' "
               "(no matching P%s row found in '%s')\n",
               sat_id, sat_id, input_sp3);
    } else {
        double duration_h = (et_last_this - et_first_this) / 3600.0;
        spody_log_eprintf(
               "sp3: '%s' (time %s) -> %zu records (sat=%s, et=%.6f..%.6f, %.3f h)\n",
               input_sp3, ts_field[0] ? ts_field : "GPS", n_records_this, sat_id,
               et_first_this, et_last_this, duration_h);
    }
    return 0;
}


int spody_convert_sp3_to_state_icrf(int n_inputs,
                                    const char *const *input_sp3_paths,
                                    const char *output_bin,
                                    const char *sat_id,
                                    const char *eop_file,
                                    const char *iau2006_dir) {
    if (n_inputs <= 0 || !input_sp3_paths || !output_bin ||
        !sat_id || !eop_file || !iau2006_dir) {
        spody_log_eprintf("sp3: NULL argument or empty input list\n");
        return 1;
    }
    for (int i = 0; i < n_inputs; ++i) {
        if (!input_sp3_paths[i]) {
            spody_log_eprintf("sp3: NULL input path at index %d\n", i);
            return 1;
        }
    }
    if (strlen(sat_id) != 3) {
        spody_log_eprintf(
               "sp3: sat_id must be 3 chars (got '%s', length %zu)\n",
               sat_id, strlen(sat_id));
        return 1;
    }

    /* --- Bring up the EOP + IAU 2006 machinery one-shot --------- */
    MappedEOPData eop_data = {0};
    if (spody_setup_MappedEOPData(&eop_data, eop_file) != 0) {
        spody_log_eprintf("sp3: cannot load EOP from '%s'\n", eop_file);
        return 1;
    }
    MappedIAU2006Data iau_data = {0};
    if (spody_setup_MappedIAU2006Data(&iau_data, iau2006_dir) != 0) {
        spody_log_eprintf("sp3: cannot load IAU 2006 tables from '%s'\n",
                  iau2006_dir);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    MappedEOP     eop_map = {0};
    MappedIAU2006 iau_map = {0};
    spody_setup_MappedEOP(&eop_map, &eop_data);
    spody_setup_MappedIAU2006(&iau_map, &iau_data);
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    /* --- Output (one binary for the whole concatenated track) --- */
    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        spody_log_eprintf("sp3: cannot open output '%s'\n", output_bin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    if (_write_sp3_out_header(fout) != 0) {
        spody_log_eprintf("sp3: cannot write output header\n");
        fclose(fout);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    /* --- Loop over input files; share et_first across all of them */
    size_t n_records_all = 0;
    double et_first_all  = 0.0;
    double et_last_all   = 0.0;
    int    rc            = 0;

    for (int i = 0; i < n_inputs; ++i) {
        const char *input_sp3 = input_sp3_paths[i];
        FILE *fin = fopen(input_sp3, "r");
        if (!fin) {
            spody_log_eprintf("sp3: cannot open input '%s'\n", input_sp3);
            rc = 1;
            break;
        }
        int file_rc = _sp3_scan_file(fin, fout, input_sp3, sat_id, &ctx,
                                      &n_records_all,
                                      &et_first_all, &et_last_all);
        fclose(fin);
        if (file_rc != 0) { rc = file_rc; break; }
    }

    /* Records are buffered: a full disk surfaces as the stream error
     * flag or at the final flush in fclose, not at the fwrite calls. */
    int write_failed = ferror(fout);
    if (fclose(fout) != 0 || write_failed) {
        spody_log_eprintf("sp3: write failed on '%s': %s\n",
                  output_bin, strerror(errno));
        rc = 1;
    }

    if (rc == 0 && n_inputs > 1 && n_records_all > 0) {
        double duration_h = (et_last_all - et_first_all) / 3600.0;
        spody_log_eprintf(
               "sp3: aggregate -> %zu records across %d files "
               "(sat=%s, et=%.6f..%.6f, %.3f h)\n",
               n_records_all, n_inputs, sat_id,
               et_first_all, et_last_all, duration_h);
    }

    spody_free_MappedIAU2006(&iau_map);
    spody_free_MappedEOP(&eop_map);
    spody_free_MappedIAU2006Data(&iau_data);
    spody_free_MappedEOPData(&eop_data);
    return rc;
}
