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
#include "spody_const.h"
#include "spody_math.h"
#include "spody_time.h"

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
    size_t n_records_this = 0;
    double et_first_this  = 0.0;
    double et_last_this   = 0.0;

    /* "P<id>" 4-byte prefix (always 4 chars; sat_id is enforced to 3
     * by the public entry point). */
    char pos_prefix[5];
    snprintf(pos_prefix, sizeof pos_prefix, "P%s", sat_id);
    size_t pos_prefix_len = strlen(pos_prefix);

    while (fgets(line, sizeof line, fin)) {
        if (line[0] == '*' && line[1] == ' ') {
            int    yy, mm, dd, hh, mn;
            double ss;
            if (sscanf(line, "* %d %d %d %d %d %lf",
                       &yy, &mm, &dd, &hh, &mn, &ss) != 6) {
                continue;
            }
            /* SP3 epochs are GPST; TT = GPST + 51.184 s exactly, no
             * leap-second table required (GPST2TT_SEC in
             * spody_const.h); TT -> TDB adds the deltet periodic
             * term (+/-1.657 ms). */
            double jd_gps = spody_greg_to_jd(yy, mm, dd, hh, mn, ss);
            double jd_tt  = jd_gps + GPST2TT_SEC / SECONDSxDAY;
            double tt_sec = ET_FROM_JD(jd_tt);
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
            fprintf(stderr,
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
        fprintf(stderr,
            "sp3: WARNING -- no records written for sat_id '%s' "
            "(no matching P%s row found in '%s')\n",
            sat_id, sat_id, input_sp3);
    } else {
        double duration_h = (et_last_this - et_first_this) / 3600.0;
        fprintf(stderr,
            "sp3: '%s' -> %zu records (sat=%s, et=%.6f..%.6f, %.3f h)\n",
            input_sp3, n_records_this, sat_id,
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
        fprintf(stderr, "sp3: NULL argument or empty input list\n");
        return 1;
    }
    for (int i = 0; i < n_inputs; ++i) {
        if (!input_sp3_paths[i]) {
            fprintf(stderr, "sp3: NULL input path at index %d\n", i);
            return 1;
        }
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
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    /* --- Output (one binary for the whole concatenated track) --- */
    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        fprintf(stderr, "sp3: cannot open output '%s'\n", output_bin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    if (_write_sp3_out_header(fout) != 0) {
        fprintf(stderr, "sp3: cannot write output header\n");
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
            fprintf(stderr, "sp3: cannot open input '%s'\n", input_sp3);
            rc = 1;
            break;
        }
        int file_rc = _sp3_scan_file(fin, fout, input_sp3, sat_id, &ctx,
                                      &n_records_all,
                                      &et_first_all, &et_last_all);
        fclose(fin);
        if (file_rc != 0) { rc = file_rc; break; }
    }

    fclose(fout);

    if (rc == 0 && n_inputs > 1 && n_records_all > 0) {
        double duration_h = (et_last_all - et_first_all) / 3600.0;
        fprintf(stderr,
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
