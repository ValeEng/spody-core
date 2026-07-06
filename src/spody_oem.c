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
 * CCSDS OEM -> SpOdy SPDYOUT_ converter. See spody_oem.h for the
 * supported subset and the wire format.
 *
 * Parse strategy: a line-oriented sweep with three modes --
 *   1. metadata (between META_START and META_STOP): collect REF_FRAME
 *      and TIME_SYSTEM, ignore everything else;
 *   2. covariance (between COVARIANCE_START and COVARIANCE_STOP):
 *      skipped entirely;
 *   3. data: any line starting with a digit is an ephemeris row
 *      "epoch x y z vx vy vz [ax ay az]"; COMMENT / blank lines and
 *      the header keywords (CCSDS_OEM_VERS, CREATION_DATE, ...) are
 *      passed over.
 * A file may carry several META blocks (multi-segment OEM); each one
 * re-validates the frame / time system before its data rows.
 */
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "spody_oem.h"
#include "spody_const.h"
#include "spody_time.h"

/* Static offsets carried by the SpOdy SPDYOUT_ format. Kept in sync
 * with spody/src/sim_run.c -- see the comment block in spody_oem.h.
 * The magic is the ASCII literal, no NUL terminator. */
#define SPODY_OEM_OUT_MAGIC       "SPDYOUT_"
#define SPODY_OEM_OUT_VERSION     1u
#define SPODY_OEM_OUT_STATE_DIM   6u

/* OEM line buffer. The spec (CCSDS 502.0-B-3, sect. 7.6) caps lines at
 * 254 characters. */
#define SPODY_OEM_LINE            512

/* TIME_SYSTEM values we accept, as a discriminator. */
typedef enum {
    OEM_TS_NONE = 0,
    OEM_TS_UTC,
    OEM_TS_TDB
} OemTimeSystem;

/* Write the 24-byte SPDYOUT_ preamble. Returns 0 on success, non-zero
 * on a short write. */
static int write_oem_out_header(FILE *fp) {
    if (fwrite(SPODY_OEM_OUT_MAGIC, 1, 8, fp) != 8) return -1;
    uint32_t hdr[4] = {
        SPODY_OEM_OUT_VERSION,
        SPODY_OEM_OUT_STATE_DIM,
        0u,
        0u
    };
    if (fwrite(hdr, sizeof(uint32_t), 4, fp) != 4) return -1;
    return 0;
}

/* Extract the value token of a "KEY = VALUE" metadata line into
 * out[out_sz], trimming surrounding whitespace and the trailing
 * newline. Returns 1 when a value was found, 0 otherwise. */
static int oem_meta_value(const char *line, char *out, size_t out_sz) {
    const char *eq = strchr(line, '=');
    if (!eq) return 0;
    const char *p = eq + 1;
    while (*p == ' ' || *p == '\t') ++p;
    size_t n = 0;
    while (p[n] && p[n] != ' ' && p[n] != '\t' &&
           p[n] != '\r' && p[n] != '\n' && n + 1 < out_sz) {
        out[n] = p[n];
        ++n;
    }
    out[n] = '\0';
    return n > 0;
}

/* Scan one open OEM file and append SPDYOUT_ state records to *fout.
 * Cross-file accumulators (et_first, et_last, counters) live in the
 * caller so the time column stays 0-anchored on the very first record
 * across all inputs and overlap-skipping works across file
 * boundaries. Returns 0 on success, non-zero on parse / write
 * failure. */
static int oem_scan_file(FILE *fin,
                         FILE *fout,
                         const char *input_oem,
                         size_t *n_written_inout,
                         size_t *n_skipped_inout,
                         double *et_first_inout,
                         double *et_last_inout) {
    char line[SPODY_OEM_LINE];

    char          ref_frame[32] = {0};
    OemTimeSystem time_system   = OEM_TS_NONE;
    int           in_meta       = 0;
    int           in_covariance = 0;
    int           frame_checked = 0;
    long          line_no       = 0;

    while (fgets(line, sizeof line, fin)) {
        ++line_no;

        if (in_covariance) {
            if (strncmp(line, "COVARIANCE_STOP", 15) == 0) in_covariance = 0;
            continue;
        }
        if (strncmp(line, "COVARIANCE_START", 16) == 0) {
            in_covariance = 1;
            continue;
        }
        if (strncmp(line, "META_START", 10) == 0) {
            in_meta        = 1;
            ref_frame[0]   = '\0';
            time_system    = OEM_TS_NONE;
            frame_checked  = 0;
            continue;
        }
        if (strncmp(line, "META_STOP", 9) == 0) {
            in_meta = 0;
            continue;
        }
        if (in_meta) {
            char value[32];
            if (strncmp(line, "REF_FRAME", 9) == 0 &&
                oem_meta_value(line, value, sizeof value)) {
                snprintf(ref_frame, sizeof ref_frame, "%s", value);
            } else if (strncmp(line, "TIME_SYSTEM", 11) == 0 &&
                       oem_meta_value(line, value, sizeof value)) {
                if      (strcmp(value, "UTC") == 0) time_system = OEM_TS_UTC;
                else if (strcmp(value, "TDB") == 0) time_system = OEM_TS_TDB;
                else {
                    fprintf(stderr,
                        "oem: '%s' line %ld: unsupported TIME_SYSTEM '%s' "
                        "(supported: UTC, TDB)\n",
                        input_oem, line_no, value);
                    return 1;
                }
            }
            continue;
        }

        /* Data section: ephemeris rows start with a year digit;
         * everything else (COMMENT, header keywords, blanks) is
         * passed over. */
        if (!isdigit((unsigned char)line[0])) continue;

        /* Validate the segment's metadata on its first data row, not
         * at META_STOP, so files listing keys we do not care about
         * still parse and the error points at real data. */
        if (!frame_checked) {
            if (time_system == OEM_TS_NONE) {
                fprintf(stderr,
                    "oem: '%s' line %ld: data row before a TIME_SYSTEM "
                    "declaration\n", input_oem, line_no);
                return 1;
            }
            if (strcmp(ref_frame, "ICRF")    != 0 &&
                strcmp(ref_frame, "EME2000") != 0 &&
                strcmp(ref_frame, "J2000")   != 0) {
                fprintf(stderr,
                    "oem: '%s' line %ld: unsupported REF_FRAME '%s' "
                    "(supported: ICRF, EME2000, J2000)\n",
                    input_oem, line_no,
                    ref_frame[0] ? ref_frame : "(none)");
                return 1;
            }
            frame_checked = 1;
        }

        int    y, mo, d, hh, mn;
        double ss, state[6];
        int n_fields = sscanf(line,
                "%d-%d-%dT%d:%d:%lf %lf %lf %lf %lf %lf %lf",
                &y, &mo, &d, &hh, &mn, &ss,
                &state[0], &state[1], &state[2],
                &state[3], &state[4], &state[5]);
        if (n_fields != 12) {
            fprintf(stderr,
                "oem: '%s' line %ld: malformed ephemeris row "
                "(parsed %d of 12 fields)\n", input_oem, line_no, n_fields);
            return 1;
        }

        /* Epoch -> ET (s past J2000 TDB). The date's midnight JD is
         * exact in double and (jd0 - JD_J2000) is a small half-integer,
         * so building seconds-past-J2000 as day-difference * 86400 +
         * seconds-of-day preserves the full resolution of the text
         * field (a full-magnitude JD would quantise epochs at ~40 us,
         * i.e. ~30 cm along a LEO track). UTC needs the leap-second
         * chain UTC -> TAI -> TT plus the TDB periodic term; TDB is
         * ET by definition. */
        double jd0      = spody_greg_to_jd(y, mo, d, 0, 0, 0.0);
        double base_sec = (jd0 - JD_J2000) * SECONDSxDAY
                        + hh * 3600.0 + mn * 60.0 + ss;
        double et;
        if (time_system == OEM_TS_UTC) {
            double tt_sec = base_sec
                          + spody_tai_minus_utc(jd0 - JD_MJD_EPOCH)
                          - TT2TAI_SEC;
            et = tt_sec + spody_tdb_minus_tt(tt_sec);
        } else {
            et = base_sec;
        }

        /* Overlap / duplicate guard: consecutive daily OEM releases
         * overlap in span, so a record that does not advance past the
         * last written epoch is dropped (first-file-wins). */
        if (*n_written_inout > 0 && et <= *et_last_inout) {
            ++(*n_skipped_inout);
            continue;
        }

        /* 0-anchored time column, same contract as sim_run.c's
         * emit_trajectory: downstream tooling aligns propagator and
         * reference by this column; the absolute epoch is carried by
         * [simulation].et_start_s in the companion TOML. */
        if (*n_written_inout == 0) *et_first_inout = et;
        double rec[7] = {
            et - *et_first_inout,
            state[0], state[1], state[2],
            state[3], state[4], state[5]
        };
        if (fwrite(rec, sizeof(double), 7, fout) != 7) {
            fprintf(stderr, "oem: short write at record %zu (et=%.6f)\n",
                    *n_written_inout, et);
            return 1;
        }
        *et_last_inout = et;
        ++(*n_written_inout);
    }

    return 0;
}


int spody_convert_oem_to_state_icrf(int n_inputs,
                                    const char *const *input_oem_paths,
                                    const char *output_bin) {
    if (n_inputs <= 0 || !input_oem_paths || !output_bin) {
        fprintf(stderr, "oem: NULL argument or empty input list\n");
        return 1;
    }
    for (int i = 0; i < n_inputs; ++i) {
        if (!input_oem_paths[i]) {
            fprintf(stderr, "oem: NULL input path at index %d\n", i);
            return 1;
        }
    }

    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        fprintf(stderr, "oem: cannot open output '%s'\n", output_bin);
        return 1;
    }
    if (write_oem_out_header(fout) != 0) {
        fprintf(stderr, "oem: cannot write output header\n");
        fclose(fout);
        return 1;
    }

    size_t n_written = 0;
    size_t n_skipped = 0;
    double et_first  = 0.0;
    double et_last   = 0.0;
    int    rc        = 0;

    for (int i = 0; i < n_inputs; ++i) {
        const char *input_oem = input_oem_paths[i];
        FILE *fin = fopen(input_oem, "r");
        if (!fin) {
            fprintf(stderr, "oem: cannot open input '%s'\n", input_oem);
            rc = 1;
            break;
        }
        int file_rc = oem_scan_file(fin, fout, input_oem,
                                    &n_written, &n_skipped,
                                    &et_first, &et_last);
        fclose(fin);
        if (file_rc != 0) { rc = file_rc; break; }
    }

    fclose(fout);

    if (rc == 0) {
        if (n_written == 0) {
            fprintf(stderr,
                "oem: WARNING -- no ephemeris records found across %d "
                "input file%s\n", n_inputs, n_inputs == 1 ? "" : "s");
        } else {
            double duration_h = (et_last - et_first) / 3600.0;
            fprintf(stderr,
                "oem: %zu records (%d file%s, et=%.6f..%.6f, %.3f h, "
                "%zu overlapping record%s skipped)\n",
                n_written, n_inputs, n_inputs == 1 ? "" : "s",
                et_first, et_last, duration_h,
                n_skipped, n_skipped == 1 ? "" : "s");
        }
    }
    return rc;
}
