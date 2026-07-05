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
 * Implementation of spody_convert_gps_to_state_icrf. See spody_gps.h
 * for the format and algorithm references (IS-GPS-200 sect. 20.3.3.4.3
 * positions, Remondi 2004 velocities, RINEX 3.05 sect. 6.10.1 nav
 * record layout, GPS-time TT bridge GPST2TT_SEC = 51.184 s exactly).
 *
 * RINEX 3.x GPS nav record layout (one PRN line + seven continuation
 * lines, each numeric field in Fortran D-format which sscanf %lf reads
 * once we swap D->E):
 *
 *   PRN row:  G<NN> YYYY MM DD HH MM SS  clk_af0 clk_af1 clk_af2
 *   L1:       <indent> IODE   Crs   delta_n   M0
 *   L2:       <indent> Cuc    e     Cus       sqrt_A
 *   L3:       <indent> Toe    Cic   Omega0    Cis
 *   L4:       <indent> i0     Crc   omega     OmegaDot
 *   L5:       <indent> iDot   codesL2 GPSweek L2P_flag
 *   L6:       <indent> SVacc  SVhealth  TGD   IODC
 *   L7:       <indent> txTime fitInterval  spare  spare
 *
 * Units (per IS-GPS-200 / RINEX 3.05):
 *   sqrt_A     m^0.5         -> a = sqrt_A^2 in metres
 *   e          dimensionless
 *   M0, omega, Omega0, i0, Cuc, Cus, Cic, Cis  radians
 *   Crs, Crc                 metres
 *   delta_n, iDot, OmegaDot  rad/s
 *   Toe                      seconds of GPS week
 *   IODE                     issue-of-data ephemeris (integer, parsed
 *                            as double; not used for propagation, kept
 *                            for stricter "same arc" sanity if needed)
 *
 * For our purpose we evaluate each broadcast record AT its own TOC
 * (the per-record epoch on the PRN row). For most messages
 * TOC seconds-of-week == Toe exactly, so tk = 0 and the propagation
 * is essentially "what the broadcast says (r, v) is right now". The
 * resulting binary therefore has one record per RINEX nav message
 * (~12/day for GPS, vs ~48/day for GLONASS).
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_gps.h"
#include "spody_eop.h"
#include "spody_earth_orientation.h"
#include "spody_forcemodels.h"
#include "spody_const.h"
#include "spody_math.h"
#include "spody_time.h"

#define SPODY_GPS_OUT_MAGIC          "SPDYOUT_"
#define SPODY_GPS_OUT_VERSION        1u
#define SPODY_GPS_OUT_STATE_DIM      6u

/* RINEX nav lines fit in 80 cols + newline + sentinel; 256 covers
 * any RINEX 3.x version with comfortable slack. */
#define SPODY_RINEX_LINE             256

/* EARTH_MU is in km^3/s^2; convert to m^3/s^2 for the broadcast
 * algorithm which works in metres. */
#define SPODY_EARTH_MU_M3S2          (EARTH_MU * 1.0e9)

/* Fortran-D exponent normaliser for sscanf. RINEX 3 emits 1.23D-04;
 * MSVC's strtod rejects D, so we re-parse via E. Re-used for any line
 * whose first attempt yields exactly 0.0 with a D in the buffer. */
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

/* Reparse one continuation line through `_parse_rinex_double` field-
 * wise when the first sscanf got 0s because of Fortran D exponents.
 * Mutates the inputs in place; returns 1 if 4 fields parsed cleanly,
 * 0 otherwise. */
static int _reparse_d_line(const char *line, double *a, double *b,
                            double *c, double *d) {
    char buf[SPODY_RINEX_LINE];
    size_t i;
    for (i = 0; i + 1 < sizeof(buf) && line[i]; ++i)
        buf[i] = (line[i] == 'D' || line[i] == 'd') ? 'E' : line[i];
    buf[i] = '\0';
    return sscanf(buf, "%lf %lf %lf %lf", a, b, c, d) == 4;
}

/* SPDYOUT_ 24-byte preamble identical to sim_run.c's writer. */
static int _write_gps_out_header(FILE *fp) {
    if (fwrite(SPODY_GPS_OUT_MAGIC, 1, 8, fp) != 8) return -1;
    uint32_t hdr[4] = {
        SPODY_GPS_OUT_VERSION,
        SPODY_GPS_OUT_STATE_DIM,
        0u, 0u
    };
    if (fwrite(hdr, sizeof(uint32_t), 4, fp) != 4) return -1;
    return 0;
}

/* Plain-old-data bundle of the 16 GPS broadcast parameters we need
 * for the propagation, plus the per-record epoch. Filled by
 * `_parse_gps_record` from one RINEX nav message, then passed to
 * `_gps_brdc_propagate`. */
typedef struct {
    double toe;          /* time of ephemeris, seconds of GPS week    */
    double sqrt_a;       /* sqrt(semi-major axis), sqrt(m)            */
    double e;            /* eccentricity                              */
    double i0;           /* inclination at Toe, rad                   */
    double Omega0;       /* RAAN at week start, rad                   */
    double omega;        /* argument of perigee, rad                  */
    double M0;           /* mean anomaly at Toe, rad                  */
    double delta_n;      /* mean-motion correction, rad/s             */
    double OmegaDot;     /* RAAN rate, rad/s                          */
    double idot;         /* inclination rate, rad/s                   */
    double Cuc, Cus;     /* arg-of-latitude harmonic corrections, rad */
    double Crc, Crs;     /* radius harmonic corrections, m            */
    double Cic, Cis;     /* inclination harmonic corrections, rad     */
} _GpsEphemeris;

/* Solve Kepler's equation E - e*sin(E) = M for E by Newton's method.
 * GPS broadcast eccentricities are bounded ~0.02 (GPS orbits are
 * near-circular by design); the iteration converges in 3-4 steps
 * to machine precision. We cap at 12 as a hard safety. */
static double _kepler_solve(double M, double e) {
    double E = M;
    for (int k = 0; k < 12; ++k) {
        double f  = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double dE = f / fp;
        E -= dE;
        if (fabs(dE) < 1.0e-15) break;
    }
    return E;
}

/* GPS broadcast propagation: evaluate (r, v) in WGS-84 ECEF at the
 * given seconds-of-GPS-week `t`. Position formulas are IS-GPS-200
 * sect. 20.3.3.4.3 verbatim; velocity formulas are the closed-form
 * analytic derivatives published by Remondi (2004). Output is in
 * SI metres / metres per second; the caller converts to km / km/s. */
static void _gps_brdc_propagate(const _GpsEphemeris *eph,
                                double t_sec_of_week,
                                double r_ecef[3],
                                double v_ecef[3]) {
    double a   = eph->sqrt_a * eph->sqrt_a;
    double n0  = sqrt(SPODY_EARTH_MU_M3S2 / (a * a * a));
    double n   = n0 + eph->delta_n;

    /* Time from Toe, with GPS-week-rollover correction. */
    double tk = t_sec_of_week - eph->toe;
    if (tk >  HALF_GPS_WEEK_SEC) tk -= GPS_WEEK_SEC;
    if (tk < -HALF_GPS_WEEK_SEC) tk += GPS_WEEK_SEC;

    double M = eph->M0 + n * tk;
    double E = _kepler_solve(M, eph->e);

    double sinE = sin(E), cosE = cos(E);
    double one_minus_ecosE = 1.0 - eph->e * cosE;
    double sqrt_one_e2 = sqrt(1.0 - eph->e * eph->e);
    double nu = atan2(sqrt_one_e2 * sinE, cosE - eph->e);

    double Phi = nu + eph->omega;
    double s2P = sin(2.0 * Phi), c2P = cos(2.0 * Phi);

    /* 2nd-order harmonic corrections. */
    double du = eph->Cus * s2P + eph->Cuc * c2P;
    double dr = eph->Crs * s2P + eph->Crc * c2P;
    double di = eph->Cis * s2P + eph->Cic * c2P;

    double u = Phi + du;
    double r = a * one_minus_ecosE + dr;
    double i = eph->i0 + eph->idot * tk + di;

    double cu = cos(u), su = sin(u);
    double xp = r * cu;
    double yp = r * su;

    /* Corrected RAAN -- referenced to the rotating ECEF frame, so
     * the -omega_earth*toe term anchors the longitude origin at the
     * week-start Greenwich meridian. */
    double Omega = eph->Omega0 + (eph->OmegaDot - EARTH_ROT_RATE_RADPS) * tk
                                - EARTH_ROT_RATE_RADPS * eph->toe;
    double cO = cos(Omega), sO = sin(Omega);
    double ci = cos(i), si = sin(i);

    r_ecef[0] = xp * cO - yp * ci * sO;
    r_ecef[1] = xp * sO + yp * ci * cO;
    r_ecef[2] = yp * si;

    /* --- Velocity: closed-form analytic derivatives (Remondi 2004). */
    double Edot   = n / one_minus_ecosE;
    double nudot  = sqrt_one_e2 * Edot / one_minus_ecosE;
    double Phidot = nudot;

    double dudot  = 2.0 * Phidot * (eph->Cus * c2P - eph->Cuc * s2P);
    double drdot  = 2.0 * Phidot * (eph->Crs * c2P - eph->Crc * s2P);
    double didot  = 2.0 * Phidot * (eph->Cis * c2P - eph->Cic * s2P);

    double udot   = Phidot + dudot;
    double rdot   = a * eph->e * sinE * Edot + drdot;
    double idotk  = eph->idot + didot;
    double Omdot  = eph->OmegaDot - EARTH_ROT_RATE_RADPS;

    double xpdot  = rdot * cu - r * udot * su;
    double ypdot  = rdot * su + r * udot * cu;

    v_ecef[0] = xpdot * cO - ypdot * ci * sO + yp * si * sO * idotk
              - (xp * sO + yp * ci * cO) * Omdot;
    v_ecef[1] = xpdot * sO + ypdot * ci * cO - yp * si * cO * idotk
              + (xp * cO - yp * ci * sO) * Omdot;
    v_ecef[2] = ypdot * si + yp * ci * idotk;
}

/* Parse the 16 broadcast parameters from the 7 continuation lines of
 * one RINEX nav record. Returns 1 on success, 0 if any line fails to
 * yield 4 numeric fields (corrupt record -- caller skips). */
static int _parse_gps_record(const char *l1, const char *l2,
                             const char *l3, const char *l4,
                             const char *l5,
                             _GpsEphemeris *out) {
    double iode, codesL2, week, l2p;
    double a1, b1, c1, d1;

    /* L1: IODE Crs delta_n M0 */
    if (sscanf(l1, "%lf %lf %lf %lf", &iode, &a1, &b1, &c1) != 4) {
        if (!_reparse_d_line(l1, &iode, &a1, &b1, &c1)) return 0;
    }
    (void)iode;
    out->Crs     = a1;
    out->delta_n = b1;
    out->M0      = c1;

    /* L2: Cuc e Cus sqrt_A */
    if (sscanf(l2, "%lf %lf %lf %lf", &a1, &b1, &c1, &d1) != 4) {
        if (!_reparse_d_line(l2, &a1, &b1, &c1, &d1)) return 0;
    }
    out->Cuc    = a1;
    out->e      = b1;
    out->Cus    = c1;
    out->sqrt_a = d1;

    /* L3: Toe Cic Omega0 Cis */
    if (sscanf(l3, "%lf %lf %lf %lf", &a1, &b1, &c1, &d1) != 4) {
        if (!_reparse_d_line(l3, &a1, &b1, &c1, &d1)) return 0;
    }
    out->toe    = a1;
    out->Cic    = b1;
    out->Omega0 = c1;
    out->Cis    = d1;

    /* L4: i0 Crc omega OmegaDot */
    if (sscanf(l4, "%lf %lf %lf %lf", &a1, &b1, &c1, &d1) != 4) {
        if (!_reparse_d_line(l4, &a1, &b1, &c1, &d1)) return 0;
    }
    out->i0       = a1;
    out->Crc      = b1;
    out->omega    = c1;
    out->OmegaDot = d1;

    /* L5: iDot codesL2 GPSweek L2P -- only iDot matters for our use. */
    if (sscanf(l5, "%lf %lf %lf %lf", &a1, &codesL2, &week, &l2p) != 4) {
        if (!_reparse_d_line(l5, &a1, &codesL2, &week, &l2p)) return 0;
    }
    (void)codesL2; (void)week; (void)l2p;
    out->idot = a1;

    return 1;
}

/* GPS seconds-of-week from a Gregorian (UTC-ish) date passed as the
 * RINEX TOC fields. The TOC is GPS time per RINEX 3.05 sect. 6.10.1,
 * so day-of-week + seconds-of-day map directly without a leap-second
 * chain. Returns the seconds-of-week value used by the Kepler-rollover
 * arithmetic in `_gps_brdc_propagate`. */
static double _gpst_seconds_of_week(int y, int m, int d, int hh, int mn,
                                     double ss) {
    /* Day-of-week: Sunday = 0, Saturday = 6. JD has Monday = 0, so
     * we shift by 1.5 days (12 h to JD-midnight + Sunday offset). */
    double jd = spody_greg_to_jd(y, m, d, 0, 0, 0.0);
    int    dow = ((int)floor(jd + 1.5)) % 7;
    return (double)dow * SECONDSxDAY
         + (double)hh  * 3600.0
         + (double)mn  *   60.0
         + ss;
}


/* Scan one open RINEX GPS nav file (header already skipped) and append
 * SPDYOUT_ records for *sat_id* to *fout*. Cross-file accumulators
 * (et_first, et_last, n_written) live in the caller so the time
 * column stays 0-anchored to the very first record across all inputs
 * and the multi-file summary is correct. *n_total_out* receives the
 * count of nav records of any satellite scanned in this file.
 * Returns 0 on success, non-zero on parse / write failure. */
static int _gps_scan_file(FILE *fin,
                          FILE *fout,
                          const char *input_rnx,
                          const char *sat_id,
                          const ForceModelContext *ctx,
                          size_t *n_total_out,
                          size_t *n_written_inout,
                          double *et_first_inout,
                          double *et_last_inout) {
    char line[SPODY_RINEX_LINE];

    /* Skip header. */
    int in_body = 0;
    while (fgets(line, sizeof line, fin)) {
        if (strstr(line, "END OF HEADER")) {
            in_body = 1;
            break;
        }
    }
    if (!in_body) {
        fprintf(stderr,
            "gps: input '%s' has no 'END OF HEADER' marker\n", input_rnx);
        return 1;
    }

    size_t n_total            = 0;
    size_t n_written_this     = 0;
    double et_first_this      = 0.0;
    double et_last_this       = 0.0;

    while (fgets(line, sizeof line, fin)) {
        /* PRN row begins with the 3-char sat id at column 0. */
        if (line[0] != 'G') continue;

        char rec_id[4] = {0};
        memcpy(rec_id, line, 3);
        ++n_total;
        int sat_match = (strcmp(rec_id, sat_id) == 0);

        int    y, mo, d, h, mi;
        double sec;
        double af0, af1, af2;
        /* Line 0 layout after the PRN: " YYYY MM DD HH MM SS  af0 af1 af2" */
        if (sscanf(line + 3, "%d %d %d %d %d %lf %lf %lf %lf",
                   &y, &mo, &d, &h, &mi, &sec, &af0, &af1, &af2) < 6) {
            continue;
        }
        (void)af0; (void)af1; (void)af2;

        /* Read the 7 continuation lines (L1..L7). We consume them
         * unconditionally so the body sweep advances even on
         * non-target satellites. */
        char l1[SPODY_RINEX_LINE], l2[SPODY_RINEX_LINE];
        char l3[SPODY_RINEX_LINE], l4[SPODY_RINEX_LINE];
        char l5[SPODY_RINEX_LINE], l6[SPODY_RINEX_LINE];
        char l7[SPODY_RINEX_LINE];
        if (!fgets(l1, sizeof l1, fin)) break;
        if (!fgets(l2, sizeof l2, fin)) break;
        if (!fgets(l3, sizeof l3, fin)) break;
        if (!fgets(l4, sizeof l4, fin)) break;
        if (!fgets(l5, sizeof l5, fin)) break;
        if (!fgets(l6, sizeof l6, fin)) break;
        if (!fgets(l7, sizeof l7, fin)) break;
        (void)l6; (void)l7;

        if (!sat_match) continue;

        _GpsEphemeris eph;
        if (!_parse_gps_record(l1, l2, l3, l4, l5, &eph)) continue;

        /* Evaluate the broadcast at its own TOC -- gives "what the
         * broadcast says (r, v) is at this epoch" with zero
         * extrapolation. tk in the propagator wraps to ~0 when
         * TOC == Toe. */
        double t_sow = _gpst_seconds_of_week(y, mo, d, h, mi, sec);
        double r_ecef[3], v_ecef[3];
        _gps_brdc_propagate(&eph, t_sow, r_ecef, v_ecef);

        /* GPST -> TT -> ET (TDB) bridge. RINEX TOC is GPST per
         * RINEX 3.05 sect. 6.10.1, and TT = GPST + 51.184 exactly;
         * TT -> TDB adds the deltet periodic term (+/-1.657 ms). */
        double jd_gpst = spody_greg_to_jd(y, mo, d, h, mi, sec);
        double jd_tt   = jd_gpst + GPST2TT_SEC / SECONDSxDAY;
        double tt_sec  = ET_FROM_JD(jd_tt);
        double et      = tt_sec + spody_tdb_minus_tt(tt_sec);

        /* Rotation ECEF -> ICRF (IAU 2006/2000A_R06 + IERS EOP). */
        double R_i2bf[3][3], R_bf2i[3][3];
        spody_bf_rotation_earth(ctx, et, R_i2bf, R_bf2i);

        /* WGS-84 m -> km, ECEF -> ICRF. */
        double r_ecef_km[3] = {
            r_ecef[0] * 1.0e-3, r_ecef[1] * 1.0e-3, r_ecef[2] * 1.0e-3
        };
        double v_ecef_km[3] = {
            v_ecef[0] * 1.0e-3, v_ecef[1] * 1.0e-3, v_ecef[2] * 1.0e-3
        };
        double r_icrf[3], v_icrf_rot[3];
        spody_rotate_vector(R_bf2i, r_ecef_km, r_icrf);
        spody_rotate_vector(R_bf2i, v_ecef_km, v_icrf_rot);

        /* omega x r in ICRF with the full ITRS z-axis (not nominal
         * z-hat). See spody_glonass.c for the J2024 ~480 arcsec
         * precession argument that motivates this. */
        double omega_icrf[3] = {
            EARTH_ROT_RATE_RADPS * R_bf2i[0][2],
            EARTH_ROT_RATE_RADPS * R_bf2i[1][2],
            EARTH_ROT_RATE_RADPS * R_bf2i[2][2],
        };
        double v_icrf[3];
        v_icrf[0] = v_icrf_rot[0]
                  + omega_icrf[1] * r_icrf[2] - omega_icrf[2] * r_icrf[1];
        v_icrf[1] = v_icrf_rot[1]
                  + omega_icrf[2] * r_icrf[0] - omega_icrf[0] * r_icrf[2];
        v_icrf[2] = v_icrf_rot[2]
                  + omega_icrf[0] * r_icrf[1] - omega_icrf[1] * r_icrf[0];

        /* Multi-file mode: et_first anchored on the FIRST written
         * record across ALL input files (caller-owned *et_first_inout). */
        if (*n_written_inout == 0) *et_first_inout = et;
        if (n_written_this == 0)   et_first_this   = et;
        double rec[7] = {
            et - *et_first_inout,
            r_icrf[0], r_icrf[1], r_icrf[2],
            v_icrf[0], v_icrf[1], v_icrf[2]
        };
        if (fwrite(rec, sizeof(double), 7, fout) != 7) {
            fprintf(stderr,
                "gps: short write at record %zu (et=%.6f)\n",
                *n_written_inout, et);
            return 1;
        }
        *et_last_inout = et;
        et_last_this   = et;
        ++(*n_written_inout);
        ++n_written_this;
    }

    if (n_written_this == 0) {
        fprintf(stderr,
            "gps: WARNING -- no records written for sat_id '%s' "
            "(scanned %zu total RINEX records in '%s')\n",
            sat_id, n_total, input_rnx);
    } else {
        double duration_h = (et_last_this - et_first_this) / 3600.0;
        fprintf(stderr,
            "gps: '%s' -> %zu records (sat=%s, et=%.6f..%.6f, "
            "%.3f h, scanned %zu nav messages)\n",
            input_rnx, n_written_this, sat_id,
            et_first_this, et_last_this, duration_h, n_total);
    }

    *n_total_out = n_total;
    return 0;
}


int spody_convert_gps_to_state_icrf(int n_inputs,
                                    const char *const *input_rnx_paths,
                                    const char *output_bin,
                                    const char *sat_id,
                                    const char *eop_file,
                                    const char *iau2006_dir) {
    if (n_inputs <= 0 || !input_rnx_paths || !output_bin ||
        !sat_id || !eop_file || !iau2006_dir) {
        fprintf(stderr, "gps: NULL argument or empty input list\n");
        return 1;
    }
    for (int i = 0; i < n_inputs; ++i) {
        if (!input_rnx_paths[i]) {
            fprintf(stderr, "gps: NULL input path at index %d\n", i);
            return 1;
        }
    }
    if (strlen(sat_id) != 3 || (sat_id[0] != 'G')) {
        fprintf(stderr,
            "gps: sat_id must be 3 chars starting with 'G' "
            "(got '%s')\n", sat_id);
        return 1;
    }

    /* Bring up EOP + IAU 2006 once, share across all input files. */
    MappedEOPData eop_data = {0};
    if (spody_setup_MappedEOPData(&eop_data, eop_file) != 0) {
        fprintf(stderr, "gps: cannot load EOP from '%s'\n", eop_file);
        return 1;
    }
    MappedIAU2006Data iau_data = {0};
    if (spody_setup_MappedIAU2006Data(&iau_data, iau2006_dir) != 0) {
        fprintf(stderr, "gps: cannot load IAU 2006 from '%s'\n",
                iau2006_dir);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    MappedEOP     eop_map = {0};
    MappedIAU2006 iau_map = {0};
    spody_setup_MappedEOP(&eop_map, &eop_data);
    spody_setup_MappedIAU2006(&iau_map, &iau_data);
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    FILE *fout = fopen(output_bin, "wb");
    if (!fout) {
        fprintf(stderr, "gps: cannot open output '%s'\n", output_bin);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    if (_write_gps_out_header(fout) != 0) {
        fprintf(stderr, "gps: cannot write output header\n");
        fclose(fout);
        spody_free_MappedIAU2006(&iau_map);
        spody_free_MappedEOP(&eop_map);
        spody_free_MappedIAU2006Data(&iau_data);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }

    size_t n_total_all   = 0;
    size_t n_written_all = 0;
    double et_first_all  = 0.0;
    double et_last_all   = 0.0;
    int    rc            = 0;

    for (int i = 0; i < n_inputs; ++i) {
        const char *input_rnx = input_rnx_paths[i];
        FILE *fin = fopen(input_rnx, "r");
        if (!fin) {
            fprintf(stderr, "gps: cannot open input '%s'\n", input_rnx);
            rc = 1;
            break;
        }
        size_t n_total_file = 0;
        int file_rc = _gps_scan_file(fin, fout, input_rnx, sat_id, &ctx,
                                      &n_total_file, &n_written_all,
                                      &et_first_all, &et_last_all);
        fclose(fin);
        n_total_all += n_total_file;
        if (file_rc != 0) { rc = file_rc; break; }
    }

    fclose(fout);

    if (rc == 0) {
        if (n_written_all == 0) {
            fprintf(stderr,
                "gps: WARNING -- no records written for sat_id '%s' "
                "across %d input file%s (scanned %zu RINEX nav messages "
                "total)\n",
                sat_id, n_inputs, n_inputs == 1 ? "" : "s", n_total_all);
        } else if (n_inputs > 1) {
            double duration_h = (et_last_all - et_first_all) / 3600.0;
            fprintf(stderr,
                "gps: aggregate -> %zu records across %d files "
                "(sat=%s, et=%.6f..%.6f, %.3f h, scanned %zu nav "
                "messages total)\n",
                n_written_all, n_inputs, sat_id,
                et_first_all, et_last_all, duration_h, n_total_all);
        }
    }

    /* Reference _parse_rinex_double here too so the file's first
     * unused-static warning is silenced -- _reparse_d_line uses an
     * inline equivalent and exposing the symbol keeps the helper
     * available for future records whose continuation lines have
     * D exponents we want to handle field-by-field. */
    (void)_parse_rinex_double;

    spody_free_MappedIAU2006(&iau_map);
    spody_free_MappedEOP(&eop_map);
    spody_free_MappedIAU2006Data(&iau_data);
    spody_free_MappedEOPData(&eop_data);
    return rc;
}
