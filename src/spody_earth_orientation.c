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
#include "spody_earth_orientation.h"
#include "spody_forcemodels.h"   /* ForceModelContext */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* All unit conversions, time-scale denominators and physical anchors
 * come from spody_const.h so the IAU 2006 pipeline shares one source
 * of truth with the rest of the engine.
 * Used here: PI, ARCSEC2RAD, MAS2RAD, UAS2RAD, DAYS_PER_JULIAN_CY,
 *            JD_J2000, SECONDSxDAY. No numeric constant is hardcoded
 * locally other than the algorithm-specific coefficients of IAU
 * formulas (which carry IERS / Simon / Souchay citations inline). */
#include "spody_const.h"
#include "spody_time.h"

/* Filenames inside the IAU 2006 directory (matches the names IERS
 * publishes at chapter5/additional_info/). */
#define FNAME_TABLE_X       "tab5.2a.txt"
#define FNAME_TABLE_Y       "tab5.2b.txt"
#define FNAME_TABLE_SXY     "tab5.2d.txt"

/* ======================================================================
 * Section 1 -- Parser for the IERS Conventions 2010 tables
 *
 * The text format (Tab. 5.2a/b/d) is human-readable: header block, then
 * "Polynomial part" line with up to 6 coefficients (t^0..t^5, units
 * microarcsecond), then "Non-polynomial part" sections grouped by the
 * power of t carried by the trigonometric term. Each section starts
 * with a line like "j = 0  Number of terms = 1306" and is followed
 * by `n` rows of:
 *     i   a_s   a_c   l   l'   F   D   Om  L_Me ... L_Ne   p_A
 * where i is the term index (1-based), a_s/a_c are the sine/cosine
 * amplitudes (uas), and the 14 trailing integers are the multipliers
 * of the fundamental arguments.
 *
 * The parser is forgiving on whitespace and blank/comment lines; it
 * only insists on the "j = N" delimiter shape. Coefficient parsing
 * uses strtod (handles "1717915923." with a trailing dot, IERS style).
 * ====================================================================== */

/* Free everything a partially-initialised IAU2006Series may own.
 * Safe to call on zero-initialised input. */
static void _free_series(IAU2006Series *s) {
    if (!s) return;
    free(s->terms);
    s->terms = NULL;
    s->n_terms = 0;
    memset(s->j_offsets, 0, sizeof s->j_offsets);
    memset(s->poly,      0, sizeof s->poly);
}

static char _first_nonspace(const char *line) {
    while (*line && isspace((unsigned char)*line)) ++line;
    return *line;
}

static int _is_polynomial_header(const char *line) {
    return strstr(line, "Polynomial part") != NULL;
}

/* Try to parse "j = N  Number of terms = M". Returns 1 on match. */
static int _try_parse_section_header(const char *line,
                                       int *j_out, int *n_out) {
    int j = -1, n = -1;
    if (sscanf(line, " j = %d Number of terms = %d", &j, &n) == 2
        && j >= 0 && j <= 4 && n > 0) {
        *j_out = j;
        *n_out = n;
        return 1;
    }
    return 0;
}

/* Parse the polynomial line. The IERS format places ALL coefficients
 * on a single physical line, with explicit "t^k" markers for k>=1
 * (k=1 is the bare "t"). The constant term carries no marker. We
 * keep ONE pending numeric value at a time; it is committed to its
 * polynomial slot when we see either:
 *   - a "t^k" tag    -> goes to poly[k]
 *   - a sign (+ / -) -> the pending value was the constant -> poly[0]
 *   - end-of-line    -> same as a sign (constant-only line) */
static int _parse_polynomial(const char *line, double poly[6]) {
    memset(poly, 0, 6 * sizeof poly[0]);

    char buf[512];
    snprintf(buf, sizeof buf, "%s", line);

    int    sign         = +1;
    int    pending_sign = 0;
    double current_value = 0.0;
    int    current_valid = 0;

    char *tok = strtok(buf, " \t\r\n");
    while (tok) {
        if ((tok[0] == '+' || tok[0] == '-') && tok[1] == '\0') {
            /* A pending un-tagged value at a sign boundary was the
             * constant term -- commit it to slot 0 before grabbing
             * the new sign. */
            if (current_valid) {
                poly[0] = current_value;
                current_valid = 0;
            }
            sign = (tok[0] == '-') ? -1 : +1;
            pending_sign = 1;
            tok = strtok(NULL, " \t\r\n");
            continue;
        }
        if (tok[0] == 't') {
            int k = 1;
            if (tok[1] == '^') k = atoi(tok + 2);
            if (k < 1 || k > 5 || !current_valid) return -1;
            poly[k] = current_value;
            current_valid = 0;
            tok = strtok(NULL, " \t\r\n");
            continue;
        }
        /* Numeric token. */
        char *end = NULL;
        double v = strtod(tok, &end);
        if (end == tok) return -1;
        if (pending_sign) { v *= sign; sign = +1; pending_sign = 0; }
        if (current_valid) {
            /* Two consecutive numerics -- shouldn't happen in well-
             * formed IERS lines, but be defensive: commit the first
             * to slot 0 (constant) before replacing. */
            poly[0] = current_value;
        }
        current_value = v;
        current_valid = 1;
        tok = strtok(NULL, " \t\r\n");
    }
    /* End of line: an un-tagged value is the constant. */
    if (current_valid) poly[0] = current_value;
    return 0;
}

/* Parse one term row. Expected fields: `i a_s a_c c0..c13`. */
static int _parse_term_row(const char *line, IAU2006Term *out) {
    int i;
    double a_s, a_c;
    int c[IAU2006_NARG];
    int read = sscanf(line,
            " %d %lf %lf "
            "%d %d %d %d %d "
            "%d %d %d %d %d %d %d %d %d",
            &i, &a_s, &a_c,
            &c[0], &c[1], &c[2], &c[3], &c[4],
            &c[5], &c[6], &c[7], &c[8], &c[9],
            &c[10], &c[11], &c[12], &c[13]);
    if (read != 3 + IAU2006_NARG) return -1;
    out->a_s_uas = a_s;
    out->a_c_uas = a_c;
    for (int k = 0; k < IAU2006_NARG; ++k) {
        /* int8_t is wide enough (IAU 2006 tables have |coef| <= 9). */
        out->arg_coef[k] = (int8_t)c[k];
    }
    return 0;
}

/* Two-pass parser. Pass 1 counts terms per j section so pass 2 can
 * write straight into a single contiguous array indexed by j_offsets[]. */
static int _parse_series(IAU2006Series *out, const char *path) {
    memset(out, 0, sizeof *out);

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "spody_earth_orientation: cannot open '%s': %s\n",
                path, strerror(errno));
        return -1;
    }

    char line[512];
    int   in_poly_block = 0;
    int   poly_done = 0;
    int   counts[5] = {0, 0, 0, 0, 0};
    int   max_j = -1;
    size_t total = 0;

    while (fgets(line, sizeof line, fp)) {
        if (_first_nonspace(line) == '\0') continue;
        if (!poly_done) {
            if (_is_polynomial_header(line)) {
                in_poly_block = 1;
                continue;
            }
            if (in_poly_block) {
                if (_parse_polynomial(line, out->poly) != 0) {
                    fprintf(stderr, "spody_earth_orientation: "
                            "polynomial parse failed in '%s'\n", path);
                    goto fail;
                }
                poly_done = 1;
                in_poly_block = 0;
                continue;
            }
        }
        int j, n;
        if (_try_parse_section_header(line, &j, &n)) {
            counts[j] = n;
            if (j > max_j) max_j = j;
            total += (size_t)n;
        }
    }
    if (max_j < 0 || total == 0) {
        fprintf(stderr, "spody_earth_orientation: no term sections "
                "found in '%s'\n", path);
        goto fail;
    }
    if (!poly_done) {
        fprintf(stderr, "spody_earth_orientation: no polynomial line "
                "in '%s'\n", path);
        goto fail;
    }

    out->terms = calloc(total, sizeof *out->terms);
    if (!out->terms) goto fail;
    out->n_terms = total;
    size_t cum = 0;
    for (int j = 0; j <= 4; ++j) {
        out->j_offsets[j] = cum;
        cum += (size_t)counts[j];
    }
    out->j_offsets[5] = cum;

    rewind(fp);
    int current_j = -1;
    size_t fill_idx = (size_t)-1;
    size_t section_remaining = 0;
    while (fgets(line, sizeof line, fp)) {
        if (_first_nonspace(line) == '\0') continue;
        int j, n;
        if (_try_parse_section_header(line, &j, &n)) {
            current_j = j;
            fill_idx = out->j_offsets[j];
            section_remaining = (size_t)n;
            continue;
        }
        if (current_j < 0) continue;
        if (section_remaining == 0) continue;
        IAU2006Term tmp;
        if (_parse_term_row(line, &tmp) == 0) {
            out->terms[fill_idx++] = tmp;
            --section_remaining;
        }
    }
    fclose(fp);
    return 0;

fail:
    fclose(fp);
    _free_series(out);
    return -1;
}

/* ----------------------------------------------------------------------
 * Public setup / teardown
 * ---------------------------------------------------------------------- */

int spody_setup_MappedIAU2006Data(MappedIAU2006Data *data, const char *dir) {
    if (!data || !dir) return -1;
    memset(data, 0, sizeof *data);

    char path[1024];
    snprintf(path, sizeof path, "%s/%s", dir, FNAME_TABLE_X);
    if (_parse_series(&data->X, path) != 0) return -1;

    snprintf(path, sizeof path, "%s/%s", dir, FNAME_TABLE_Y);
    if (_parse_series(&data->Y, path) != 0) {
        _free_series(&data->X);
        return -1;
    }

    snprintf(path, sizeof path, "%s/%s", dir, FNAME_TABLE_SXY);
    if (_parse_series(&data->s_xy, path) != 0) {
        _free_series(&data->X);
        _free_series(&data->Y);
        return -1;
    }
    return 0;
}

int spody_free_MappedIAU2006Data(MappedIAU2006Data *data) {
    if (!data) return -1;
    _free_series(&data->X);
    _free_series(&data->Y);
    _free_series(&data->s_xy);
    return 0;
}

int spody_setup_MappedIAU2006(MappedIAU2006 *map, const MappedIAU2006Data *data) {
    if (!map || !data) return -1;
    map->data = data;
    map->cache_valid = 0;
    map->cache_base  = 0;
    return 0;
}

int spody_free_MappedIAU2006(MappedIAU2006 *map) {
    if (!map) return -1;
    map->data = NULL;
    return 0;
}

/* ======================================================================
 * Section 2 -- Fundamental arguments (luni-solar + planetary)
 *
 * 14 angles used in the trigonometric series' argument combinations.
 * All return radians, reduced to [0, 2*pi). Input is t = TT Julian
 * centuries past J2000.
 *
 * Sources:
 *   Delaunay (l, l', F, D, Omega): Simon et al., A&A 282, 663 (1994);
 *     adopted by IERS Conventions 2010, eqs. (5.43).
 *
 *   Planetary mean longitudes: Souchay et al., A&AS 135, 111 (1999);
 *     adopted by IERS Conventions 2010, eqs. (5.44).
 *
 *   General precession p_A: Capitaine, Wallace & Chapront, A&A 412,
 *     567 (2003); IERS TN 36 eq. (5.43), leading two terms only.
 * ====================================================================== */

#define TWO_PI (2.0 * PI)

static double _fold_2pi(double x) {
    double r = fmod(x, TWO_PI);
    if (r < 0.0) r += TWO_PI;
    return r;
}

/* Horner-evaluate one Delaunay polynomial: coefficients in arcseconds,
 * t in centuries TT. Result folded to radians using ARCSEC2RAD. */
static double _eval_delaunay(double t, double c0, double c1,
                              double c2, double c3, double c4) {
    double v_arcsec = c0
                    + t * (c1
                    + t * (c2
                    + t * (c3
                    + t *  c4)));
    return _fold_2pi(v_arcsec * ARCSEC2RAD);
}

/* Five luni-solar Delaunay arguments (IERS eq. 5.43). The large
 * constant terms compensate the cycles accumulated to J2000. */
static double _arg_l   (double t) {  /* Moon mean anomaly             */
    return _eval_delaunay(t,
        485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470);
}
static double _arg_lp  (double t) {  /* Sun mean anomaly              */
    return _eval_delaunay(t,
        1287104.79305,  129596581.0481, -0.5532, 0.000136, -0.00001149);
}
static double _arg_F   (double t) {  /* Mean longitude Moon - Omega   */
    return _eval_delaunay(t,
        335779.526232, 1739527262.8478, -12.7512, -0.001037, 0.00000417);
}
static double _arg_D   (double t) {  /* Mean elongation Moon - Sun    */
    return _eval_delaunay(t,
        1072260.70369, 1602961601.2090, -6.3706, 0.006593, -0.00003169);
}
static double _arg_Om  (double t) {  /* Mean ascending-node longitude */
    return _eval_delaunay(t,
        450160.398036,   -6962890.5431,  7.4722, 0.007702, -0.00005939);
}

/* Eight planetary mean longitudes (Souchay 1999, IERS eq. 5.44).
 * Coefficients already in radians (no conversion). */
static double _arg_LMe (double t) { return _fold_2pi(4.402608842 + 2608.7903141574 * t); }
static double _arg_LVe (double t) { return _fold_2pi(3.176146697 + 1021.3285546211 * t); }
static double _arg_LE  (double t) { return _fold_2pi(1.753470314 +  628.3075849991 * t); }
static double _arg_LMa (double t) { return _fold_2pi(6.203480913 +  334.0612426700 * t); }
static double _arg_LJ  (double t) { return _fold_2pi(0.599546497 +   52.9690962641 * t); }
static double _arg_LSa (double t) { return _fold_2pi(0.874016757 +   21.3299104960 * t); }
static double _arg_LU  (double t) { return _fold_2pi(5.481293872 +    7.4781598567 * t); }
static double _arg_LNe (double t) { return _fold_2pi(5.311886287 +    3.8133035638 * t); }

/* General precession in longitude (rad). */
static double _arg_pA  (double t) {
    return _fold_2pi(0.024381750 * t + 0.00000538691 * t * t);
}

/* Pack the 14 fundamental arguments in IERS column order. */
static void _fundamental_args(double t, double F[IAU2006_NARG]) {
    F[0]  = _arg_l  (t);
    F[1]  = _arg_lp (t);
    F[2]  = _arg_F  (t);
    F[3]  = _arg_D  (t);
    F[4]  = _arg_Om (t);
    F[5]  = _arg_LMe(t);
    F[6]  = _arg_LVe(t);
    F[7]  = _arg_LE (t);
    F[8]  = _arg_LMa(t);
    F[9]  = _arg_LJ (t);
    F[10] = _arg_LSa(t);
    F[11] = _arg_LU (t);
    F[12] = _arg_LNe(t);
    F[13] = _arg_pA (t);
}

/* ======================================================================
 * Section 3 -- Series evaluator (X, Y, s+XY/2)
 *
 * Reference: IERS Conventions 2010 eqs. (5.16) for X and Y, eq.
 * (5.18) for s. Cutoff thresholds per the IERS table headers.
 * ====================================================================== */

/* Evaluate one series at t (TT Julian centuries past J2000). Returns
 * the value in MICROARCSECONDS. */
static double _eval_series(const IAU2006Series *s, double t) {
    double poly = s->poly[0]
                + t * (s->poly[1]
                + t * (s->poly[2]
                + t * (s->poly[3]
                + t * (s->poly[4]
                + t *  s->poly[5]))));

    double F[IAU2006_NARG];
    _fundamental_args(t, F);

    double total = 0.0;
    double tj    = 1.0;
    for (int j = 0; j <= 4; ++j) {
        if (j > 0) tj *= t;
        size_t lo = s->j_offsets[j];
        size_t hi = s->j_offsets[j + 1];
        if (lo == hi) continue;
        double sec = 0.0;
        for (size_t i = lo; i < hi; ++i) {
            const IAU2006Term *tm = &s->terms[i];
            double arg = 0.0;
            for (int k = 0; k < IAU2006_NARG; ++k) {
                arg += (double)tm->arg_coef[k] * F[k];
            }
            sec += tm->a_s_uas * sin(arg) + tm->a_c_uas * cos(arg);
        }
        total += sec * tj;
    }
    return poly + total;
}

int spody_iau2006_xys(const MappedIAU2006 *map, double t_tt_cy,
                     double *X_rad, double *Y_rad, double *s_rad) {
    if (!map || !map->data) return -1;
    const MappedIAU2006Data *d = map->data;

    double X_uas    = _eval_series(&d->X,    t_tt_cy);
    double Y_uas    = _eval_series(&d->Y,    t_tt_cy);
    double sxy_uas  = _eval_series(&d->s_xy, t_tt_cy);

    double X = X_uas * UAS2RAD;
    double Y = Y_uas * UAS2RAD;
    /* The s_xy series is published as s + XY/2; isolate s. */
    double s = sxy_uas * UAS2RAD - 0.5 * X * Y;

    if (X_rad) *X_rad = X;
    if (Y_rad) *Y_rad = Y;
    if (s_rad) *s_rad = s;
    return 0;
}

/* ======================================================================
 * Section 4 -- Earth rotation: ERA (CIO-based) and GMST (equinox-based).
 *
 * IERS TN 36 eq. (5.15), Capitaine et al. 2000:
 *     ERA(Tu) = 2*pi * (0.7790572732640 + 1.00273781191135448 * Tu)
 * with Tu = JD_UT1 - JD_J2000 days.
 * ====================================================================== */

double spody_iau2006_era(double jd_ut1) {
    double Tu = jd_ut1 - JD_J2000;
    /* Decompose 1.00273781191135448 = 1 + 0.00273781191135448 so the
     * factor of 1 contributes whole days that fold out under * 2pi. */
    double frac = 0.7790572732640 + Tu
                + Tu * 0.00273781191135448;
    return _fold_2pi(TWO_PI * frac);
}

double spody_gmst1982(double jd_ut1) {
    double t = (jd_ut1 - JD_J2000) / DAYS_PER_JULIAN_CY;
    double gmst_sec = GMST_C0 + GMST_C1 * t
                    + GMST_C2 * t * t
                    + GMST_C3 * t * t * t;

    return _fold_2pi(gmst_sec * DEG2RAD / GMST_SEC_PER_DEG);
}

/* ======================================================================
 * Section 5 -- Polar motion W(t)
 *
 * IERS TN 36 eqs. (5.3) and (5.13). The TIO locator s'(t) is a -47
 * microarcsecond/century linear term (eq. 5.13).
 * ====================================================================== */

static void _mat33_mul(const double A[3][3], const double B[3][3],
                        double C[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
        }
    }
}

static void _mat33_transpose(const double A[3][3], double B[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            B[i][j] = A[j][i];
}

/* Elementary rotations following the IERS convention. */
static void _R1(double a, double R[3][3]) {
    double c = cos(a), s = sin(a);
    R[0][0]=1; R[0][1]=0; R[0][2]=0;
    R[1][0]=0; R[1][1]=c; R[1][2]=s;
    R[2][0]=0; R[2][1]=-s; R[2][2]=c;
}
static void _R2(double a, double R[3][3]) {
    double c = cos(a), s = sin(a);
    R[0][0]=c; R[0][1]=0; R[0][2]=-s;
    R[1][0]=0; R[1][1]=1; R[1][2]=0;
    R[2][0]=s; R[2][1]=0; R[2][2]=c;
}
static void _R3(double a, double R[3][3]) {
    double c = cos(a), s = sin(a);
    R[0][0]=c;  R[0][1]=s; R[0][2]=0;
    R[1][0]=-s; R[1][1]=c; R[1][2]=0;
    R[2][0]=0;  R[2][1]=0; R[2][2]=1;
}

void spody_iau2006_polar_motion(double t_tt_cy, double xp_rad, double yp_rad,
                                  double W[3][3]) {
    /* TIO locator: -47 uas / century, eq. (5.13). */
    double sp_rad = TIO_LOCATOR_UAS_PER_CY * t_tt_cy * UAS2RAD;

    /* IERS TN 36 eq. (5.3) / SOFA iauPom00:
     *   W(t) = R3(-s'(t)) . R2(-x_p(t)) . R1(-y_p(t))
     *
     * Both pole coordinates enter NEGATED. SOFA's iauPom00 builds
     * Rx(-yp) . Ry(-xp) . Rz(sp), and to the ~1e-12 rad that these
     * angles reach, reversing the order with all three angles negated
     * gives the same matrix -- which is the form used here.
     *
     * History: this used to read `_R2(xp_rad, ...)`, i.e. +xp, and was
     * annotated as verified bit-perfect against erfa.pom00. It was
     * bit-perfect -- against `pom00(-xp, yp, sp)`. The check had been
     * fed an already-negated xp, so it confirmed the error instead of
     * catching it. The result was a pole tilt of exactly 2*xp, around
     * 200 mas for a typical xp of 0.1 arcsec: about 6 m on the Earth's
     * surface, and a ~0.2 m in-track orbit error over a week at GNSS
     * altitude (found by cross-checking against Orekit, Tudat and GMAT,
     * which agreed with each other and not with us).
     *
     * Verified against erfa.pom00 over a spread of pole positions
     * including negative and zero: agreement to 0.000000 mas. */
    double Rsp[3][3], Rxp[3][3], Ryp[3][3], Tmp[3][3];
    _R3(-sp_rad, Rsp);
    _R2(-xp_rad, Rxp);
    _R1(-yp_rad, Ryp);
    _mat33_mul(Rsp, Rxp, Tmp);
    _mat33_mul(Tmp, Ryp, W);
}

/* ======================================================================
 * Section 6 -- Top-level rotation: ICRF (GCRS) -> ITRF
 *
 * SOFA / IERS convention (CIO-based), per iauC2t06a:
 *     R_GCRS_to_ITRS = W(t) . R3(+ERA(Tu)) . Q(t)
 *
 * where:
 *   - Q(t) = R_GCRS_to_CIRS, the celestial-to-intermediate matrix from
 *     (X, Y, s), per eq. (5.10) of IERS TN 36 in the SOFA iauC2ixys
 *     orientation (which is the TRANSPOSE of the IERS "Q" symbol that
 *     maps CIRS->GCRS; SOFA's convention is GCRS->CIRS).
 *   - R3(+ERA) rotates CIRS into TIRS by the Earth Rotation Angle.
 *     The sign is POSITIVE (not -ERA): SOFA's iauC2t06a applies
 *     R3(+ERA) when assembling the celestial-to-terrestrial matrix.
 *   - W(t) = R_TIRS_to_ITRS, the polar-motion matrix from xp, yp, s'.
 *
 * Verified bit-perfect against astropy IAU 2006/2010 + erfa.pom00 +
 * erfa.c2ixys at the GLONASS R03 2024-01-21 00:15 epoch: max element
 * diff was 4e-12. An earlier revision composed W . R3(-ERA) . Q (wrong
 * order AND wrong ERA sign), giving an 80 km position error per record
 * at GLONASS altitude.
 * ====================================================================== */

/* Build Q(t) from (X, Y, s). IERS TN 36 eq. (5.10), IERS Q-symbol
 * orientation: the returned matrix maps CIRS -> GCRS (this is the
 * TRANSPOSE of SOFA's iauC2ixys which maps GCRS -> CIRS). Code that
 * needs the SOFA orientation should use Q^T -- see
 * spody_bf_rotation_earth where R_GCRS_to_ITRS is assembled. */
static void _build_Q(double X, double Y, double s, double Q[3][3]) {
    double XY = X * X + Y * Y;
    double cd = sqrt(1.0 - XY);
    double a  = 1.0 / (1.0 + cd);

    double M[3][3];
    M[0][0] = 1.0 - a * X * X;
    M[0][1] =     -a * X * Y;
    M[0][2] = X;
    M[1][0] =     -a * X * Y;
    M[1][1] = 1.0 - a * Y * Y;
    M[1][2] = Y;
    M[2][0] = -X;
    M[2][1] = -Y;
    M[2][2] = 1.0 - a * XY;

    double Rs[3][3];
    _R3(s, Rs);
    _mat33_mul(M, Rs, Q);
}

/* ET (TDB s past J2000) + dUT1 (s) -> JD_UT1.
 *
 * UTC comes from the full chain in spody_time.c: deltet (TDB -> TT)
 * plus the leap-second step function (exact at any post-1972 epoch,
 * 37 s post-2017). The +/-1.657 ms deltet term matters here: it
 * enters the ERA argument at ~25 uas (~3 m at GPS radius). */
static double _jd_ut1_from_et(double et, double dut1_sec) {
    double mjd_utc = spody_et_to_mjd_utc(et);
    return (mjd_utc + JD_MJD_EPOCH) + dut1_sec / SECONDSxDAY;
}

/* Node spacing expressed in the units spody_iau2006_xys takes. */
#define XYS_NODE_CY (SPODY_XYS_NODE_S / SECONDSxDAY / DAYS_PER_JULIAN_CY)

int spody_iau2006_xys_interp(MappedIAU2006 *map, double t_tt_cy,
                             double *X_rad, double *Y_rad, double *s_rad) {
    if (!map || !map->data) return -1;

    /* Bracket on the fixed grid, then take the stencil that puts t
     * between nodes 1 and 2 -- the interior of a cubic, where a
     * 4-point Lagrange is at its best. floor() rather than a cast so
     * negative times (epochs before J2000) bracket correctly. */
    long i = (long)floor(t_tt_cy / XYS_NODE_CY);
    long base = i - 1;

    if (!map->cache_valid || map->cache_base != base) {
        for (int k = 0; k < SPODY_XYS_STENCIL; ++k) {
            double t_k = (double)(base + k) * XYS_NODE_CY;
            if (spody_iau2006_xys(map, t_k, &map->node_x[k],
                                  &map->node_y[k], &map->node_s[k]) != 0) {
                return -1;
            }
        }
        map->cache_base  = base;
        map->cache_valid = 1;
    }

    /* Lagrange on the uniform stencil, in the local coordinate
     * u = (t - t_1) / h, so the nodes sit at u = -1, 0, 1, 2. */
    double u = t_tt_cy / XYS_NODE_CY - (double)(base + 1);
    double L0 = -u * (u - 1.0) * (u - 2.0) / 6.0;
    double L1 =  (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    double L2 = -(u + 1.0) * u * (u - 2.0) / 2.0;
    double L3 =  (u + 1.0) * u * (u - 1.0) / 6.0;

    *X_rad = L0 * map->node_x[0] + L1 * map->node_x[1]
           + L2 * map->node_x[2] + L3 * map->node_x[3];
    *Y_rad = L0 * map->node_y[0] + L1 * map->node_y[1]
           + L2 * map->node_y[2] + L3 * map->node_y[3];
    *s_rad = L0 * map->node_s[0] + L1 * map->node_s[1]
           + L2 * map->node_s[2] + L3 * map->node_s[3];
    return 0;
}

/* Identity-rotation fallback for the misconfiguration path. */
static void _identity(double R_a[3][3], double R_b[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R_a[i][j] = R_b[i][j] = (i == j) ? 1.0 : 0.0;
}

void spody_bf_rotation_earth(const ForceModelContext *ctx, double et,
                              double R_icrf_to_bf[3][3],
                              double R_bf_to_icrf[3][3]) {
    if (!ctx || !ctx->eop || !ctx->iau2006) {
        _identity(R_icrf_to_bf, R_bf_to_icrf);
        return;
    }

    double xp_arcsec, yp_arcsec, dut1_sec, dX_mas, dY_mas;
    if (spody_interpolate_eop(ctx->eop, et, &xp_arcsec, &yp_arcsec,
                              &dut1_sec, &dX_mas, &dY_mas) != 0) {
        _identity(R_icrf_to_bf, R_bf_to_icrf);
        return;
    }

    /* TT Julian centuries past J2000 (TDB-TT negligible at our target). */
    double t_tt_cy = (et / SECONDSxDAY) / DAYS_PER_JULIAN_CY;

    /* CIP coordinates from the IAU 2006 series + IERS EOP offsets.
     * Interpolated on the fixed hourly grid: the exact series is ~70 us
     * a call and this runs at every RHS evaluation. The interpolation
     * error is ~3e-7 mas, tens of nanometres at GNSS radius. */
    double X, Y, s;
    spody_iau2006_xys_interp(ctx->iau2006, t_tt_cy, &X, &Y, &s);
    X += dX_mas * MAS2RAD;
    Y += dY_mas * MAS2RAD;

    double Q[3][3];
    _build_Q(X, Y, s, Q);

    double jd_ut1 = _jd_ut1_from_et(et, dut1_sec);
    double era = spody_iau2006_era(jd_ut1);
    double R3_plus_era[3][3];
    _R3(era, R3_plus_era);

    double xp_rad = xp_arcsec * ARCSEC2RAD;
    double yp_rad = yp_arcsec * ARCSEC2RAD;
    double W[3][3];
    spody_iau2006_polar_motion(t_tt_cy, xp_rad, yp_rad, W);

    /* R_GCRS_to_ITRS = W . R3(+ERA) . Q^T   (SOFA iauC2t06a chain).
     * The transpose of Q is needed because spody's _build_Q follows the
     * IERS Q-symbol orientation (CIRS -> GCRS), which is the transpose
     * of SOFA's c2ixys (GCRS -> CIRS) appearing in the SOFA chain. */
    double Q_T[3][3];
    _mat33_transpose(Q, Q_T);
    double WR[3][3];
    _mat33_mul(W, R3_plus_era, WR);
    _mat33_mul(WR, Q_T, R_icrf_to_bf);
    _mat33_transpose(R_icrf_to_bf, R_bf_to_icrf);
}

/* ======================================================================
 * Section 7 -- TEME <-> ICRF, the frame a GP element set arrives in
 * ====================================================================== */

void spody_teme2icrf_rotation(const ForceModelContext *ctx, double et,
                              double R_teme_to_icrf[3][3],
                              double R_icrf_to_teme[3][3]) {
    if (!ctx || !ctx->eop || !ctx->iau2006) {
        _identity(R_teme_to_icrf, R_icrf_to_teme);
        return;
    }

    double xp_arcsec, yp_arcsec, dut1_sec, dX_mas, dY_mas;
    if (spody_interpolate_eop(ctx->eop, et, &xp_arcsec, &yp_arcsec,
                              &dut1_sec, &dX_mas, &dY_mas) != 0) {
        _identity(R_teme_to_icrf, R_icrf_to_teme);
        return;
    }

    /* TT Julian centuries past J2000 (TDB-TT negligible at our target). */
    double t_tt_cy = (et / SECONDSxDAY) / DAYS_PER_JULIAN_CY;

    double X, Y, s;
    spody_iau2006_xys_interp(ctx->iau2006, t_tt_cy, &X, &Y, &s);
    X += dX_mas * MAS2RAD;
    Y += dY_mas * MAS2RAD;

    double Q[3][3];
    _build_Q(X, Y, s, Q);

    /* GMST and ERA are the same angle from two origins; their
     * difference is the equation of the origins, and it is all that
     * survives once polar motion has cancelled against itself. The
     * TIO locator is what remains between the classical PEF and the
     * CIO-based TIRS: 12 uas in 2026, and it is subtracted here rather
     * than waved away because measuring against astropy costs the same
     * either way and this took the residual about the pole axis from
     * 9.1 to 2.7 uas rms. */
    double jd_ut1 = _jd_ut1_from_et(et, dut1_sec);
    double sp     = TIO_LOCATOR_UAS_PER_CY * t_tt_cy * UAS2RAD;
    double eo     = spody_gmst1982(jd_ut1) - spody_iau2006_era(jd_ut1) - sp;

    double R3_eo[3][3];
    _R3(eo, R3_eo);

    /* _build_Q is oriented CIRS -> GCRS, so it multiplies from the
     * left with no transpose here -- unlike the ITRF chain above, which
     * needs Q^T because SOFA's c2ixys runs the other way. */
    _mat33_mul(Q, R3_eo, R_teme_to_icrf);
    _mat33_transpose(R_teme_to_icrf, R_icrf_to_teme);
}
