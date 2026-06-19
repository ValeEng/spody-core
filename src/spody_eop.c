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
 * spody_eop.c -- IERS finals2000A.all parser + interpolator.
 *
 * The on-disk format is a daily-stepped text table; parsing is one
 * pass through the file, allocating one EOPRecord per non-empty line
 * that carries at least Bulletin A polar motion + UT1-UTC. Column
 * widths follow the readme.finals2000A spec from USNO:
 *
 *   https://maia.usno.navy.mil/ser7/readme.finals2000A
 *
 * We always prefer Bulletin B (final, IERS reanalysis at ~1 month
 * latency) when present and fall back to Bulletin A (rapid service +
 * predictions) otherwise. The has_bulletin_b flag is set accordingly.
 *
 * Time scale: the file's MJD column is UTC (per IERS spec). To use
 * it from ET (the SPODY-internal time), we run a textbook chain:
 *     ET  = seconds (TDB) past J2000      (input to interp)
 *     TT  ~= ET                            (TDB-TT < 2 ms, ignored)
 *     TAI = TT  - 32.184 s
 *     UTC = TAI - leap_seconds(UTC)        (integer-second jumps)
 * The leap-second table below covers 1972-01-01..present (37 leap
 * seconds added). Adding a new leap second is one row at the bottom.
 */
#include "spody_eop.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_const.h"

/* ----------------------------------------------------------------------
 * Time-scale plumbing: ET (TDB) -> UTC MJD
 * ----------------------------------------------------------------------
 *
 * Leap-second table: (MJD_UTC at which leap second is INTRODUCED,
 * cumulative TAI - UTC in seconds after that introduction).
 *
 * Source: IERS Bulletin C history,
 *   https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat
 *
 * The last leap second was inserted at the end of 2016-12-31; no
 * new one has been scheduled since (BIPM CGPM 2022 resolution
 * recommends phasing out leap seconds by/around 2035). Until then
 * this table must be updated when IERS Bulletin C announces a new
 * leap-second insertion.
 */
typedef struct { double mjd_utc; double tai_minus_utc; } LeapEntry;

static const LeapEntry _leap_table[] = {
    { 41317.0, 10.0 },  /* 1972-01-01 -- first IERS-defined entry  */
    { 41499.0, 11.0 },  /* 1972-07-01 */
    { 41683.0, 12.0 },  /* 1973-01-01 */
    { 42048.0, 13.0 },  /* 1974-01-01 */
    { 42413.0, 14.0 },  /* 1975-01-01 */
    { 42778.0, 15.0 },  /* 1976-01-01 */
    { 43144.0, 16.0 },  /* 1977-01-01 */
    { 43509.0, 17.0 },  /* 1978-01-01 */
    { 43874.0, 18.0 },  /* 1979-01-01 */
    { 44239.0, 19.0 },  /* 1980-01-01 */
    { 44786.0, 20.0 },  /* 1981-07-01 */
    { 45151.0, 21.0 },  /* 1982-07-01 */
    { 45516.0, 22.0 },  /* 1983-07-01 */
    { 46247.0, 23.0 },  /* 1985-07-01 */
    { 47161.0, 24.0 },  /* 1988-01-01 */
    { 47892.0, 25.0 },  /* 1990-01-01 */
    { 48257.0, 26.0 },  /* 1991-01-01 */
    { 48804.0, 27.0 },  /* 1992-07-01 */
    { 49169.0, 28.0 },  /* 1993-07-01 */
    { 49534.0, 29.0 },  /* 1994-07-01 */
    { 50083.0, 30.0 },  /* 1996-01-01 */
    { 50630.0, 31.0 },  /* 1997-07-01 */
    { 51179.0, 32.0 },  /* 1999-01-01 */
    { 53736.0, 33.0 },  /* 2006-01-01 */
    { 54832.0, 34.0 },  /* 2009-01-01 */
    { 56109.0, 35.0 },  /* 2012-07-01 */
    { 57204.0, 36.0 },  /* 2015-07-01 */
    { 57754.0, 37.0 },  /* 2017-01-01 -- current value as of 2024  */
};
static const size_t _leap_table_n =
        sizeof _leap_table / sizeof _leap_table[0];

/* Return TAI-UTC at the given UTC MJD. Step-function: piecewise
 * constant between consecutive leap-second insertions. */
static double _tai_minus_utc(double mjd_utc) {
    /* For MJDs before the first table entry we return 10 -- the
     * pre-1972 TAI-UTC offset history is fractional and irrelevant
     * for GNSS/Earth-orbit propagation. */
    if (mjd_utc < _leap_table[0].mjd_utc) return _leap_table[0].tai_minus_utc;
    for (size_t i = _leap_table_n; i-- > 0; ) {
        if (mjd_utc >= _leap_table[i].mjd_utc) {
            return _leap_table[i].tai_minus_utc;
        }
    }
    return _leap_table[0].tai_minus_utc;
}

/* ET (TDB s past J2000) -> UTC MJD.
 * Single-iteration fixed-point: TAI-UTC is a step function of UTC,
 * so we first estimate UTC from TT-32.184s (treating it as TAI), then
 * subtract the leap-second offset to get the actual UTC. The estimate
 * is wrong by at most ~37 seconds (current leap offset), which never
 * straddles a leap-second boundary unless the user picks an ET within
 * 37 seconds of midnight UTC on a leap day -- in which case one more
 * iteration would fix it. We do the second iteration to be safe.
 *
 * TDB-TT correction is <2 ms over the EOP table's century-scale
 * coverage and well below the daily EOP sampling, so we treat
 * TT == ET.*/
static double _et_to_mjd_utc(double et) {
    double mjd_tt  = (et / SECONDSxDAY) + (JD_J2000 - 2400000.5);
    double mjd_tai = mjd_tt - 32.184 / SECONDSxDAY;
    /* First pass: guess UTC = TAI - leap(at TAI). */
    double leap = _tai_minus_utc(mjd_tai);
    double mjd_utc = mjd_tai - leap / SECONDSxDAY;
    /* Second pass with the better UTC argument to leap(). */
    leap = _tai_minus_utc(mjd_utc);
    mjd_utc = mjd_tai - leap / SECONDSxDAY;
    return mjd_utc;
}

/* ----------------------------------------------------------------------
 * Text parser for finals2000A.all
 *
 * Spec: https://maia.usno.navy.mil/ser7/readme.finals2000A
 *
 * Column widths (1-based, inclusive):
 *   8-15    F8.2   MJD (UTC)
 *   18-27   F10.6  Bulletin A PM-x  (")
 *   38-46   F9.6   Bulletin A PM-y  (")
 *   59-68   F10.7  Bulletin A UT1-UTC (s)
 *   98-106  F9.3   Bulletin A dX wrt IAU2000A  (mas)
 *   118-126 F9.3   Bulletin A dY wrt IAU2000A  (mas)
 *   135-144 F10.6  Bulletin B PM-x  (")
 *   145-154 F10.6  Bulletin B PM-y  (")
 *   155-165 F11.7  Bulletin B UT1-UTC (s)
 *   166-175 F10.3  Bulletin B dX  (mas)
 *   176-185 F10.3  Bulletin B dY  (mas)
 *
 * Bulletin B columns are blank-filled when not yet available
 * (i.e. for prediction days). The parser detects blank by testing
 * isspace() on the column's first non-leading character.
 * ---------------------------------------------------------------------- */

/* Read a fixed-width F-style numeric field. `col` is 1-based per the
 * IERS spec; `width` is the length. Returns 1 on success (value
 * written), 0 when the field is blank (value left untouched). */
static int _read_f_field(const char *line, int col, int width,
                          double *out) {
    char buf[24];
    if (width <= 0 || width >= (int)sizeof buf) return 0;
    int line_len = (int)strlen(line);
    if (col - 1 + width > line_len) return 0;
    memcpy(buf, line + col - 1, (size_t)width);
    buf[width] = '\0';
    /* Blank field detection: all whitespace. */
    int all_blank = 1;
    for (int i = 0; i < width; ++i) {
        if (!isspace((unsigned char)buf[i])) { all_blank = 0; break; }
    }
    if (all_blank) return 0;
    char *end = NULL;
    double v = strtod(buf, &end);
    if (end == buf) return 0;
    *out = v;
    return 1;
}

int spody_setup_MappedEOPData(MappedEOPData *med, const char *filename) {
    if (!med || !filename) return -1;
    memset(med, 0, sizeof *med);

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "spody_eop: cannot open '%s': %s\n",
                filename, strerror(errno));
        return -1;
    }

    /* Two-pass: first count valid lines so we allocate exactly. The
     * file is ~3 MB / ~15k lines -- a second pass costs sub-ms. */
    char line[512];
    size_t n = 0;
    while (fgets(line, sizeof line, fp)) {
        double mjd, xa, ya, dut1a;
        if (_read_f_field(line, 8,  8, &mjd)   &&
            _read_f_field(line, 18, 10, &xa)   &&
            _read_f_field(line, 38, 9, &ya)    &&
            _read_f_field(line, 59, 10, &dut1a)) {
            ++n;
        }
    }
    if (n == 0) {
        fprintf(stderr, "spody_eop: no usable records in '%s'\n", filename);
        fclose(fp);
        return -1;
    }

    med->records = calloc(n, sizeof *med->records);
    if (!med->records) {
        fclose(fp);
        return -1;
    }
    med->n_records = n;

    rewind(fp);
    size_t k = 0;
    double mjd_last_b = -INFINITY;
    while (k < n && fgets(line, sizeof line, fp)) {
        EOPRecord rec;
        memset(&rec, 0, sizeof rec);
        double xa, ya, dut1a, dxa, dya;
        double xb, yb, dut1b, dxb, dyb;
        if (!_read_f_field(line, 8, 8, &rec.mjd))         continue;
        if (!_read_f_field(line, 18, 10, &xa))            continue;
        if (!_read_f_field(line, 38, 9, &ya))             continue;
        if (!_read_f_field(line, 59, 10, &dut1a))         continue;

        /* Bulletin A: dX/dY are mandatory for IAU 2006 modelling but
         * absent from very old records. We accept the row even if
         * the dX/dY columns are blank (set to 0, equivalent to
         * "no CIP correction"). */
        if (!_read_f_field(line, 98, 9, &dxa))   dxa = 0.0;
        if (!_read_f_field(line, 118, 9, &dya))  dya = 0.0;

        /* Bulletin B columns: blank on prediction rows. Use them when
         * present, otherwise keep the Bulletin A values. */
        int has_b = 1;
        has_b &= _read_f_field(line, 135, 10, &xb);
        has_b &= _read_f_field(line, 145, 10, &yb);
        has_b &= _read_f_field(line, 155, 11, &dut1b);
        /* dX/dY in Bulletin B may also be blank even when xp/yp/UT1
         * are present -- accept the row and fall back to Bulletin A. */
        int b_has_dxy = _read_f_field(line, 166, 10, &dxb) &&
                         _read_f_field(line, 176, 10, &dyb);

        if (has_b) {
            rec.xp_arcsec = xb;
            rec.yp_arcsec = yb;
            rec.dut1_sec  = dut1b;
            rec.dx_mas    = b_has_dxy ? dxb : dxa;
            rec.dy_mas    = b_has_dxy ? dyb : dya;
            rec.has_bulletin_b = 1;
            if (rec.mjd > mjd_last_b) mjd_last_b = rec.mjd;
        } else {
            rec.xp_arcsec = xa;
            rec.yp_arcsec = ya;
            rec.dut1_sec  = dut1a;
            rec.dx_mas    = dxa;
            rec.dy_mas    = dya;
            rec.has_bulletin_b = 0;
        }

        med->records[k++] = rec;
    }
    fclose(fp);
    /* Resize down if our two-pass count overcounted (defensive: the
     * second pass uses the same predicate, so this should never fire
     * unless the file was modified between passes). */
    if (k < med->n_records) med->n_records = k;
    if (k == 0) {
        free(med->records);
        med->records = NULL;
        return -1;
    }

    med->mjd_first           = med->records[0].mjd;
    med->mjd_last_predicted  = med->records[med->n_records - 1].mjd;
    med->mjd_last_observed   = isfinite(mjd_last_b) ? mjd_last_b
                                                     : med->mjd_first;
    return 0;
}

int spody_free_MappedEOPData(MappedEOPData *med) {
    if (!med) return -1;
    free(med->records);
    med->records = NULL;
    med->n_records = 0;
    return 0;
}

int spody_setup_MappedEOP(MappedEOP *map, const MappedEOPData *med) {
    if (!map || !med) return -1;
    map->med           = med;
    map->cached_idx    = 0;
    map->cached_valid  = 0;
    return 0;
}

int spody_free_MappedEOP(MappedEOP *map) {
    if (!map) return -1;
    map->med = NULL;
    map->cached_valid = 0;
    return 0;
}

/* Binary search for the largest index `i` such that records[i].mjd <= mjd.
 * Records are uniformly-stepped daily so a linear hop from the cached
 * index handles the common integrator case (consecutive queries within
 * a few days) without log-overhead. */
static size_t _bracket_idx(const MappedEOPData *med, size_t hint, double mjd) {
    if (med->n_records == 0) return 0;
    /* Fast path: hint already brackets. */
    if (hint < med->n_records - 1 &&
        med->records[hint].mjd <= mjd &&
        med->records[hint + 1].mjd > mjd) {
        return hint;
    }
    /* Linear hop forward (typical integrator step is <1 day). */
    if (hint < med->n_records && med->records[hint].mjd <= mjd) {
        size_t i = hint;
        while (i + 1 < med->n_records && med->records[i + 1].mjd <= mjd) ++i;
        return i;
    }
    /* Fall back to binary search. */
    size_t lo = 0, hi = med->n_records;
    while (lo + 1 < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (med->records[mid].mjd <= mjd) lo = mid; else hi = mid;
    }
    return lo;
}

int spody_interpolate_eop(MappedEOP *map, double et,
                          double *xp_arcsec, double *yp_arcsec,
                          double *dut1_sec,
                          double *dx_mas,    double *dy_mas) {
    if (!map || !map->med || map->med->n_records == 0) return -1;
    const MappedEOPData *med = map->med;

    double mjd = _et_to_mjd_utc(et);
    if (mjd < med->mjd_first || mjd > med->mjd_last_predicted) {
        return -1;   /* out of coverage */
    }

    size_t i = _bracket_idx(med, map->cached_valid ? map->cached_idx : 0, mjd);
    map->cached_idx = i;
    map->cached_valid = 1;

    /* Boundary case: mjd lands exactly on the last record. */
    if (i + 1 >= med->n_records) {
        const EOPRecord *r = &med->records[i];
        if (xp_arcsec) *xp_arcsec = r->xp_arcsec;
        if (yp_arcsec) *yp_arcsec = r->yp_arcsec;
        if (dut1_sec)  *dut1_sec  = r->dut1_sec;
        if (dx_mas)    *dx_mas    = r->dx_mas;
        if (dy_mas)    *dy_mas    = r->dy_mas;
        return 0;
    }

    const EOPRecord *lo = &med->records[i];
    const EOPRecord *hi = &med->records[i + 1];
    double dmjd = hi->mjd - lo->mjd;
    /* Records are nominally daily-stepped; guard against duplicate
     * MJDs (defensive, shouldn't happen with finals2000A.all). */
    double frac = (dmjd > 0.0) ? (mjd - lo->mjd) / dmjd : 0.0;

    if (xp_arcsec) *xp_arcsec = lo->xp_arcsec + frac * (hi->xp_arcsec - lo->xp_arcsec);
    if (yp_arcsec) *yp_arcsec = lo->yp_arcsec + frac * (hi->yp_arcsec - lo->yp_arcsec);
    if (dut1_sec)  *dut1_sec  = lo->dut1_sec  + frac * (hi->dut1_sec  - lo->dut1_sec);
    if (dx_mas)    *dx_mas    = lo->dx_mas    + frac * (hi->dx_mas    - lo->dx_mas);
    if (dy_mas)    *dy_mas    = lo->dy_mas    + frac * (hi->dy_mas    - lo->dy_mas);
    return 0;
}

double spody_eop_last_observed_mjd(const MappedEOPData *med) {
    return med ? med->mjd_last_observed : 0.0;
}

double spody_eop_last_predicted_mjd(const MappedEOPData *med) {
    return med ? med->mjd_last_predicted : 0.0;
}
