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
#include <math.h>
#include <stddef.h>
#include "spody_time.h"
#include "spody_const.h"

/* ----------------------------------------------------------------------
 * Leap-second table: (MJD_UTC at which the leap second is INTRODUCED,
 * TAI-UTC in effect from that MJD on). Post-1972 IERS Bulletin C
 * chain; 37 s since 2017-01-01, none scheduled since (BIPM CGPM 2022
 * recommends phasing leap seconds out by/around 2035).
 *
 * THE single C-side copy. Python mirror: python/spopy/time.py.
 * ---------------------------------------------------------------------- */
typedef struct { double mjd_utc; double tai_minus_utc; } LeapEntry;

static const LeapEntry leap_table[] = {
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
    { 57754.0, 37.0 },  /* 2017-01-01 -- current value */
};
static const size_t leap_table_n = sizeof leap_table / sizeof leap_table[0];

double spody_greg_to_jd(int y, int m, int d, int hh, int mn, double ss) {
    if (m <= 2) { y -= 1; m += 12; }
    int A = y / 100;
    int B = 2 - A + (A / 4);
    double jd_midnight = floor(365.25 * (double)(y + 4716))
                        + floor(30.6001 * (double)(m + 1))
                        + (double)d + (double)B - 1524.5;
    double day_frac = ((double)hh * 3600.0 + (double)mn * 60.0 + ss)
                      / SECONDSxDAY;
    return jd_midnight + day_frac;
}

void spody_mjd_to_doy(double mjd, int *year_out, int *doy_out,
                      double *sec_of_day_out) {
    /* Day / time-of-day split first: everything below is integer. */
    double day = floor(mjd);
    if (sec_of_day_out) *sec_of_day_out = (mjd - day) * SECONDSxDAY;

    /* Fliegel-Van Flandern inverse (JDN -> Gregorian civil date).
     * JDN is noon-based, so the civil date containing this MJD's
     * midnight has JDN = floor(mjd) + 2400001. All-integer
     * arithmetic; positive for any MJD >= 0 (post-1858), where C
     * truncation and Python floor division agree -- keep it that way
     * or the spopy twin breaks. */
    long l = (long)day + 2400001L + 68569L;
    long n = (4L * l) / 146097L;
    long i, j, d, m, y;
    l = l - (146097L * n + 3L) / 4L;
    i = (4000L * (l + 1L)) / 1461001L;
    l = l - (1461L * i) / 4L + 31L;
    j = (80L * l) / 2447L;
    d = l - (2447L * j) / 80L;
    l = j / 11L;
    m = j + 2L - 12L * l;
    y = 100L * (n - 49L) + i + l;

    if (year_out) *year_out = (int)y;
    if (doy_out) {
        static const int cum_days[12] =
            {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
        int leap = ((y % 4 == 0 && y % 100 != 0) || y % 400 == 0);
        *doy_out = cum_days[m - 1] + (int)d + ((m > 2 && leap) ? 1 : 0);
    }
}

double spody_tai_minus_utc(double mjd_utc) {
    if (mjd_utc < leap_table[0].mjd_utc) return leap_table[0].tai_minus_utc;
    for (size_t i = leap_table_n; i-- > 0; ) {
        if (mjd_utc >= leap_table[i].mjd_utc) return leap_table[i].tai_minus_utc;
    }
    return leap_table[0].tai_minus_utc;
}

double spody_tdb_minus_tt(double et) {
    double m = DELTET_M0 + DELTET_M1 * et;
    double e = m + DELTET_EB * sin(m);
    return DELTET_K * sin(e);
}

double spody_et_to_mjd_utc(double et) {
    /* TT2TAI_SEC = TAI - TT, so adding it walks TT -> TAI. The first
     * pass evaluates the step function at TAI (off by <=37 s from
     * UTC); the second pass re-evaluates at the estimated UTC, which
     * only matters within 37 s of a leap-second boundary. */
    double tt_sec  = et - spody_tdb_minus_tt(et);
    double mjd_tt  = (tt_sec / SECONDSxDAY) + (JD_J2000 - JD_MJD_EPOCH);
    double mjd_tai = mjd_tt + TT2TAI_SEC / SECONDSxDAY;
    double leap    = spody_tai_minus_utc(mjd_tai);
    double mjd_utc = mjd_tai - leap / SECONDSxDAY;
    leap    = spody_tai_minus_utc(mjd_utc);
    mjd_utc = mjd_tai - leap / SECONDSxDAY;
    return mjd_utc;
}
