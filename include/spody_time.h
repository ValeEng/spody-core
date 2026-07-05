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
#ifndef SPODY_TIME_H
#define SPODY_TIME_H

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * Shared calendar / time-scale helpers.
 *
 * Single home for the conversions that used to be copy-pasted across
 * the GNSS converters (gps, glonass, sp3) and the UTC consumers
 * (eop, atmosphere, earth_orientation):
 *
 *   - Gregorian civil date -> Julian Date (Meeus)
 *   - TAI-UTC leap-second step function (IERS Bulletin C chain)
 *   - TDB-TT periodic term (SPICE deltet)
 *   - ET (TDB s past J2000) -> UTC MJD
 *
 * The leap-second table lives in spody_time.c and NOWHERE else on the
 * C side. When IERS announces a new leap second, update that table
 * (one row) and the mirror in python/spopy/time.py.
 * ---------------------------------------------------------------------- */

/* Gregorian civil date -> Julian Date (JD at midnight + day fraction).
 * Standard Meeus algorithm (Astronomical Algorithms, ch. 7); valid for
 * any Gregorian-calendar date (post 1582-10-15). The y/m/d/hh/mn/ss
 * fields are wall-clock values in whatever time scale the caller is
 * working in (GPST, UTC, ...); the returned JD is in that same scale. */
double spody_greg_to_jd(int y, int m, int d, int hh, int mn, double ss);

/* TAI - UTC in seconds at the given UTC MJD. Piecewise-constant step
 * function over the post-1972 IERS Bulletin C insertions (10 s at
 * 1972-01-01 up to 37 s since 2017-01-01). For MJDs before the first
 * table entry returns 10.0 -- the pre-1972 fractional TAI-UTC history
 * is irrelevant for GNSS / Earth-orbit propagation. */
double spody_tai_minus_utc(double mjd_utc);

/* TDB - TT in seconds at ET (TDB seconds past J2000): the SPICE
 * `deltet` periodic term, K*sin(E) with a one-step Kepler solve
 * (constants DELTET_* in spody_const.h, from the NAIF LSK kernel).
 * Amplitude +/-1.657 ms. The argument may be ET or TT seconds
 * interchangeably: swapping them changes the result by less than a
 * picosecond, far below the formula's own fidelity, so both the
 * ET->UTC and the UTC->ET direction call this without iterating. */
double spody_tdb_minus_tt(double et);

/* ET (TDB seconds past J2000) -> UTC MJD. Full chain:
 * TT = ET - deltet, TAI = TT - 32.184 s, UTC = TAI - leap(UTC).
 * TAI-UTC is a step function of UTC itself, so a two-iteration fixed
 * point evaluates it first at TAI, then at the estimated UTC (the
 * second pass only matters within 37 s of a leap boundary). Matches
 * the SPICE ET->UTC conversion to sub-us. */
double spody_et_to_mjd_utc(double et);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_TIME_H */
