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
 *   - ET (TDB s past J2000) -> UTC MJD
 *
 * The leap-second table lives in spody_time.c and NOWHERE else on the
 * C side. When IERS announces a new leap second, update that table
 * (one row) and the mirror in python/spopy/eop.py.
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

/* ET (TDB seconds past J2000) -> UTC MJD.
 * Two-iteration fixed point around the leap-second step function:
 * TAI-UTC is a function of UTC itself, so a first estimate treats
 * TT - 32.184 s as UTC, then the offset is re-evaluated at the
 * estimated UTC. TDB-TT (<2 ms) is neglected: far below the daily
 * sampling of every UTC-indexed table this feeds (EOP, space weather). */
double spody_et_to_mjd_utc(double et);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_TIME_H */
