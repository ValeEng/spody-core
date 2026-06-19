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
 * Earth Orientation Parameters (EOP) loader and query.
 *
 * Reads the IERS Bulletin A/B combined file `finals2000A.all` and
 * exposes interpolated EOP values to the body-fixed rotation code
 * (`spody_bf_rotation_earth`, P2.2). The loader follows the same
 * "parse text once, hold in heap, query thread-safely" pattern as
 * HarmonicGravityData: no memory-mapped binary cache, no per-thread
 * mutation -- a single MappedEOPData is shared by every worker.
 *
 * SOURCE
 * ------
 *   IERS finals2000A.all
 *     official: https://datacenter.iers.org/data/9/finals2000A.all
 *     mirror  : https://maia.usno.navy.mil/ser7/finals2000A.all
 *   Format spec:
 *     https://maia.usno.navy.mil/ser7/readme.finals2000A
 *   The file is a daily-stepped text table mixing IERS Bulletin A
 *   (Rapid Service, daily updates with ~365 days of predictions)
 *   and Bulletin B (Final values, ~1 month latency). Bulletin B
 *   columns supersede Bulletin A on the days where both are
 *   present; predictions use only Bulletin A.
 *
 * MODEL CONVENTION
 * ----------------
 *   The file ships the corrections for the IAU 2000A precession-
 *   nutation theory (dX, dY -- celestial pole offsets in CIP
 *   coordinates). They feed directly into IAU 2006/2010 ICRF<->ITRF
 *   rotation per IERS TN 36, Chapter 5. The legacy IAU 1980 dPsi,
 *   dEps columns (older finals.all naming) are NOT consumed.
 *
 *   Polar motion (xp, yp) are in arcseconds, UT1-UTC in seconds,
 *   dX/dY in milliarcseconds. Units match the file as-is; callers
 *   convert when feeding the rotation pipeline.
 */
#ifndef SPODY_EOP_H
#define SPODY_EOP_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* One parsed EOP record. Mirrors the columns of finals2000A.all that
 * we actually need; the file ships extra columns (errors, LOD) that
 * the rotation model does not consume. Values that are blank in the
 * source (predictions for Bulletin B, or LOD not yet computed) are
 * stored as 0 with the `has_bulletin_b` flag cleared so the caller
 * can distinguish observed from predicted at a glance. */
typedef struct {
    double mjd;             /* Modified Julian Date (UTC), 7 decimals    */
    double xp_arcsec;       /* polar motion x, " (Bulletin B if avail)   */
    double yp_arcsec;       /* polar motion y, " (Bulletin B if avail)   */
    double dut1_sec;        /* UT1-UTC, s (Bulletin B if avail)          */
    double dx_mas;          /* CIP correction dX, mas (Bulletin B/A)     */
    double dy_mas;          /* CIP correction dY, mas (Bulletin B/A)     */
    int    has_bulletin_b;  /* 1 = Bulletin B (observed final), 0 = pred */
} EOPRecord;

/* Shared, read-only EOP table. `records` is heap-allocated, sorted
 * by ascending MJD, contiguous (no gaps -- the file is daily-stepped
 * with no missing rows). Thread-safe to share across workers because
 * no mutating access is exposed: every query is a binary search +
 * linear interpolation that touches the table read-only.
 *
 * `mjd_last_observed` is the MJD of the last record that carries
 * Bulletin B final values; everything after is IERS prediction. The
 * staleness gate (`spody_eop_last_observed_mjd`) reports this so
 * the app can warn the user when the run window crosses into the
 * predicted region. `mjd_last_predicted` is the hard upper bound:
 * queries beyond it return an error. */
typedef struct {
    EOPRecord *records;
    size_t     n_records;
    double     mjd_first;
    double     mjd_last_observed;
    double     mjd_last_predicted;
} MappedEOPData;

/* Per-thread query handle. Holds the shared MappedEOPData pointer
 * plus a tiny last-hit cache (record-index) so consecutive queries
 * at nearby ETs (typical integrator pattern) skip the bsearch. */
typedef struct {
    const MappedEOPData *med;
    size_t cached_idx;          /* last bracketing-record-low index */
    int    cached_valid;
} MappedEOP;

/* Parse `filename` (an IERS finals2000A.all text file) into *med.
 * Heap-allocates med->records; the caller must spody_free_MappedEOPData
 * on shutdown. Returns 0 on success, -1 on parse / IO failure. */
int spody_setup_MappedEOPData(MappedEOPData *med, const char *filename);

/* Release med->records. Safe on a zero-initialised struct. */
int spody_free_MappedEOPData(MappedEOPData *med);

/* Bind a per-thread MappedEOP handle to a shared MappedEOPData.
 * Does not allocate; just stores the pointer and resets the cache. */
int spody_setup_MappedEOP(MappedEOP *map, const MappedEOPData *med);

/* Release the per-thread handle. Today a no-op; provided for
 * symmetry with the other Mapped* objects. */
int spody_free_MappedEOP(MappedEOP *map);

/* Interpolate EOP values at Ephemeris Time `et` (seconds past J2000
 * TDB). Linear interpolation between the two bracketing daily
 * records is sufficient for the IAU 2006/2010 model's mas-level
 * target -- the EOP series is itself daily-sampled and IERS does
 * not publish sub-day variations.
 *
 * Internally converts ET -> UTC MJD (the EOP table's time scale)
 * using the standard TT-TAI = +32.184s and TAI-UTC leap-second
 * table; the leap-second corrections live in `spody_const.h`.
 *
 * Returns 0 on success, -1 when `et` falls outside the table's
 * coverage [mjd_first, mjd_last_predicted] (caller decides whether
 * to abort or extrapolate). Out parameters are written only on
 * success; on failure they are untouched. Any of them may be NULL
 * if the caller doesn't need that quantity. */
int spody_interpolate_eop(MappedEOP *map, double et,
                          double *xp_arcsec, double *yp_arcsec,
                          double *dut1_sec,
                          double *dx_mas,    double *dy_mas);

/* MJD of the last record carrying Bulletin B (observed) data. Past
 * this MJD the table is IERS prediction; precision degrades from
 * mas (observed) to ~10 mas at +1 month, ~100 mas at +6 months,
 * then several arcseconds towards the +365 day prediction horizon.
 *
 * The app uses this to warn the user when the run window crosses
 * into the predicted region -- still usable, but the user should
 * know. */
double spody_eop_last_observed_mjd(const MappedEOPData *med);

/* MJD of the last record in the table (hard upper bound, including
 * IERS predictions). Queries beyond this MJD fail. The app uses
 * this to refuse a run that would integrate past the table's
 * predicted horizon, prompting the user to update the file. */
double spody_eop_last_predicted_mjd(const MappedEOPData *med);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_EOP_H */
