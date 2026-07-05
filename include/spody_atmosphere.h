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
 * Atmospheric density model interface + space weather data.
 *
 * Two things live in this module:
 *
 *   1. The contract every atmosphere model implements
 *      (`SpodyAtmosphereDensityFn` callback + `SpodyAtmosphere` handle).
 *      The drag force (`spody_force_drag`, spody_forcemodels.c) calls
 *      the callback
 *      generically -- it does not know whether the underlying model is
 *      NRLMSISE-00 (Earth), JB2008 (Earth, future), MCD (Mars, future)
 *      or a 4-parameter exponential. New models are added by writing
 *      a wrapper that conforms to the callback and registering it on
 *      the central body's spec at the application layer.
 *
 *   2. The space weather loader for the CelesTrak combined file
 *      (`sw19571001.txt`). This is the universal solar-activity input
 *      consumed by virtually every empirical thermosphere model on
 *      every body (F10.7 is also used as an EUV proxy for Mars MCD),
 *      so it lives here once rather than inside a specific model's
 *      wrapper. The parser follows the same "parse text once, hold in
 *      heap, query thread-safely via a per-thread handle" pattern as
 *      MappedEOPData.
 *
 * Time scale: queries take ET (seconds past J2000 TDB), same as the
 * integrator. Internally converted to UTC MJD (the CelesTrak file's
 * native time tag) via the same leap-second chain as `spody_eop.c`.
 *
 * SOURCE
 * ------
 *   CelesTrak combined space weather:
 *     https://celestrak.org/SpaceData/SW-Last5Years.csv  (recent)
 *     https://celestrak.org/SpaceData/sw19571001.txt    (full history)
 *   Format spec:
 *     https://celestrak.org/SpaceData/SpaceWx-format.php
 *
 *   The file is a daily-stepped text table from 1957-10-01 onward, with
 *   ~45-day predictions appended at the tail. Columns we consume are
 *   the daily F10.7 (observed + 81-day centered average) and the
 *   geomagnetic Ap (daily mean + eight 3-hour values). Kp values are
 *   ignored -- NRLMSISE-00 and the planned JB2008 wrapper both work
 *   off Ap. The header lines + "BEGIN OBSERVED" / "BEGIN PREDICTED"
 *   markers separate the observed and predicted regions; we record the
 *   boundary so callers can tell when a run window crosses into the
 *   predicted tail.
 */
#ifndef SPODY_ATMOSPHERE_H
#define SPODY_ATMOSPHERE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declaration: the atmosphere callback receives the full
 * ForceModelContext (mirrors the get_bf_rotation pattern) so each
 * model can pull whichever runtime inputs it needs (the shared
 * space weather table, an Earth-orientation handle for geodetic
 * conversion, etc.) without forcing the callback signature to grow
 * a new parameter every time a new model arrives. */
struct ForceModelContext;
typedef struct ForceModelContext ForceModelContext;

/* One parsed space weather record. Mirrors the CelesTrak columns we
 * actually need; the file ships extra columns (Cp, C9, ISN, F10.7-adj
 * variants) that no consumer reads today. F10.7 values are in solar
 * flux units (10^-22 W/m^2/Hz), Ap is the 3-hourly planetary index
 * in nominal 2-nT units, MJD is UTC. */
typedef struct {
    double mjd;                /* Modified Julian Date (UTC)                */
    double f107_obs;           /* daily F10.7 observed, sfu                 */
    double f107_obs_avg81;     /* 81-day centered avg of f107_obs, sfu      */
    double ap_daily;           /* daily mean planetary Ap                   */
    double ap_3h[8];           /* 3-hour Ap values (00,03,06,09,12,15,18,21)*/
} SpaceWeatherRecord;

/* Shared, read-only space weather table. Daily-stepped, contiguous (no
 * gaps in the source file). Thread-safe to share across workers because
 * every query is a binary search + indexed lookup that touches the
 * table read-only.
 *
 * `mjd_last_observed` is the MJD of the last row from the "OBSERVED"
 * section; everything after is CelesTrak prediction. Queries past
 * `mjd_last_predicted` fail (caller decides whether to abort or
 * extrapolate -- typical mission-design propagations beyond the
 * predicted horizon should switch to a climatological mean rather
 * than trust the tail). */
typedef struct {
    SpaceWeatherRecord *records;
    size_t              n_records;
    double              mjd_first;
    double              mjd_last_observed;
    double              mjd_last_predicted;
} MappedSpaceWeatherData;

/* Per-thread query handle. Holds the shared MappedSpaceWeatherData
 * pointer plus a tiny last-hit cache (record-index) so consecutive
 * queries at nearby ETs (typical integrator pattern) skip the
 * bsearch. */
typedef struct {
    const MappedSpaceWeatherData *msw;
    size_t cached_idx;
    int    cached_valid;
} MappedSpaceWeather;

/* Parse `filename` (a CelesTrak `sw19571001.txt` text file) into *msw.
 * Heap-allocates msw->records; caller must spody_free_MappedSpaceWeatherData
 * on shutdown. Returns 0 on success, -1 on parse / IO failure. */
int spody_setup_MappedSpaceWeatherData(MappedSpaceWeatherData *msw,
                                        const char *filename);

/* Release msw->records. Safe on a zero-initialised struct. */
int spody_free_MappedSpaceWeatherData(MappedSpaceWeatherData *msw);

/* Bind a per-thread MappedSpaceWeather handle to a shared
 * MappedSpaceWeatherData. Does not allocate; stores the pointer and
 * resets the cache. */
int spody_setup_MappedSpaceWeather(MappedSpaceWeather *map,
                                    const MappedSpaceWeatherData *msw);

/* Release the per-thread handle. No-op today; provided for symmetry. */
int spody_free_MappedSpaceWeather(MappedSpaceWeather *map);

/* Interpolate space weather inputs at Ephemeris Time `et` (seconds
 * past J2000 TDB). Daily values (f107_obs, f107_obs_avg81, ap_daily)
 * are linearly interpolated between bracketing days; the 3-hour Ap
 * values are picked by nearest-bin (NRLMSISE-00 prescribes the
 * current 3h Ap + four lagged bins separately, so a smoothing
 * interpolation would distort the model's intended use).
 *
 * Out parameters are written only on success; any may be NULL when
 * the caller doesn't need that quantity. `ap_3h_out` must point to
 * a 8-element array when non-NULL.
 *
 * Returns 0 on success, -1 when `et` falls outside the table's
 * coverage [mjd_first, mjd_last_predicted]. */
int spody_interpolate_space_weather(MappedSpaceWeather *map, double et,
                                     double *f107_obs_out,
                                     double *f107_obs_avg81_out,
                                     double *ap_daily_out,
                                     double  ap_3h_out[8]);

/* Space weather inputs in the exact layout NRLMSISE-00 prescribes
 * (see spody_nrlmsise00.h): daily values are taken from the UT-day
 * records WITHOUT interpolation, because that is how the model was
 * fit.
 *
 *   *f107_prev_out : observed daily F10.7 of the day BEFORE the
 *                    query day, sfu
 *   *f107a_out     : observed 81-day centered F10.7 average of the
 *                    query day, sfu
 *   ap_msis_out[7] : [0] daily Ap of the query day,
 *                    [1] 3h ap of the bin containing `et`,
 *                    [2..4] bins 3/6/9 hours earlier,
 *                    [5] mean of the eight bins 12..33 h earlier,
 *                    [6] mean of the eight bins 36..57 h earlier
 *                    (this is Fortran AP(1:7) with SW(9) = -1).
 *
 * Any output pointer may be NULL. Returns 0 on success, -1 when `et`
 * is outside the table or closer than 3 days to its start (the ap
 * history reaches back 57 hours). */
int spody_space_weather_msis_inputs(MappedSpaceWeather *map, double et,
                                     double *f107_prev_out,
                                     double *f107a_out,
                                     double  ap_msis_out[7]);

/* MJD of the last row in the OBSERVED section. Past this MJD the
 * table is CelesTrak's ~45-day prediction; the F10.7 / Ap predictions
 * are climatological estimates with no claim of accuracy beyond
 * ~10 days. The app uses this to warn the user when the run window
 * crosses into predictions. */
double spody_space_weather_last_observed_mjd(const MappedSpaceWeatherData *msw);

/* MJD of the last row (hard upper bound, predictions included).
 * Queries beyond this MJD fail. The app uses this to refuse a run
 * that would integrate past the table's predicted horizon, prompting
 * the user to update the file. */
double spody_space_weather_last_predicted_mjd(const MappedSpaceWeatherData *msw);

/* Atmosphere density callback. Returns mass density at the satellite
 * position in kg/m^3.
 *
 * Inputs:
 *   ctx       : ForceModelContext (read-only). The callback pulls the
 *               shared space weather handle, the Earth-orientation
 *               pipeline (for geodetic lat/lon conversion), etc.
 *               from ctx fields as the specific model needs them.
 *   et        : Ephemeris Time (seconds past J2000 TDB).
 *   r_bf_km   : satellite position in the central body's body-fixed
 *               frame, in kilometres.
 *
 * Outputs:
 *   *rho_kg_m3_out : mass density at the satellite location.
 *
 * Returns 0 on success, non-zero when the model could not be
 * evaluated (out-of-range altitude, missing space weather data, ...);
 * the drag force then defaults to zero acceleration. */
typedef int (*SpodyAtmosphereDensityFn)(
        const ForceModelContext *ctx, double et,
        const double r_bf_km[3], double *rho_kg_m3_out);

/* Atmosphere model handle. Wired onto ForceModelContext when the
 * central body has an atmosphere AND the user enabled the drag
 * force. A NULL pointer on the context means "no atmosphere":
 * `spody_force_drag` returns zero without complaint, so the drag
 * toggle can be left on the context unconditionally.
 *
 * `state` is an opaque pointer for model-specific data (cached
 * coefficient tables, lazy-loaded tabulated grids, ...). For
 * stateless models like the NRLMSISE-00 wrapper it stays NULL --
 * the wrapper reads everything it needs from ctx. */
typedef struct {
    SpodyAtmosphereDensityFn density;
    void                    *state;
} SpodyAtmosphere;

#ifdef __cplusplus
}
#endif

#endif /* SPODY_ATMOSPHERE_H */
