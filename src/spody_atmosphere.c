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
 * spody_atmosphere.c -- CelesTrak space weather parser + ET lookup.
 *
 * Atmosphere callback / SpodyAtmosphere struct have no code here --
 * they are pure interface, consumed by `spody_force_drag`
 * (spody_forcemodels.c) and implemented by per-model wrappers at the
 * application layer (the density model itself lives in the engine:
 * spody_nrlmsise00.c).
 *
 * Parser targets the CelesTrak CSV combined file (`SW-All.csv` /
 * `SW-Last5Years.csv`). The CSV header is stable and self-describing,
 * unlike the fixed-width legacy `sw19571001.txt` whose column offsets
 * have shifted across CelesTrak revisions. Format:
 *
 *   DATE,BSRN,ND,KP1..KP8,KP_SUM,AP1..AP8,AP_AVG,CP,C9,ISN,
 *   F10.7_OBS,F10.7_ADJ,F10.7_DATA_TYPE,F10.7_OBS_CENTER81,
 *   F10.7_OBS_LAST81,F10.7_ADJ_CENTER81,F10.7_ADJ_LAST81
 *
 * Spec: https://celestrak.org/SpaceData/SpaceWx-format.php
 *
 * We consume DATE -> MJD, AP1..AP8 -> ap_3h[8], AP_AVG -> ap_daily,
 * F10.7_OBS -> f107_obs, F10.7_OBS_CENTER81 -> f107_obs_avg81. The
 * F10.7_DATA_TYPE column ("OBS" / "INT" / "PRD") marks the boundary
 * between observed and predicted rows.
 *
 * Time scale: queries take ET (s past J2000 TDB); converted to UTC
 * MJD via spody_et_to_mjd_utc (spody_time.c), which owns the single
 * C-side leap-second table.
 */
#include "spody_atmosphere.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spody_const.h"
#include "spody_time.h"

/* Gregorian YYYY-MM-DD -> MJD (UTC midnight) via the shared Meeus
 * chain in spody_time.c. Returns -1.0 on malformed input. */
static double date_to_mjd(int year, int month, int day) {
    if (month < 1 || month > 12 || day < 1 || day > 31) return -1.0;
    return spody_greg_to_jd(year, month, day, 0, 0, 0.0) - JD_MJD_EPOCH;
}

/* Split a CSV line into up to `cap` field pointers (in-place: commas
 * are replaced with NULs). Returns the number of fields produced.
 * Empty fields are kept as zero-length strings, not skipped. */
static int split_csv(char *line, char *fields[], int cap) {
    int n = 0;
    char *p = line;
    fields[n++] = p;
    while (*p && n < cap) {
        if (*p == ',') { *p = '\0'; fields[n++] = p + 1; }
        ++p;
    }
    return n;
}

/* Parse a trimmed CSV field as a double; returns 1 on success, 0 on
 * blank / unparseable (out value untouched). */
static int csv_to_double(const char *s, double *out) {
    while (*s && isspace((unsigned char)*s)) ++s;
    if (!*s) return 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if (end == s) return 0;
    *out = v;
    return 1;
}

/* CSV column indices for the CelesTrak SW-All.csv schema. Keep
 * these in lockstep with the header line; if CelesTrak reorders
 * columns we'd want to switch to header-name-based lookup, but the
 * schema has been stable for >5 years. */
enum {
    COL_DATE     = 0,
    COL_AP1      = 12,
    /* AP2..AP8 = 13..19 */
    COL_AP_AVG   = 20,
    COL_F107_OBS         = 24,
    COL_F107_DATA_TYPE   = 26,
    COL_F107_OBS_CENTER81 = 27,
    COL_N_MIN    = 28,   /* minimum number of fields we require */
};

int spody_setup_MappedSpaceWeatherData(MappedSpaceWeatherData *msw,
                                        const char *filename) {
    if (!msw || !filename) return -1;
    memset(msw, 0, sizeof *msw);

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "spody_atmosphere: cannot open '%s': %s\n",
                filename, strerror(errno));
        return -1;
    }

    /* Two-pass: count rows that parse, then allocate exactly. The
     * file is ~3 MB / ~25k lines, a second pass is sub-ms. */
    char line[2048];
    size_t n = 0;
    int header_skipped = 0;
    while (fgets(line, sizeof line, fp)) {
        if (!header_skipped) { header_skipped = 1; continue; }
        if (line[0] == '\0' || line[0] == '\n' || line[0] == '\r') continue;
        ++n;
    }
    if (n == 0) {
        fprintf(stderr, "spody_atmosphere: no rows in '%s'\n", filename);
        fclose(fp);
        return -1;
    }

    msw->records = calloc(n, sizeof *msw->records);
    if (!msw->records) { fclose(fp); return -1; }
    msw->n_records = n;

    rewind(fp);
    header_skipped = 0;
    size_t k = 0;
    double mjd_last_obs = -INFINITY;
    while (k < n && fgets(line, sizeof line, fp)) {
        if (!header_skipped) { header_skipped = 1; continue; }
        /* Strip trailing CR/LF in place so the last field is clean. */
        size_t ll = strlen(line);
        while (ll && (line[ll-1] == '\n' || line[ll-1] == '\r')) line[--ll] = '\0';
        if (ll == 0) continue;

        char *fields[64];
        int nf = split_csv(line, fields, 64);
        if (nf < COL_N_MIN) continue;

        /* DATE: YYYY-MM-DD */
        int yy, mm, dd;
        if (sscanf(fields[COL_DATE], "%d-%d-%d", &yy, &mm, &dd) != 3) continue;
        double mjd = date_to_mjd(yy, mm, dd);
        if (mjd < 0.0) continue;

        SpaceWeatherRecord rec;
        memset(&rec, 0, sizeof rec);
        rec.mjd = mjd;

        for (int j = 0; j < 8; ++j) {
            double v = 0.0;
            if (csv_to_double(fields[COL_AP1 + j], &v)) rec.ap_3h[j] = v;
        }
        csv_to_double(fields[COL_AP_AVG],            &rec.ap_daily);
        csv_to_double(fields[COL_F107_OBS],          &rec.f107_obs);
        csv_to_double(fields[COL_F107_OBS_CENTER81], &rec.f107_obs_avg81);

        /* "OBS" rows are observed; "INT" (interpolated) and "PRD"
         * (predicted) are everything past the data cutoff. CelesTrak
         * emits the boundary in order so we can track it as we read. */
        const char *dt = fields[COL_F107_DATA_TYPE];
        while (*dt && isspace((unsigned char)*dt)) ++dt;
        if (dt[0] == 'O') {
            if (mjd > mjd_last_obs) mjd_last_obs = mjd;
        }

        msw->records[k++] = rec;
    }
    fclose(fp);

    if (k < msw->n_records) msw->n_records = k;
    if (k == 0) {
        free(msw->records);
        msw->records = NULL;
        msw->n_records = 0;
        return -1;
    }

    msw->mjd_first          = msw->records[0].mjd;
    msw->mjd_last_predicted = msw->records[k - 1].mjd;
    msw->mjd_last_observed  = isfinite(mjd_last_obs) ? mjd_last_obs : msw->mjd_first;
    return 0;
}

int spody_free_MappedSpaceWeatherData(MappedSpaceWeatherData *msw) {
    if (!msw) return -1;
    free(msw->records);
    msw->records   = NULL;
    msw->n_records = 0;
    return 0;
}

int spody_setup_MappedSpaceWeather(MappedSpaceWeather *map,
                                    const MappedSpaceWeatherData *msw) {
    if (!map || !msw) return -1;
    map->msw = msw;
    map->cached_idx = 0;
    map->cached_valid = 0;
    return 0;
}

int spody_free_MappedSpaceWeather(MappedSpaceWeather *map) {
    if (!map) return -1;
    map->msw = NULL;
    map->cached_idx = 0;
    map->cached_valid = 0;
    return 0;
}

/* Binary search: returns the index `i` such that records[i].mjd <= mjd
 * < records[i+1].mjd. Falls back to the last index when mjd >=
 * records[last].mjd and to 0 when mjd <= records[0].mjd. */
static size_t find_bracketing_index(const MappedSpaceWeatherData *msw,
                                    double mjd) {
    if (msw->n_records == 0) return 0;
    if (mjd <= msw->records[0].mjd) return 0;
    if (mjd >= msw->records[msw->n_records - 1].mjd) return msw->n_records - 1;
    size_t lo = 0, hi = msw->n_records - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (msw->records[mid].mjd <= mjd) lo = mid;
        else hi = mid;
    }
    return lo;
}

int spody_interpolate_space_weather(MappedSpaceWeather *map, double et,
                                     double *f107_obs_out,
                                     double *f107_obs_avg81_out,
                                     double *ap_daily_out,
                                     double  ap_3h_out[8]) {
    if (!map || !map->msw) return -1;
    const MappedSpaceWeatherData *msw = map->msw;
    if (msw->n_records == 0) return -1;

    double mjd = spody_et_to_mjd_utc(et);
    if (mjd < msw->mjd_first || mjd > msw->mjd_last_predicted) return -1;

    /* Use the per-thread cache first: the integrator typically queries
     * monotonically-increasing ET values, so the previous bracket is
     * still valid most of the time. */
    size_t i;
    if (map->cached_valid && map->cached_idx + 1 < msw->n_records
            && msw->records[map->cached_idx].mjd     <= mjd
            && msw->records[map->cached_idx + 1].mjd >  mjd) {
        i = map->cached_idx;
    } else {
        i = find_bracketing_index(msw, mjd);
        map->cached_idx = i;
        map->cached_valid = 1;
    }

    const SpaceWeatherRecord *r0 = &msw->records[i];
    const SpaceWeatherRecord *r1 = (i + 1 < msw->n_records) ? &msw->records[i + 1] : r0;
    double t = (r1 == r0) ? 0.0
                          : (mjd - r0->mjd) / (r1->mjd - r0->mjd);

    if (f107_obs_out)        *f107_obs_out        = r0->f107_obs        + t * (r1->f107_obs        - r0->f107_obs);
    if (f107_obs_avg81_out)  *f107_obs_avg81_out  = r0->f107_obs_avg81  + t * (r1->f107_obs_avg81  - r0->f107_obs_avg81);
    if (ap_daily_out)        *ap_daily_out        = r0->ap_daily        + t * (r1->ap_daily        - r0->ap_daily);

    /* 3-hour Ap: return all eight bins of the day where `et` falls
     * (the nearer of the two bracketing rows). NRLMSISE-00 wants a
     * specific lagged-window arrangement (current bin + four 3h
     * lags + two 12h/24h running means) that crosses day boundaries;
     * building that layout is the model wrapper's job -- this layer
     * just exposes the raw rows. */
    if (ap_3h_out) {
        const SpaceWeatherRecord *rd = (t < 0.5) ? r0 : r1;
        for (int j = 0; j < 8; ++j) ap_3h_out[j] = rd->ap_3h[j];
    }

    return 0;
}

int spody_space_weather_msis_inputs(MappedSpaceWeather *map, double et,
                                     double *f107_prev_out,
                                     double *f107a_out,
                                     double  ap_msis_out[7]) {
    if (!map || !map->msw) return -1;
    const MappedSpaceWeatherData *msw = map->msw;
    if (msw->n_records == 0) return -1;

    double mjd = spody_et_to_mjd_utc(et);
    /* No interpolation here (the model prescription is built on daily
     * values and UT-day-aligned 3h bins), so the whole last predicted
     * day is usable, not just its midnight. */
    if (mjd < msw->mjd_first || mjd >= msw->mjd_last_predicted + 1.0)
        return -1;

    double day = floor(mjd);
    size_t idx = (size_t)(day - msw->mjd_first);
    /* The prescription reaches back 57 h: require 3 full prior days
     * and verify the window is gap-free (daily-stepped records). */
    if (idx < 3 || idx >= msw->n_records) return -1;
    if (msw->records[idx].mjd != day || msw->records[idx - 3].mjd != day - 3.0)
        return -1;

    if (f107_prev_out) *f107_prev_out = msw->records[idx - 1].f107_obs;
    if (f107a_out)     *f107a_out     = msw->records[idx].f107_obs_avg81;

    if (ap_msis_out) {
        /* Global 3h-bin index since mjd_first; bin b of day i is
         * records[i].ap_3h[b] (UT-aligned: 00-03, 03-06, ...). */
        long g = (long)idx * 8 + (int)((mjd - day) * 8.0);
        double s;
        int k;
        ap_msis_out[0] = msw->records[idx].ap_daily;
        for (k = 0; k < 4; ++k) {   /* current bin + 3h/6h/9h lags */
            long b = g - k;
            ap_msis_out[1 + k] = msw->records[b / 8].ap_3h[b % 8];
        }
        s = 0.0;                    /* average of the 12..33 h lags  */
        for (k = 4; k <= 11; ++k) {
            long b = g - k;
            s += msw->records[b / 8].ap_3h[b % 8];
        }
        ap_msis_out[5] = s / 8.0;
        s = 0.0;                    /* average of the 36..57 h lags  */
        for (k = 12; k <= 19; ++k) {
            long b = g - k;
            s += msw->records[b / 8].ap_3h[b % 8];
        }
        ap_msis_out[6] = s / 8.0;
    }
    return 0;
}

double spody_space_weather_last_observed_mjd(const MappedSpaceWeatherData *msw) {
    if (!msw) return NAN;
    return msw->mjd_last_observed;
}

double spody_space_weather_last_predicted_mjd(const MappedSpaceWeatherData *msw) {
    if (!msw) return NAN;
    return msw->mjd_last_predicted;
}
