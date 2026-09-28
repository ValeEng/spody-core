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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "spody_io.h"
#include "spody_time.h"
#include "spody_const.h"

/* ============================================================
 * Text log mirror (see spody_io.h)
 * ============================================================ */
static FILE *g_log_mirror = NULL;

int spody_log_open_mirror(const char *path) {
    if (g_log_mirror) { fclose(g_log_mirror); g_log_mirror = NULL; }
    g_log_mirror = fopen(path, "w");
    return g_log_mirror ? 0 : -1;
}

void spody_log_close_mirror(void) {
    if (!g_log_mirror) return;
    /* The log is a copy of what was already printed: a failed final
     * flush (full disk) truncates the file but loses no result, so it
     * is reported on stderr and does not change the exit code. */
    int write_failed = ferror(g_log_mirror);
    if (fclose(g_log_mirror) != 0 || write_failed)
        fprintf(stderr, "warning: write failed on the log file (%s): "
                        "the log is incomplete\n", strerror(errno));
    g_log_mirror = NULL;
}

/* Write `fmt+args` to `term` always, and to the mirror if set. Splits
 * the va_list with va_copy so each consumer sees a fresh iterator. */
static void tee_vprintf(FILE *term, const char *fmt, va_list ap) {
    va_list ap2;
    va_copy(ap2, ap);
    vfprintf(term, fmt, ap);
    if (g_log_mirror) {
        vfprintf(g_log_mirror, fmt, ap2);
        /* No fflush per call -- spody_log_close_mirror flushes at the end. */
    }
    va_end(ap2);
}

void spody_log_printf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    tee_vprintf(stdout, fmt, ap);
    va_end(ap);
}

void spody_log_eprintf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    tee_vprintf(stderr, fmt, ap);
    va_end(ap);
}

/* ET -> "YYYY-DDDThh:mm:ss.ffffff" (ISO 8601 ordinal date, UTC), for
 * reading. Built from seconds past J2000 (resolution ~0.1 us), not
 * from the UTC MJD (a double of order 6e4 days resolves only ~0.6 us);
 * the MJD serves only for the leap offset. The ET printed next to it
 * is the exact value. */
static void utc_ordinal(double et, char out[32]) {
    double tt    = et - spody_tdb_minus_tt(et);
    double utc_s = tt + TT2TAI_SEC - spody_tai_minus_utc(spody_et_to_mjd_utc(et));
    double s0    = utc_s + 0.5 * SECONDSxDAY;       /* from 2000-01-01 00:00 */
    double day   = floor(s0 / SECONDSxDAY);
    double sec   = s0 - day * SECONDSxDAY;
    int year = 0, doy = 0;
    spody_mjd_to_doy((JD_J2000 - JD_MJD_EPOCH - 0.5) + day, &year, &doy, NULL);
    long long us = (long long)floor(sec * 1e6 + 0.5);
    int hh = (int)(us / 3600000000LL);
    int mm = (int)((us / 60000000LL) % 60);
    int ss = (int)((us / 1000000LL) % 60);
    int ff = (int)(us % 1000000LL);
    snprintf(out, 32, "%04d-%03dT%02d:%02d:%02d.%06d",
             year, doy, hh, mm, ss, ff);
}

void spody_log_time_anchor(const char *who, double et_first,
                           double et_last, size_t n_records) {
    char u0[32], u1[32];
    utc_ordinal(et_first, u0);
    utc_ordinal(et_last, u1);
    spody_log_printf("%s: time anchor t0 = ET %.17g s past J2000 TDB "
                     "[%a] = %s UTC\n", who, et_first, et_first, u0);
    spody_log_printf("%s: time column t = ET - t0; last record t = %.17g s, "
                     "ET %.17g [%a] = %s UTC; %zu records\n", who,
                     et_last - et_first, et_last, et_last, u1, n_records);
}
