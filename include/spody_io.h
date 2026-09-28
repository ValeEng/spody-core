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
#ifndef SPODY_IO_H
#define SPODY_IO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdio.h>

/* ============================================================
 * Text log mirror: tee stdout / stderr to a file
 *
 * Every message the library or its host prints goes through
 * spody_log_printf (progress, notices -> stdout) or spody_log_eprintf
 * (warnings, errors -> stderr). Both always write the terminal; after
 * spody_log_open_mirror they also write the same text to the mirror
 * file, so a run's saved log holds the library's own diagnoses (a
 * damaged ephemeris, a malformed harmonics row) next to the host's
 * lines. Without a mirror the terminal output is exactly a bare
 * printf / fprintf(stderr, ...).
 *
 * One mirror per process. Open and close it from a single thread
 * (before and after any parallel region); the printf calls themselves
 * are safe from worker threads, since each one is a single vfprintf
 * per stream and the C runtime locks a FILE per call, so lines never
 * interleave mid-line. Debug-build traces (#if DEBUG_*) stay bare
 * printf on purpose.
 * ============================================================ */

/* Open `path` for write (truncate). Closes any previous mirror first.
 * Returns 0 on success, -1 on fopen failure. */
int  spody_log_open_mirror(const char *path);

/* Close the current mirror, if any. Idempotent. The final flush happens
 * here: if it fails (full disk) a warning goes to stderr -- the log is a
 * copy of the terminal, so no result is lost and nothing else changes. */
void spody_log_close_mirror(void);

/* gcc/clang check every call's arguments against its format string, as
 * they do for printf itself; MSVC has no equivalent attribute. */
#if defined(__MINGW32__)
/* MinGW's `printf` archetype is the legacy msvcrt one (no %zu); the
 * UCRT runtime the build links does handle C99 formats. */
#define SPODY_PRINTF_FMT __attribute__((format(gnu_printf, 1, 2)))
#elif defined(__GNUC__) || defined(__clang__)
#define SPODY_PRINTF_FMT __attribute__((format(printf, 1, 2)))
#else
#define SPODY_PRINTF_FMT
#endif

/* printf-style write to stdout (and to the mirror, if open). */
void spody_log_printf (const char *fmt, ...) SPODY_PRINTF_FMT;

/* printf-style write to stderr (and to the mirror, if open). */
void spody_log_eprintf(const char *fmt, ...) SPODY_PRINTF_FMT;

/* Time anchor of a binary with a relative time column (SPDYOUT_):
 * the ET of t = 0 and of the last record, each as %.17g (round-trips
 * to the same double) and as a hex float (the bits themselves), plus
 * the UTC reading to the microsecond. A converter calls it once on
 * success, so its log is enough to recover every absolute epoch as
 * ET = t0 + t. `who` prefixes the lines ("sp3", "gps", ...). */
void spody_log_time_anchor(const char *who, double et_first,
                           double et_last, size_t n_records);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_IO_H */
