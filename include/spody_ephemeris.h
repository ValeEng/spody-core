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
#ifndef SPODY_EPHEMERIS_H
#define SPODY_EPHEMERIS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include "spody_mapping.h"
#include "spody_const.h"

    //system libreries
#ifdef _WIN32
    #include <windows.h>
    #include <direct.h>
    #define MAKE_DIR(path) _mkdir(path)
    #define STAT_STRUCT struct _stat // stat is for file poperties 
    #define STAT_FUNC _stat
#else
    #include <unistd.h>
    #define MAKE_DIR(path) mkdir(path, 0755) //POSIX : we need to add the acces permission (0755 is standard)
    #define STAT_STRUCT struct stat
    #define STAT_FUNC stat
#endif

    //debug
#define DEBUG_CHEBYSHEV 0 // 0 = no debug | 1 = debug |---> CODE TESTING
#define DEBUG_EPHEMERIS 0 // 0 = no debug | 1 = debug |---> CODE TESTING

#define BUFFER_SIZE_EPH 256
#define EPH_CACHE_SLOTS 15

/* On-disk magic string (exactly 8 bytes). The .spody binary format stores
 * this at the very first bytes of the file so the loader can reject
 * mismatching formats with a clear error. Naming follows the project-wide
 * convention "SPDY" + 2-letter module ("EP" = ephemeris) + 2-letter variant.
 *   "SPDYEPET" -> epochs are Ephemeris Time (seconds past J2000 TDB)
 *   "SPDYEPJD" -> epochs are Julian Date (legacy, no longer produced) */
#define SPODY_EPH_MAGIC_LEN     8
#define SPODY_EPH_MAGIC_ET      "SPDYEPET"
#define SPODY_EPH_MAGIC_JD      "SPDYEPJD"
#define SPODY_EPH_FORMAT_VERSION 1u

/* On-disk header. The first 16 bytes (magic + version + reserved) are at
 * fixed offsets and let any reader identify the file kind without parsing
 * the rest. Epoch fields are interpreted according to the magic:
 *   - magic == SPDYEPET -> start_epoch, end_epoch in ET seconds past J2000
 *                          seconds_per_record = (DE440 days_per_record) * 86400
 *   - magic == SPDYEPJD -> start_epoch, end_epoch in Julian Date
 *                          seconds_per_record = days_per_record (in days, legacy)
 * The current spody_createfile_MappedEphemerisData writes only SPDYEPET. */
typedef struct {                                       /* DE440 has 15 triplets */
    char     magic[SPODY_EPH_MAGIC_LEN];               /* "SPDYEPET" / "SPDYEPJD" */
    uint32_t format_version;                           /* SPODY_EPH_FORMAT_VERSION */
    uint32_t reserved;                                 /* must be 0           */

    double start_epoch;                                /* ET seconds (SPDYEPET) or JD (SPDYEPJD) */
    double end_epoch;                                  /* same convention as start_epoch         */
    int seconds_per_record;                            /* duration of one record (s for SPDYEPET, d for SPDYEPJD) */
    int bytes_per_record;
    int number_coefficients_per_record;
    int location[15];
    int number_coefficients_per_component[15];
    int number_complete_sets_coefficients_per_record[15];
}EphemerisFile_Header;

typedef struct {
    int record_number;
    int number_coefficients_per_record;
    double start_epoch;                                /* ET seconds past J2000 (SPDEET format) */
    double end_epoch;
    double record[];
}EphemerisFile_Record;

/* Shared, read-only ephemeris data. The header is always a private heap
 * copy (its coverage epochs are reconciled against the records for subset
 * files -- the mapping itself is read-only); records point into the
 * memory-mapped file (or, for the partial setup, into private heap
 * copies). Neither is mutated by query calls, so a single
 * MappedEphemerisData can be safely shared across threads.
 *
 * Time scale: epochs in the file are already in ET (Ephemeris Time, seconds
 * past J2000 TDB) thanks to the SPDYEPET file format. The query API takes
 * ET directly as its time argument; no runtime conversion is needed. */
typedef struct {
    MappedFile mf;
    EphemerisFile_Header *header;
    EphemerisFile_Record **records;
    size_t num_records;
} MappedEphemerisData;

/* Per-thread query handle: holds a pointer to the shared (read-only)
 * MappedEphemerisData plus a private single-shot cache keyed by DE440
 * body index. If the same (idx, et) is requested again on the same
 * handle, the cached position is returned without re-evaluating Chebyshev.
 * Velocity queries fill cache_vel too; position-only queries leave the
 * velocity slot invalid (cache_vel_valid tracks it separately).
 *
 * Threading model:
 *   - Setup the MappedEphemerisData once (e.g. on the main thread).
 *   - Each worker thread declares its own MappedEphemeris bound to that
 *     shared MappedEphemerisData and calls the query functions with it.
 *   - Never share the same MappedEphemeris across threads. */
typedef struct {
    const MappedEphemerisData *med;
    double cache_jd[EPH_CACHE_SLOTS];
    double cache_pos[EPH_CACHE_SLOTS][3];
    double cache_vel[EPH_CACHE_SLOTS][3];
    int cache_valid[EPH_CACHE_SLOTS];
    int cache_vel_valid[EPH_CACHE_SLOTS];
} MappedEphemeris;

int spody_createfile_MappedEphemerisData(const char *path, const char **file_names, const int n_files, const char *de);
int spody_setup_MappedEphemerisData(MappedEphemerisData *med, const char *filename);
int spody_setup_partialMappedEphemerisData(MappedEphemerisData *med, const char *filename, double in_start, double in_end);
int spody_setup_MappedEphemeris(MappedEphemeris *map, const MappedEphemerisData *med);
int spody_free_MappedEphemeris(MappedEphemeris *map);
int spody_free_MappedEphemerisData(MappedEphemerisData *med);

/* Time argument: et = seconds past J2000 (TDB), the SPICE Ephemeris Time
 * convention. Convert from a Julian Date with ET_FROM_JD(jd). */
int spody_get_ephposition(MappedEphemeris *map, int central_idx, int target_idx, double et, double result[3]);
/* Velocity of target relative to central at et, ICRF km/s. Exact: the
 * analytic derivative of the position Chebyshev series (no finite
 * differences, no extra data). Returns 0 on success, -1 on unsupported
 * body ids. */
int spody_get_ephvelocity(MappedEphemeris *map, int central_idx, int target_idx, double et, double result[3]);
/* Full state [x,y,z,vx,vy,vz] of target relative to central at et,
 * ICRF km and km/s. The position half is bit-identical to
 * spody_get_ephposition; the velocity half to spody_get_ephvelocity.
 * Returns 0 on success, -1 on unsupported body ids. */
int spody_get_ephstate(MappedEphemeris *map, int central_idx, int target_idx, double et, double result[6]);
/* Batch query: writes n_targets positions into a flat buffer of 3*n_targets
 * doubles, laid out as [x0,y0,z0, x1,y1,z1, ...]. Caller owns the buffer. */
int spody_get_ephposition_batch(MappedEphemeris *map, int central_idx, const int *target_idx_array, int n_targets, double et, double *result);
int spody_get_lunarlibrationangles(MappedEphemeris *map, double et, double result[3]);
    
/* ICRF (J2000) -> Lunar body-fixed (PA frame)
 *
 * From:
    -The JPL Planetary and Lunar Ephemerides DE440 and DE441
    Ryan S. Park, William M. Folkner, James G. Williams, and Dale H. Boggs
    Published 2021 February 8 • © 2021 The Author(s). Published by The American Astronomical Society.
    The Astronomical Journal, Volume 161, Number 3
    Citation: Ryan S. Park et al 2021 AJ 161 105
    DOI: 10.3847/1538-3881/abd414
 *
 * Convention:
 *   C = Rz(psi) * Rx(theta) * Rz(phi)
 *
 * such that:
 *   r_PA = C * r_ICRF
 *
 * Angles in radians.
 */
void spody_getrotmatrix_icrf2moonpa(double phi, double theta, double psi, double R[3][3]);

/* Inverse transformation (PA -> ICRF)
 * Obtained as matrix transpose.
 */
void spody_getrotmatrix_moonpa2icrf(double phi, double theta, double psi, double R[3][3]);



#ifdef __cplusplus
}
#endif

#endif // SPODY_EPHEMERIS_H
