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

typedef struct { // DE440 has 15 triplets
    double start_epoch;
    double end_epoch;
    int days_per_record;
    int bytes_per_record;
    int number_coefficients_per_record;
    int location[15];
    int number_coefficients_per_component[15];
    int number_complete_sets_coefficients_per_record[15];
}EphemerisFile_Header;

typedef struct { 
    int record_number;
    int number_coefficients_per_record;
    double start_epoch;
    double end_epoch;
    double record[];
}EphemerisFile_Record;

typedef struct {
    MappedFile mf;
    EphemerisFile_Header *header;
    EphemerisFile_Record **records;
    size_t num_records;
    // single-shot cache per DE440 body index: if the same (idx, jd_epoch)
    // is requested again, return the cached position without re-evaluating Chebyshev.
    double cache_jd[EPH_CACHE_SLOTS];
    double cache_pos[EPH_CACHE_SLOTS][3];
    int cache_valid[EPH_CACHE_SLOTS];
} MappedEphemeris;

int spody_createfile_MappedEphemeris(const char *path, const char **file_names, const int n_files, const char *de);
int spody_setup_MappedEphemeris(MappedEphemeris *map, const char *filename);
int spody_setup_partialMappedEphemeris(MappedEphemeris *map, const char *filename, double in_start, double in_end);
int spody_get_ephposition(MappedEphemeris *map, int central_idx ,int target_idx, double jd_epoch, double result[3]);
/* Batch query: writes n_targets positions into a flat buffer of 3*n_targets
 * doubles, laid out as [x0,y0,z0, x1,y1,z1, ...]. Caller owns the buffer. */
int spody_get_ephposition_batch(MappedEphemeris *map, int central_idx, const int *target_idx_array, int n_targets, double jd_epoch, double *result);
int spody_get_lunarlibrationangles(MappedEphemeris *map, double jd_epoch, double result[3]);
    
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
