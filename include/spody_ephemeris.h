#ifndef SPODY_EPHEMERIS_H
#define SPODY_EPHEMERIS_H

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

#define BUFFER_SIZE_EPH 256

    //debug
#define DEBUG_CHEBYSHEV 0
#define DEBUG_MYEPH 0

    // Conversion units
#define SECONDSxDAY 86400 


    //gravitational parameters KM^3/s^2
#define EARTH_MU 398600.4415 //398600.435507 //KM -->//(GRAV_CONST * EARTH_MASS) // m^3 s^-2
#define MOON_MU 4902.8005821478 //4902.798815861232 // 4902.800118 KM -->// (GRAV_CONST * MOON_MASS) // SUN_MU / ( 328900.56 * ( 1.0 + 81.30059 ) ) m^3 s^-2
#define SUN_MU  132712440017.99 //1.32712440018e11 //KM -->//(GRAV_CONST * SUN_MASS) // m^3 s^-2
#define MARS_MU 42828.314258067 //42828.31425807143 // 42828.375214 KM -->//(GRAV_CONST * MARS_MASS) // SUN_MU / 3098708.0 m^3 s^-2
#define MERCURY_MU 22032.080486418 //22032.08048642008//22032.080486196 KM -->//(GRAV_CONST * MERCURY_MASS) // SUN_MU / 6023600.0 m^3 s^-2
#define VENUS_MU 324858.59882646 //324858.5988264916//324858.592 KM -->//(GRAV_CONST * VENUS_MASS) // SUN_MU / 408523.71 m^3 s^-2
#define JUPITER_MU 126712767.85780 //1.267127678578078e8//1.26686534e8 KM -->//(GRAV_CONST * JUPITER_MASS) // SUN_MU / 1047.3486 m^3 s^-2
#define SATURN_MU 37940626.061137
#define URANUS_MU 5794549.0070719
#define NEPTUNE_MU 6836534.0638793
#define PLUTO_MU 981.60088770700

#define EMRAT 0.813005682214972154E+02 //TBD from JPL DE440

#ifdef __cplusplus
extern "C" {
#endif

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
} MappedEphemeris;

/*my epheermis JPL reader*/

    //file maniupaltion
int read_ephemeris_file_header(FILE *file, EphemerisFile_Header *ep );
int read_record_block(FILE *fp, EphemerisFile_Header *ep, EphemerisFile_Record *eprec);
int create_binary_ephemeris_file(EphemerisFile_Header *ep, int *old_epoch, const char *ascp_filename, const char *bin_filename);
double fds2cd(char *str);
    
    //position calculation
double chebyshev_evaluate(double time_scaled, const double *coefficients, int n_coeffs);
int calculate_body_position(MappedEphemeris *map, int target_idx, double jd_epoch, double result[3]);
    
    //mapping functions
int ephemeris_map_file(MappedEphemeris *map, const char *filename);
int ephemeris_unmap_file(MappedEphemeris *map);

    //wrapper
int spody_createfile_MappedEphemeris(const char *path, const char **file_names, const int n_files, const char *de);
int spody_setup_MappedEphemeris(MappedEphemeris *map, const char *filename);
int spody_setup_partialMappedEphemeris(MappedEphemeris *map, const char *filename, double in_start, double in_end);
int spody_get_position(MappedEphemeris *map, int central_idx ,int target_idx, double jd_epoch, double result[3]);
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
