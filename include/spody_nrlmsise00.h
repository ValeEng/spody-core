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
 * NRLMSISE-00 empirical neutral-atmosphere model (surface to exosphere).
 *
 * Native C port, translated for SPODY directly from the official NRL
 * Fortran distribution (a U.S. Government work, public domain). The
 * translation keeps the Fortran's exact operation order and coefficient
 * values so the output is reproducible against the reference test
 * driver bundled with the original distribution; it differs from the
 * original only in being double precision and fully re-entrant (no
 * COMMON-block state: every call works on stack scratch, so the model
 * is safe to evaluate concurrently from OpenMP batch workers).
 *
 * SOURCE
 * ------
 *   Official Fortran distribution (NRLMSISE-00.FOR, "MSISE-00 01-FEB-02",
 *   includes the reference test driver used for validation):
 *     https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSISE-00/
 *   Model paper:
 *     J.M. Picone, A.E. Hedin, D.P. Drob, A.C. Aikin (2002),
 *     "NRLMSISE-00 empirical model of the atmosphere: Statistical
 *     comparisons and scientific issues", J. Geophys. Res. 107(A12),
 *     1468, doi:10.1029/2002JA009430
 *
 * SCOPE OF THE PORT
 * -----------------
 *   - MASS=48 only (all species + total mass density): the only mode
 *     the drag force needs, and the mode all reference test cases use.
 *     Single-species selection, the MASS=0 temperature-only mode, the
 *     GHP7 pressure-grid solver and the METERS unit toggle were not
 *     ported.
 *   - Switches: the model always runs with the standard configuration
 *     (all 25 variations ON). The one runtime choice kept from TSELEC
 *     is SW(9) = +1 / -1: daily Ap vs the 7-element 3-hour Ap history
 *     (`use_ap_array`). Research-grade per-variation switching would
 *     have to reintroduce the SW/SWC arrays.
 *
 * UNITS: native model units, identical to the Fortran default (IMR=0):
 * number densities in cm^-3, total mass density in g/cm^3, temperatures
 * in K. The drag wrapper converts rho to kg/m^3 (x 1000) at the boundary.
 */
#ifndef SPODY_NRLMSISE00_H
#define SPODY_NRLMSISE00_H

#ifdef __cplusplus
extern "C" {
#endif

/* Model inputs. UT, local solar time and longitude are independent
 * model coordinates; for a physically consistent query keep
 * lst_hr = sec/3600 + glon_deg/15 (the wrapper does).
 *
 * F10.7 conventions (from the model definition, easy to get wrong):
 *   f107  = observed (at Earth) daily F10.7 of the PREVIOUS day,
 *   f107a = observed 81-day average CENTERED on the query day,
 * both in solar flux units. Below ~80 km the model wants
 * f107 = f107a = 150 and ap = 4.
 *
 * ap_array layout when use_ap_array != 0 (Fortran AP(1:7)):
 *   [0] daily Ap
 *   [1] 3-hr ap for the current time
 *   [2] 3-hr ap 3 hours before
 *   [3] 3-hr ap 6 hours before
 *   [4] 3-hr ap 9 hours before
 *   [5] average of eight 3-hr ap, 12-33 hours prior
 *   [6] average of eight 3-hr ap, 36-57 hours prior */
typedef struct {
    int    doy;          /* day of year 1..366 (the year is not a model input) */
    double sec;          /* UT seconds of day                                  */
    double alt_km;       /* geodetic altitude, km                              */
    double glat_deg;     /* geodetic latitude, deg                             */
    double glon_deg;     /* geodetic longitude, deg                            */
    double lst_hr;       /* local apparent solar time, hours                   */
    double f107a;        /* 81-day centered average F10.7, sfu                 */
    double f107;         /* previous-day daily F10.7, sfu                      */
    double ap;           /* daily Ap (used when use_ap_array == 0)             */
    double ap_array[7];  /* 3-hr Ap history (used when use_ap_array != 0)      */
    int    use_ap_array; /* 0: daily Ap; 1: ap_array (Fortran SW(9) = -1)      */
} SpodyNrlmsise00Input;

/* d[0] He, d[1] O, d[2] N2, d[3] O2, d[4] Ar, d[6] H, d[7] N,
 * d[8] anomalous O -- number densities in cm^-3.
 * d[5] total mass density in g/cm^3: WITHOUT anomalous oxygen from
 * gtd7, WITH it ("effective total mass density for drag") from gtd7d.
 * t[0] exospheric temperature, t[1] temperature at alt_km, K. */
typedef struct {
    double d[9];
    double t[2];
} SpodyNrlmsise00Output;

/* Standard evaluation (Fortran GTD7). Thread-safe, no global state. */
void spody_nrlmsise00_gtd7(const SpodyNrlmsise00Input *inp,
                           SpodyNrlmsise00Output *out);

/* Drag variant (Fortran GTD7D): d[5] additionally includes the
 * anomalous-oxygen contribution, relevant above ~500 km. This is the
 * entry point the atmosphere callback uses. */
void spody_nrlmsise00_gtd7d(const SpodyNrlmsise00Input *inp,
                            SpodyNrlmsise00Output *out);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_NRLMSISE00_H */
