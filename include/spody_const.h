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
#ifndef SPODY_CONST_H
#define SPODY_CONST_H

#ifdef __cplusplus
extern "C" {
#endif


    // Constants
#define JD_JAN_2000 2451545.0 // Julian Date at J2000 epoch (alias: JD_J2000)
#define JD_J2000    2451545.0 // explicit name for the SPICE-style ET origin
#define PI 3.14159265358979323846264338327950288419716939937511
#define SOLAR_LUMINOSITY_4PIC 1.016111489628894e+08 // Watts divided 4pi*c*10^9 | we can use r in km and am in m^2/kg

    // SPICE-style Ephemeris Time:
    //   ET = seconds (TDB) past J2000 (1 Jan 2000 12:00:00 TDB).
    // ET is the canonical internal time scale for SPODY: same units as the
    // integrator (seconds), centred near today (magnitude ~8e8 in 2026 ->
    // ULP ~183 ns, ~250x better than full JD's 47 us). Negative ET values
    // are valid (epoch before J2000) and pose no precision concerns.
#define ET_FROM_JD(jd) ( ((jd) - JD_J2000) * SECONDSxDAY )
#define JD_FROM_ET(et) ( JD_J2000 + (et) / SECONDSxDAY )


    //Radiation pressure coefficients (from Montenbruck and Gill)
#define CR_ALUMINUM_COATED_MYLAR 1.88 
#define CR_SOLARPANEL 1.21
#define CR_ALUMNINUM 1.31


    // Conversion units
#define SECONDSxDAY 86400 
#define KM2M 1000.0
#define M2KM 0.001
#define M2AU 6.684587122268445e-12
#define AU2M 149597870700
#define KM2AU 6.684587122268445e-9
#define AU2KM 149597870.7
#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)


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
}
#endif

#endif /* SPODY_CONST_H */