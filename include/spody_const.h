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
#define PI 3.14159265358979323846264338327950288419716939937511

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
}
#endif

#endif /* SPODY_CONST_H */