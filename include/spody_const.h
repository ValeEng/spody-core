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
#define JD_J2000     2451545.0 // Julian Date at the J2000 epoch (SPICE-style ET origin)
#define JD_MJD_EPOCH 2400000.5 // JD of MJD 0 (1858-11-17 00:00): MJD = JD - JD_MJD_EPOCH
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
#define CR_ALUMINUM 1.31


    // Conversion units
#define SECONDSxDAY 86400.0
#define KM2M 1000.0
#define M2KM 0.001
#define M2AU 6.684587122268445e-12
#define AU2M 149597870700.0
#define KM2AU 6.684587122268445e-9
#define AU2KM 149597870.7
#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)

    // Angle conversions (sub-degree). Used by the IAU 2006 Earth
    // orientation pipeline where amplitudes are published in
    // microarcseconds and EOP corrections in milliarcseconds.
#define ARCSEC2RAD (PI / (180.0 * 3600.0))         // 1" -> rad
#define MAS2RAD    (ARCSEC2RAD * 1.0e-3)           // 1 mas -> rad
#define UAS2RAD    (ARCSEC2RAD * 1.0e-6)           // 1 uas -> rad

    // Time scale: number of days in a Julian century (exact by
    // definition). Used to convert ET (s past J2000) to TT Julian
    // centuries, the argument expected by the IAU 2006 fundamental-
    // argument polynomials and the precession-nutation series.
#define DAYS_PER_JULIAN_CY 36525.0

    // GNSS time-scale offsets. GPS time is locked to TAI - 19 s
    // exactly since 1980-01-06, so the bridge to TT (= TAI + 32.184)
    // is a single constant -- no leap-second table required. GLONASS
    // broadcast TOC is UTC and needs the (TAI - UTC) leap chain,
    // hosted by spody_time.c (spody_tai_minus_utc).
#define GPST2TT_SEC      51.184          // TT - GPST (= 19 + 32.184)
#define TT2TAI_SEC       (-32.184)       // TAI - TT (exact)

    // TDB - TT periodic term (SPICE `deltet` algorithm):
    //   M = M0 + M1*ET, E = M + EB*sin(M), TDB-TT = K*sin(E).
    // The four DELTET/* values published in NAIF's leap-seconds
    // kernel (naif0012.tls); amplitude +/-1.657 ms. Implemented by
    // spody_tdb_minus_tt (spody_time.c).
#define DELTET_K    1.657e-3             // s     -- periodic amplitude
#define DELTET_EB   1.671e-2             //       -- eccentricity of the
                                         //          Earth-Moon barycenter orbit
#define DELTET_M0   6.239996             // rad   -- mean anomaly at J2000
#define DELTET_M1   1.99096871e-7        // rad/s -- mean-anomaly rate

    // GPS week structure. Week 0 starts 1980-01-06 00:00:00 GPST;
    // the half-week drives the week-rollover correction on tk.
#define GPS_WEEK_SEC        604800.0
#define HALF_GPS_WEEK_SEC   302400.0
#define GPS_EPOCH_JD_GPST   2444244.5    // 1980-01-06 00:00 GPST

    // Earth rotation rate. IERS Conventions 2010 sec. 1.4 nominal
    // sidereal value; identical to WGS-84 / GLONASS-ICD / GPS-ICD
    // omega_e at 1e-12. Polar-motion corrections to the instantaneous
    // axis are sub-arcsec and ignored where this constant is used
    // (GNSS broadcast ECEF -> ICRF velocity).
#define EARTH_ROT_RATE_RADPS 7.2921151467e-5


    // Gravitational parameters GM [km^3/s^2].
    // Alternative published values (DE440 headers, IAU 2015 nominal)
    // differ in the last digits; the ones below are the values the
    // engine has always integrated with. Do not swap without
    // re-running the validation suite against SPICE.
#define EARTH_MU   398600.4415        // alt: 398600.435507 (DE440)
#define MOON_MU    4902.8005821478    // alt: 4902.798815861232 (DE440)
#define SUN_MU     132712440017.99    // alt: 1.32712440018e11
#define MARS_MU    42828.314258067    // alt: 42828.375214 (system GM)
#define MERCURY_MU 22032.080486418
#define VENUS_MU   324858.59882646
#define JUPITER_MU 126712767.85780    // alt: 1.26686534e8 (planet-only GM)
#define SATURN_MU  37940626.061137
#define URANUS_MU  5794549.0070719
#define NEPTUNE_MU 6836534.0638793
#define PLUTO_MU   981.60088770700

    // Mean radii [km] from SPICE pck00011.
#define MOON_RADIUS    1737.4
#define EARTH_RADIUS   6378.1366
#define SUN_RADIUS     695700.0
#define MARS_RADIUS    3376.20
#define MERCURY_RADIUS 2440.53
#define VENUS_RADIUS   6051.8
#define JUPITER_RADIUS 71492.0
#define SATURN_RADIUS  60268.0
#define URANUS_RADIUS  25559.0
#define NEPTUNE_RADIUS 24764.0
#define PLUTO_RADIUS   1195.0

    // WGS-84 reference ellipsoid (NGA.STND.0036). Used for the
    // ECEF -> geodetic conversion feeding the Earth atmosphere model
    // (NRLMSISE-00 takes geodetic latitude/altitude). Distinct from
    // EARTH_RADIUS (SPICE mean radius, used for impact/eclipse
    // geometry) on purpose -- do not merge them.
#define WGS84_A_KM  6378.137
#define WGS84_INV_F 298.257223563

#define EMRAT 0.813005682214972154E+02 //TBD from JPL DE440

    // Primary-pair mean separations used by the CR3BP model.
    // The standard CR3BP assumes a fixed circular orbit between the
    // two primaries; L is the radius of that circle. Values are the
    // canonical literature ones for Earth-Moon (sub-AU mean separation
    // averaged over a synodic month). Sun-X pairs use AU2KM scaled by
    // the planet's semi-major axis in AU and are added here as needed.
#define EARTH_MOON_DISTANCE_KM 384400.0

#ifdef __cplusplus
}
#endif

#endif /* SPODY_CONST_H */