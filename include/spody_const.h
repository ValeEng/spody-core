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
#define MINUTESxDAY 1440.0
#define KM2M 1000.0
#define M2KM 0.001
#define M2AU 6.684587122268445e-12
#define AU2M 149597870700.0
#define KM2AU 6.684587122268445e-9
#define AU2KM 149597870.7
#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)
#define REVDAY2RADMIN (PI/720.0) // 2*PI/(24*60) = 2*PI/(1440) = PI/720

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

    // Greenwich Mean Sidereal Time, IAU 1982 (Astronomical Almanac):
    // seconds of time as a cubic in Julian centuries of UT1 past
    // J2000, 240 seconds of time to the degree. The equinox-based
    // sibling of the Earth Rotation Angle, which is the same physical
    // angle measured from the CIO instead; the two differ by the
    // equation of the origins. Evaluated flat rather than nested:
    // Horner's form moves the result by 1.8e-12 rad, and the
    // implementations this is checked against do not use it.
#define GMST_C0           67310.54841
#define GMST_C1           (876600.0 * 3600.0 + 8640184.812866)
#define GMST_C2           0.093104
#define GMST_C3          -6.2e-6
#define GMST_SEC_PER_DEG  240.0

    // GNSS time-scale offsets. GPS time is locked to TAI - 19 s
    // exactly since 1980-01-06, so the bridge to TT (= TAI + 32.184)
    // is a single constant -- no leap-second table required. GLONASS
    // broadcast TOC is UTC and needs the (TAI - UTC) leap chain,
    // hosted by spody_time.c (spody_tai_minus_utc).
#define GPST2TT_SEC      51.184          // TT - GPST (= 19 + 32.184)
#define TT2TAI_SEC       (-32.184)       // TAI - TT (exact)
    // BeiDou time (BDT) started at 2006-01-01 00:00:00 UTC, when
    // GPST - UTC was 14 s: GPST = BDT + 14 s exactly. GLONASS system
    // time is UTC(SU) + 3 h, leap seconds included (RINEX / SP3-d).
#define BDT2GPST_SEC     14.0            // GPST - BDT (exact)
#define GLOT_MINUS_UTC_SEC 10800.0       // GLONASS time - UTC (3 h)

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
    // omega_e at 1e-12. Used where the ICD prescribes it (the GPS
    // broadcast orbit algorithm) and as the spin of the co-rotating
    // atmosphere in the drag; an ECEF -> ICRF velocity instead takes the
    // true rotation of the ITRS chain (spody_bf_angular_velocity_icrf):
    // polar motion tilts the axis by tenths of an arcsecond, ~1 mm/s.
#define EARTH_ROT_RATE_RADPS 7.2921151467e-5

    // TIO locator s'(t) = -47 uas per Julian century (IERS TN 36 eq.
    // 5.13). It appears twice: inside the polar-motion matrix, and as
    // the one term separating the classical pseudo-earth-fixed frame
    // from the CIO-based terrestrial intermediate one.
#define TIO_LOCATOR_UAS_PER_CY  (-47.0)

    // Gravitational parameters GM [km^3/s^2].
    // Alternative published values (DE440 headers, IAU 2015 nominal)
    // differ in the last digits; the ones below are the values the
    // engine has always integrated with. Do not swap without
    // re-running the validation suite against SPICE.
    // A central body with a harmonics file does not use these: its
    // two-body term takes the GM of the file, the one its normalized
    // coefficients were estimated with. The values below serve the
    // central body without harmonics, the third bodies and CR3BP.
#define EARTH_MU   398600.4415        // = EIGEN-6C4 GM (TT-compatible); alt: 398600.435507 (DE440)
#define MOON_MU    4902.8001184575496 // DE440 BODY301_GM (GRAIL era); was 4902.8005821478 (DE405, GMAT default)
#define SUN_MU     132712440017.99    // alt: 1.32712440018e11
#define MARS_MU    42828.314258067    // alt: 42828.375214 (system GM)
#define MERCURY_MU 22032.080486418
#define VENUS_MU   324858.59882646
#define JUPITER_MU 126712767.85780    // alt: 1.26686534e8 (planet-only GM)
#define SATURN_MU  37940626.061137
#define URANUS_MU  5794549.0070719
#define NEPTUNE_MU 6836534.0638793
#define PLUTO_MU   981.60088770700

    // Radii [km] from SPICE pck00011 (BODYnnn_RADII): the EQUATORIAL
    // radius of every body (Earth mean 6371.0084). They size the bodies
    // drawn by the GUI and, with force_model.body_shape =
    // "equatorial_sphere", the occulting discs, the IMPACT sphere and
    // the altitude reference.
#define MOON_RADIUS    1737.4
#define EARTH_RADIUS   6378.1366
#define SUN_RADIUS     695700.0
#define MARS_RADIUS    3396.19
#define MERCURY_RADIUS 2440.53
#define VENUS_RADIUS   6051.8
#define JUPITER_RADIUS 71492.0
#define SATURN_RADIUS  60268.0
#define URANUS_RADIUS  25559.0
#define NEPTUNE_RADIUS 24764.0
#define PLUTO_RADIUS   1195.0

    // Polar radii [km], pck00011 (third value of BODYnnn_RADII), for the
    // bodies the kernel models as oblate spheroids. With body_shape =
    // "ellipsoid" (the default) shadow, IMPACT and altitude use the
    // spheroid of equatorial x polar radius. The Moon, Venus, the Sun
    // and Pluto are spheres in pck00011 and have no entry.
#define EARTH_RADIUS_POLAR   6356.7519
#define MARS_RADIUS_POLAR    3376.20
#define MERCURY_RADIUS_POLAR 2438.26
#define JUPITER_RADIUS_POLAR 66854.0
#define SATURN_RADIUS_POLAR  54364.0
#define URANUS_RADIUS_POLAR  24973.0
#define NEPTUNE_RADIUS_POLAR 24341.0

    // Spin-axis direction in ICRF (IAU WGCCRE, as in pck00011
    // BODYnnn_POLE_RA / _POLE_DEC): RA = RA0 + RA1*T, DEC = DEC0 + DEC1*T,
    // degrees, T in Julian centuries TDB from J2000. It orients the
    // spheroids above. The periodic terms of Neptune (0.7 deg) and the
    // tiny ones of Jupiter are left out; the Earth's pole so obtained
    // misses nutation (9 arcsec), which moves its limb by < 1 m.
#define EARTH_POLE_RA0      0.0
#define EARTH_POLE_RA1     -0.641
#define EARTH_POLE_DEC0    90.0
#define EARTH_POLE_DEC1    -0.557
#define MARS_POLE_RA0     317.269202
#define MARS_POLE_RA1      -0.10927547
#define MARS_POLE_DEC0     54.432516
#define MARS_POLE_DEC1     -0.05827105
#define MERCURY_POLE_RA0  281.0103
#define MERCURY_POLE_RA1   -0.0328
#define MERCURY_POLE_DEC0  61.4155
#define MERCURY_POLE_DEC1  -0.0049
#define JUPITER_POLE_RA0  268.056595
#define JUPITER_POLE_RA1   -0.006499
#define JUPITER_POLE_DEC0  64.495303
#define JUPITER_POLE_DEC1   0.002413
#define SATURN_POLE_RA0    40.589
#define SATURN_POLE_RA1    -0.036
#define SATURN_POLE_DEC0   83.537
#define SATURN_POLE_DEC1   -0.004
#define URANUS_POLE_RA0   257.311
#define URANUS_POLE_RA1     0.0
#define URANUS_POLE_DEC0  -15.175
#define URANUS_POLE_DEC1    0.0
#define NEPTUNE_POLE_RA0  299.36
#define NEPTUNE_POLE_RA1    0.0
#define NEPTUNE_POLE_DEC0  43.46
#define NEPTUNE_POLE_DEC1   0.0

    // NAIF id of the Sun. Spelled out here rather than inline because
    // two unrelated places special-case it: the SRP source direction
    // (ephemeris target) and the SRP occulter list, which must never
    // contain the light source itself.
#define SUN_NAIF 10

    // NAIF id of the Earth.
#define EARTH_NAIF 399

    // General relativity, Schwarzschild term of the central body (IERS
    // Conventions 2010 sec. 10.3, eq. 10.12, first line). The PPN
    // parameters are 1 in general relativity.
#define SPEED_OF_LIGHT_KMS  299792.458      // km/s, exact (SI definition)
#define PPN_BETA            1.0
#define PPN_GAMMA           1.0

    // Earth radiation pressure: sunlight reflected by the Earth (albedo)
    // and its thermal emission (infrared), Knocke, Ries & Tapley (1988),
    // "Earth radiation pressure effects on satellites", AIAA 88-4292.
    // Latitude phi and season t (seconds from KNOCKE_T0_ET):
    //   albedo     a = A0 + (C0 + C1 cos wt + C2 sin wt) P1(sin phi) + A2 P2(sin phi)
    //   emissivity e = E0 + (K0 + K1 cos wt + K2 sin wt) P1(sin phi) + E2 P2(sin phi)
    // with w = 2 pi / Julian year. The original paper was not
    // accessible; values as in Orekit's KnockeRediffusedForceModel,
    // which match those quoted in the literature citing Knocke.
#define KNOCKE_A0   0.34
#define KNOCKE_C0   0.0
#define KNOCKE_C1   0.10
#define KNOCKE_C2   0.0
#define KNOCKE_A2   0.29
#define KNOCKE_E0   0.68
#define KNOCKE_K0   0.0
#define KNOCKE_K1   (-0.07)
#define KNOCKE_K2   0.0
#define KNOCKE_E2   (-0.18)
    // Reference epoch of the seasonal term: 1981-12-22, taken at 00:00
    // UTC as Orekit does (the paper gives the day only); ET from the
    // SpOdy time chain (spopy.time.utc_to_et).
#define KNOCKE_T0_ET        (-568900747.8163546)
#define JULIAN_YEAR_S       31557600.0      // 365.25 d
    // Quadrature over the Earth's disk as seen from the satellite:
    // dA cos(theta)/d^2 is the solid angle of the element, so the
    // irradiance is (M/pi) * integral of dOmega. Gauss-Legendre in the
    // cosine of the nadir angle (6 nodes; the 3 positive ones and their
    // weights, the others are symmetric) x 12 uniform azimuths = 72 rays.
    // A uniformly bright sphere comes out exact; a realistic albedo with
    // the terminator in view to ~3e-4 at every altitude, the same as
    // 8 x 16 rays (the terminator step sets the floor).
#define EARTHRAD_N_NADIR    6               // Gauss-Legendre nodes, GL6_* below
#define EARTHRAD_N_AZIMUTH  12
#define GL6_X1  0.23861918608319693
#define GL6_X2  0.66120938646626448
#define GL6_X3  0.93246951420315205
#define GL6_W1  0.46791393457269137
#define GL6_W2  0.36076157304813894
#define GL6_W3  0.17132449237916975

    // NAIF id of the Moon: a tide-raising body for the Earth, and the
    // central body whose own tide the Earth raises.
#define MOON_NAIF 301

    // Solid-body tides, frequency-independent part (IERS Conventions
    // 2010, sec. 6.2.1, eq. 6.6-6.7, Table 6.3, anelastic Earth).
    // k_nm = RE + i*IM; the imaginary parts carry the anelastic lag.
    // Degree 3 is the Moon's tide only in practice (the Sun's is
    // (R/r)^1 = 4e-5 of its degree 2); K2P_ are the k^(+)_2m that let a
    // degree-2 forcing deform degree 4 through the Earth's flattening
    // and rotation.
#define TIDE_EARTH_K20_RE  0.30190
#define TIDE_EARTH_K21_RE  0.29830
#define TIDE_EARTH_K21_IM  (-0.00144)
#define TIDE_EARTH_K22_RE  0.30102
#define TIDE_EARTH_K22_IM  (-0.00130)
#define TIDE_EARTH_K30     0.093
#define TIDE_EARTH_K31     0.093
#define TIDE_EARTH_K32     0.093
#define TIDE_EARTH_K33     0.094
#define TIDE_EARTH_K2P0    (-0.00089)
#define TIDE_EARTH_K2P1    (-0.00080)
#define TIDE_EARTH_K2P2    (-0.00057)
    // Permanent part of the degree-2 zonal tide, (dC20)perm = A0*H0*k20
    // (IERS 2010 eq. 6.14, k20 = 0.30190): what a zero-tide field
    // already holds and the step-1 correction must not add again.
#define TIDE_EARTH_A0H0    (-1.39141e-8)    // 4.4228e-8 * (-0.31460)
    // Lunar Love number k2 (degree 2, one value for every order):
    // GRGM1200B, GRAIL, estimated together with the field (PDS label
    // gggrx_1200b_sha.lbl: 0.024116 +/- 0.000108). That field is
    // tide-free: the full tide has to be added (LRO fits, 2026-09).
#define TIDE_MOON_K2       0.024116

    // Half-interval of the rotation R(t+h) R(t-h)^T from which
    // spody_bf_angular_velocity_icrf reads a body's angular velocity
    // (axis-angle form, exact for a fixed axis). 60 s: long enough that
    // the rounding of the matrix entries (the Earth's UT1 date resolves
    // ~0.6 us, 4.6e-11 rad of ERA) stays below 1e-15 rad/s after the
    // division by 2h, short enough that omega barely changes across it
    // (h = 30 and 300 s agree with 60 s to 1e-16 rad/s for the Earth).
#define SPODY_BF_OMEGA_FD_STEP_S 60.0

    // Node spacing of the IAU 2006 (X, Y, s) interpolation grid in
    // spody_iau2006_xys_interp, in seconds.
    //
    // The series costs ~70 us to evaluate -- 3084 terms, each a 14-term
    // linear combination plus a sine and a cosine -- and the Earth force
    // model needs it at every RHS evaluation, 78k times over a week-long
    // GNSS propagation. X and Y are smooth, moving about 650 and 340 mas
    // respectively over eight days, so a cubic through hourly nodes
    // reproduces them to 3e-7 mas: roughly 40 nanometres at GNSS radius,
    // six orders of magnitude below the millimetre these runs care about
    // and two orders below the uncertainty of the EOP inputs themselves.
    //
    // ERA and polar motion are deliberately NOT interpolated. ERA advances
    // 15 arcsec per second, so interpolating it would be the one thing
    // here that actually costs accuracy, and both are cheap closed forms
    // anyway.
#define SPODY_XYS_NODE_S 3600.0

    // WGS-84 reference ellipsoid (NGA.STND.0036). Used for the
    // ECEF -> geodetic conversion feeding the Earth atmosphere model
    // (NRLMSISE-00 takes geodetic latitude/altitude). Distinct from
    // EARTH_RADIUS (SPICE equatorial radius, used for impact/eclipse
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

    // `spody calibrate` sliding-window density-scale fit (app-side
    // subcommand; values centralised here per the all-constants-in-
    // spody_const.h rule). The window default trades locality of the
    // k(t) estimate against per-window drag signal; 24 h keeps a full
    // day of in-track drift in each least-squares fit. Windows whose
    // sample count falls below the minimum are folded into their
    // predecessor. The delta-rms floor rejects fits where the drag
    // on/off trajectories differ by less than ~1 mm rms in-track --
    // there k is unobservable and the normal equation is noise.
#define SPODY_CAL_WINDOW_DEFAULT_H    24.0    // fit-window length, hours
#define SPODY_CAL_MIN_WINDOW_SAMPLES  6       // fewer -> fold into previous window
#define SPODY_CAL_MIN_DELTA_RMS_KM    1.0e-6  // drag-signal rms floor (1 mm)

    // Adaptive harmonics degree (spody_hg_adaptive_degree). The rule
    // is N(r) = LN_INV_EPS / ln(r / R_ref), i.e. the degree at which
    // the geometric decay (R_ref/r)^n drops the term below a relative
    // accuracy of exp(-LN_INV_EPS).
    //
    // 38 sits deliberately above ln(1/DBL_EPSILON) = 36.04, the value
    // "indistinguishable in double precision" would suggest. The rule
    // ignores the decay of the coefficients themselves, so it is an
    // upper bound only as long as that omission dominates; the extra
    // margin keeps it one at radii where the two effects are close.
    // It costs a few percent of the achievable saving, which is the
    // right trade for a truncation that must never remove signal.
    //
    // The step margin scales the |v|*h excursion used to bound the
    // lowest radius reachable within one integrator step, so the
    // degree is picked for that radius rather than for the one the
    // step starts from. 1.0 = the full straight-line excursion.
#define SPODY_HG_ADAPTIVE_LN_INV_EPS  38.0
#define SPODY_HG_ADAPTIVE_STEP_MARGIN 1.0

    // Discontinuity stops (spody_next_force_discontinuity): a step lands
    // this far before a known jump of the force model and a step of at
    // most twice this crosses it, so the jump falls inside an interval
    // where the step error it causes, ~ |jump| * eps^2, is negligible.
#define SPODY_DISC_STOP_EPS_S 1.0e-3

    // SGP4 constants. (WGS-72)
    // Source: Hoots & Roehrich, Spacetrack Report No. 3 (1980)
#define SGP4_WGS72_MU   398600.8         // km^3/s^2
#define SGP4_WGS72_RE   6378.135         // km
#define SGP4_WGS72_J2   0.001082616
#define SGP4_WGS72_J3  -0.00000253881
#define SGP4_WGS72_J4  -0.00000165597
    // Reference mean motion [Earth radii^1.5 / min], defined as
    // 60/sqrt(SGP4_WGS72_RE^3/SGP4_WGS72_MU). The literal below is that
    // expression in double precision, not the nine-digit
    // 0.0743669161 printed in the report: rounding a derived quantity
    // would leave ke, mu and Re mutually inconsistent. Recomputing
    // the formula and recovering those nine digits is the check that
    // this whole block was transcribed correctly.
#define SGP4_WGS72_KE   0.07436691613317342
#define SGP4_WGS72_K2     (0.5 * SGP4_WGS72_J2)
#define SGP4_WGS72_K4     (-0.375 * SGP4_WGS72_J4)
#define SGP4_WGS72_J3OJ2  (SGP4_WGS72_J3 / SGP4_WGS72_J2)
    // Lane & Cranford power-law density profile: the two anchor
    // altitudes above the surface [km]. The quantities actually used,
    // ((q0-s)/RE)^4 and 1+s/RE, are perigee-dependent and therefore
    // computed during initialisation, not here.
#define SGP4_Q0_KM  120.0
#define SGP4_S_KM   78.0

#define SGP4_DEEP_SPACE_PERIOD_MIN   225.0
#define SGP4_DEEP_SPACE_PERIOD_DAY  0.15625

#define SGP4_SIMPLE_DRAG_PERIGEE_KM  220.0
#define SGP4_S_MOD_PERIGEE_KM        156.0
#define SGP4_S_LOW_PERIGEE_KM         98.0
#define SGP4_S_FLOOR_KM               20.0

#define SGP4_KEPLER_TOL       1.0e-12  /* STR#6 value; STR#3 used 1e-6 */
#define SGP4_KEPLER_MAX_ITER  10       /* DO 130 I=1,10                */
#define SGP4_KEPLER_MAX_STEP  0.95     /* Crawford (1995)              */

    // Perturbed eccentricity, once drag has driven it past zero. Below the
    // trap the model stops; between the trap and zero it is held at the
    // floor. Both are fitted to the published vectors, which bracket the
    // trap between -9.2e-4 and -1.21e-3 and put a four-decade minimum on
    // the floor at 1e-6.
#define SGP4_ECC_TRAP        -0.001
#define SGP4_ECC_FLOOR        1.0e-6

    // Below this, the two drag coefficients that divide by e0 (C3 and
    // xmcof) are zeroed: the perigee of a circular orbit is not defined,
    // and zero is the limit of what they contribute, not an approximation
    // of it. Unlike the two above, this threshold is declared by the
    // source rather than measured -- only one published case (28057,
    // e0 = 8.84e-5) falls below it.
#define SGP4_LOW_ECC  1.0e-4

    // Refusal threshold on 1 + cos i0, the denominator of xlcof. Not a
    // physical limit: cos i0 carries an absolute error of about one ulp of
    // 1.0, so below this the difference has fewer than four significant
    // digits and is noise however it is computed. It sits 0.29 arcsec from
    // 180 deg, nine orders of magnitude past the most inclined object in
    // any catalogue.
#define SGP4_MIN_1P_COSI  1.0e-12 

    // Deep space: the epochs the lunar and solar mean elements are
    // tabulated against. DS50 counts days from 1950 Jan 0.0; the
    // polynomials below count from 1900 Jan 0.5.
#define SGP4_MJD_1950          33281.0
#define SGP4_DAY_1950_TO_1900  18261.5

    // Lunar and solar mean elements, linear in days since 1900 Jan 0.5
    // [rad, rad/day]. Every rate is a period that checks on its own:
    // the node regresses in 18.61 y, the mean longitude runs the
    // 27.322 d sidereal month, perigee circulates in 8.85 y, and the
    // solar anomaly is the 365.26 d year. Longitude minus perigee
    // reproduces ZNL to eight digits, which is how this block was
    // verified without leaving the page.
#define SGP4_MOON_NODE_0     4.5236020
#define SGP4_MOON_NODE_DOT  (-9.2422029e-4)
#define SGP4_MOON_LON_0      4.7199672
#define SGP4_MOON_LON_DOT    0.22997150
#define SGP4_MOON_PERI_0     5.8351514
#define SGP4_MOON_PERI_DOT   0.0019443680
#define SGP4_SUN_ANOM_0      6.2565837
#define SGP4_SUN_ANOM_DOT    0.017201977

    // Fixed geometry. cos i of the Moon on the equator oscillates as
    // 0.91375164 - 0.03568096 cos(node), i.e. 18.30 to 28.59 deg;
    // 0.089683511 is sin of its 5.145 deg tilt to the ecliptic, and
    // the obliquity pair is sin/cos of 23.4441 deg.
#define SGP4_MOON_COSI_0    0.91375164
#define SGP4_MOON_COSI_AMP  0.03568096
#define SGP4_SIN_I_MOON     0.089683511
#define SGP4_SIN_OBLIQ      0.39785416
#define SGP4_COS_OBLIQ      0.91744867

    // Third-body disturbing amplitudes and mean motions. CC x ZN is the
    // amplitude actually used: 3.568e-11 for the Sun against 7.596e-11
    // for the Moon, a ratio of 2.13 against the 2.18 of GM/r^3. These
    // are fitted coefficients rather than derived ones, but the ratio
    // is the cheapest check that neither pair was mistyped.
#define SGP4_ZNS    1.19459e-5       // [rad/min]
#define SGP4_C1SS   2.9864797e-6
#define SGP4_ZES    0.01675          // Earth orbit eccentricity
#define SGP4_ZNL    1.5835218e-4     // [rad/min]
#define SGP4_C1L    4.7968065e-7
#define SGP4_ZEL    0.05490          // lunar orbit eccentricity

    // The Sun's argument of perigee: atan2 of the pair is 281.22 deg,
    // the longitude of perihelion.
#define SGP4_SUN_SING  (-0.98088458)
#define SGP4_SUN_COSG    0.1945905

    // Below three degrees the node rate is dropped: it divides by
    // sin i0, and a near-equatorial node carries nothing to divide.
#define SGP4_LOW_INCL   5.2359877e-2   // [rad] = 3.0000 deg

    // Below this the node and the argument of perigee are not
    // separately meaningful, and the lunisolar perturbation is applied
    // through Lyddane's non-singular pair instead. Tested against the
    // inclination at epoch, never the perturbed one: otherwise a
    // satellite can change theory halfway through a propagation.
#define SGP4_LYDDANE_INCL   0.2      // [rad] = 11.4592 deg

    // Resonance bands, on the recovered mean motion [rad/min]: the 24 h
    // window is 0.800 to 1.200 rev/day, the 12 h one 1.893 to 2.118
    // rev/day, and the latter additionally needs e >= 0.5.
#define SGP4_RES_24H_LO   0.0034906585
#define SGP4_RES_24H_HI   0.0052359877
#define SGP4_RES_12H_LO   8.26e-3
#define SGP4_RES_12H_HI   9.24e-3
#define SGP4_RES_12H_ECC  0.5

    // Eccentricity branch points of the 12 h G polynomials.
#define SGP4_G_ECC_1  0.65
#define SGP4_G_ECC_2  0.715
#define SGP4_G_ECC_3  0.7

    // Tesseral harmonic amplitudes. Q22 and ROOT22 carry the same
    // number under two names in the report; kept apart because they
    // enter two different resonances and nothing guarantees they stay
    // equal if either is ever refined.
#define SGP4_Q22     1.7891679e-6
#define SGP4_Q31     2.1460748e-6
#define SGP4_Q33     2.2123015e-7
#define SGP4_ROOT22  1.7891679e-6
#define SGP4_ROOT32  3.7393792e-7
#define SGP4_ROOT44  7.3636953e-9
#define SGP4_ROOT52  1.1428639e-7
#define SGP4_ROOT54  2.1765803e-9

    // Phases of the three synchronous tesseral terms [rad].
#define SGP4_FASX2  0.13130908
#define SGP4_FASX4  2.8843198
#define SGP4_FASX6  0.37448087

    // Earth rotation rate [rad/min] = 7.29211514668855e-5 rad/s, the
    // nominal rate; the report rounds it to eight digits as 4.3752691e-3
    // and that rounding was invisible for as long as the angle at epoch
    // came from the report too. It is not invisible now. Both this rate
    // and that angle enter sin(xli - fasx2), where they set the phase of
    // the resonance forcing, and a phase error there is integrated
    // rather than cancelled: eight digits cost a factor of four hundred
    // on a 24 h resonant element set. The eight-digit form survives as
    // the thing to check this against, not as the thing to compute with.
#define SGP4_THDT  4.37526908801129966e-3

    // Resonance integration: a twelve-hour Taylor step from epoch, and
    // the half-square that goes with it. STEP2 is written out rather
    // than computed because the report writes it out.
#define SGP4_RES_STEP   720.0
#define SGP4_RES_STEP2  259200.0      // 720^2 / 2

    // Phases of the resonant tesseral terms [rad].
#define SGP4_G22  5.7686396
#define SGP4_G32  0.95240898
#define SGP4_G44  1.8014998
#define SGP4_G52  1.0508330
#define SGP4_G54  4.4108898

#ifdef __cplusplus
}
#endif

#endif /* SPODY_CONST_H */