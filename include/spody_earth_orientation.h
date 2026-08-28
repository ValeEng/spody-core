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
 * Earth orientation: ICRF (GCRS) <-> ITRF rotation via the IAU 2006
 * precession + IAU 2000A_R06 nutation model, CIO-based formulation
 * of Capitaine, with IERS EOP corrections (xp, yp, dUT1, dX, dY).
 *
 * REFERENCES
 * ----------
 *   IERS Conventions (2010), IERS Technical Note 36
 *     Gerard Petit and Brian Luzum, eds. ISBN 3-89888-989-6
 *     <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html>
 *   Chapter 5 (Transformation between the ITRS and the GCRS), eqs.
 *   (5.10), (5.11), (5.15) for ERA, (5.13) for polar motion.
 *
 *   Coefficient tables tab5.2a, tab5.2b, tab5.2d distributed by
 *   IERS Conventions Centre:
 *     <https://iers-conventions.obspm.fr/content/chapter5/additional_info/>
 *   Tab 5.2a: 1306+253+36+4+1 terms for the X coordinate of the CIP
 *   Tab 5.2b: 962+277+30+5+1 terms for the Y coordinate of the CIP
 *   Tab 5.2d: 33+3+25+4+1 terms for the s(t)+XY/2 quantity
 *
 *   Capitaine, N. and Wallace, P.T. (2006), "High precision methods
 *   for locating the celestial intermediate pole and origin",
 *   A&A 450, 855-872.
 *
 *   Mathews, P.M., Herring, T.A., Buffett, B.A. (2002), "Modeling of
 *   nutation and precession: New nutation series for nonrigid Earth
 *   and insights into the Earth's interior", JGR 107(B4), 2068.
 *   (The MHB2000 nutation model underlying IAU 2000A.)
 *
 * MODEL SUMMARY
 * -------------
 *     [r_ITRF] = W(t) . R3(-ERA(Tu)) . Q(t) . [r_GCRS]
 *
 *   Q(t)         GCRS -> Celestial Intermediate Reference System (CIRS),
 *                built from the CIP coordinates X, Y and the CIO
 *                locator s (eqs. 5.10, 5.11). The X, Y, s series are
 *                IAU 2006 precession + IAU 2000A_R06 nutation
 *                evaluated at TT.
 *
 *   R3(-ERA)     CIRS -> Terrestrial Intermediate Reference System
 *                (TIRS). The Earth Rotation Angle ERA is Capitaine's
 *                strictly-linear function of UT1 (eq. 5.15) -- no
 *                equation-of-equinoxes drift, unlike the legacy GMST.
 *
 *   W(t)         TIRS -> International Terrestrial Reference Frame
 *                (ITRF). Encodes polar motion (xp, yp from IERS EOP)
 *                and the negligible TIO locator s'(t) (eq. 5.13).
 *
 * IERS EOP corrections applied to the CIP coordinates:
 *     X_actual = X_model + dX_iers
 *     Y_actual = Y_model + dY_iers
 *   with dX, dY interpolated from finals2000A.all at the run epoch.
 *   Polar motion xp, yp are taken directly from the same table.
 *   UT1 - UTC (dUT1) drives the ERA, also from finals2000A.all.
 */
#ifndef SPODY_EARTH_ORIENTATION_H
#define SPODY_EARTH_ORIENTATION_H

#include <stddef.h>
#include <stdint.h>

#include "spody_eop.h"   /* MappedEOP (consumed by the rotation function) */

#ifdef __cplusplus
extern "C" {
#endif

/* One term of an IAU 2006 trigonometric series. The amplitudes are
 * in microarcseconds (uas) as published in the IERS tables; we
 * convert to radians at the end of the series evaluation.
 *
 * `arg_coef[k]` is the integer multiplier of the k-th fundamental
 * argument in the trigonometric argument:
 *
 *   ARG = sum_{k=0}^{13} arg_coef[k] * F_k(t)
 *
 * with the F_k in the conventional IERS order (eqs. 5.43, 5.44):
 *
 *   k=0..4  : Delaunay arguments (luni-solar): l, l', F, D, Omega
 *   k=5..12 : planetary mean longitudes: L_Me, L_Ve, L_E, L_Ma,
 *             L_J, L_Sa, L_U, L_Ne
 *   k=13    : general precession p_A
 *
 * int8_t is wide enough for the IAU 2006 tables (max |coef| <= 9).
 */
#define IAU2006_NARG 14

typedef struct {
    double a_s_uas;                       /* sine amplitude, uas    */
    double a_c_uas;                       /* cosine amplitude, uas  */
    int8_t arg_coef[IAU2006_NARG];        /* fundamental-arg multipliers */
} IAU2006Term;

/* One IAU 2006 series (X, Y, or s+XY/2). Polynomial part is
 * 6 coefficients (constant to t^5), units of microarcseconds.
 * The non-polynomial part is a flat array of terms grouped by
 * the power of t they carry (j=0..4). j_offsets[j] is the
 * starting index of section j and j_offsets[5] == n_terms is the
 * sentinel for the last section. */
typedef struct {
    double       poly[6];           /* p_0 + p_1*t + ... + p_5*t^5,  uas */
    IAU2006Term *terms;             /* heap, sorted by j section          */
    size_t       n_terms;
    size_t       j_offsets[6];      /* section delimiters, j_offsets[5]=n */
} IAU2006Series;

/* Shared, read-only IAU 2006 series triplet. Parsed once from the
 * IERS-distributed text tables in `tab5.2a.txt`, `tab5.2b.txt`,
 * `tab5.2d.txt`. Thread-safe to share: every evaluation is read-only
 * over the term arrays. */
typedef struct MappedIAU2006Data {
    IAU2006Series X;     /* CIP X coordinate                  (tab5.2a) */
    IAU2006Series Y;     /* CIP Y coordinate                  (tab5.2b) */
    IAU2006Series s_xy;  /* CIO locator quantity s(t)+XY/2    (tab5.2d) */
} MappedIAU2006Data;

/* Node spacing of the (X, Y, s) interpolation grid, in seconds.
 *
 * The series costs ~70 us to evaluate -- 3084 terms, each a 14-term
 * linear combination plus a sine and a cosine -- and the Earth force
 * model needs it at every RHS evaluation, 78k times over a week-long
 * GNSS propagation. X and Y are smooth, moving about 650 and 340 mas
 * respectively over eight days, so a cubic through hourly nodes
 * reproduces them to 3e-7 mas: roughly 40 nanometres at GNSS radius,
 * six orders of magnitude below the millimetre these runs care about
 * and two orders below the uncertainty of the EOP inputs themselves.
 *
 * ERA and polar motion are deliberately NOT interpolated. ERA advances
 * 15 arcsec per second, so interpolating it would be the one thing
 * here that actually costs accuracy, and both are cheap closed forms
 * anyway. */
#define SPODY_XYS_NODE_S 3600.0

/* Points in the interpolation stencil. Four gives a cubic. */
#define SPODY_XYS_STENCIL 4

/* Per-thread handle: the shared data pointer plus this thread's
 * interpolation cache. One handle per thread, as with
 * MappedEphemeris / MappedEOP. */
typedef struct MappedIAU2006 {
    const MappedIAU2006Data *data;

    /* Cache for spody_iau2006_xys_interp.
     *
     * The nodes sit on a FIXED grid anchored at J2000, not on a window
     * around the current time. That is what keeps the result
     * reproducible: the interpolated X(t) is then a function of t
     * alone, identical no matter which step sequence the integrator
     * took to reach it. A window that followed the current time would
     * make the force depend on the path, and two runs of the same case
     * could legitimately differ. */
    int    cache_valid;
    long   cache_base;          /* grid index of the first stencil node */
    double node_x[SPODY_XYS_STENCIL];
    double node_y[SPODY_XYS_STENCIL];
    double node_s[SPODY_XYS_STENCIL];
} MappedIAU2006;

/* Load the three IAU 2006 series from a directory containing the
 * IERS tables `tab5.2a.txt`, `tab5.2b.txt`, `tab5.2d.txt`
 * (exact filenames). Heap-allocates the term arrays; the caller
 * frees with `spody_free_MappedIAU2006Data`.
 *
 * Returns 0 on success, -1 on parse / IO failure. */
int spody_setup_MappedIAU2006Data(MappedIAU2006Data *data, const char *dir);

/* Release the heap arrays. Safe on a zero-initialised struct. */
int spody_free_MappedIAU2006Data(MappedIAU2006Data *data);

/* Bind a per-thread handle to the shared data. No allocation. */
int spody_setup_MappedIAU2006(MappedIAU2006 *map, const MappedIAU2006Data *data);

/* Release the per-thread handle. Today a no-op. */
int spody_free_MappedIAU2006(MappedIAU2006 *map);

/* Evaluate the IAU 2006 (X, Y, s) at TT Julian centuries past J2000.
 *
 *   t = (jd_TT - 2451545.0) / 36525.0
 *
 * Outputs in radians. X and Y are the CIP coordinates in the GCRS;
 * s is the CIO locator after subtracting the XY/2 quantity stored
 * in the table (i.e. s = s_xy_table - X*Y/2 per IERS convention).
 *
 * Returns 0 on success, -1 on a NULL map. */
int spody_iau2006_xys(const MappedIAU2006 *map, double t_tt_cy,
                     double *X_rad, double *Y_rad, double *s_rad);

/* Same quantities, cubic-interpolated on the fixed hourly grid.
 *
 * This is what the force model calls. It exists because the exact
 * evaluation dominated the cost of an Earth-centred propagation --
 * ~73 us of the ~85 us spent per RHS evaluation at harmonic degree 70,
 * so the price of the gravity field was almost entirely the price of
 * re-deriving the Earth's orientation for it.
 *
 * Accuracy relative to `spody_iau2006_xys`: 3e-7 mas, i.e. tens of
 * nanometres of orbit error at GNSS radius. Use the exact function
 * when validating the series itself; use this one to propagate.
 *
 * `map` is mutated (the node cache), so one handle per thread.
 *
 * Returns 0 on success, -1 on a NULL map. */
int spody_iau2006_xys_interp(MappedIAU2006 *map, double t_tt_cy,
                             double *X_rad, double *Y_rad, double *s_rad);

/* Earth Rotation Angle, Capitaine's strictly-linear function of UT1
 * (TN 36 eq. 5.15):
 *
 *   ERA(Tu) = 2*pi * (0.7790572732640 + 1.00273781191135448 * Tu)
 *
 * where Tu = JD(UT1) - 2451545.0 (days, fractional).
 *
 * Returns ERA in radians, reduced to [0, 2*pi). */
double spody_iau2006_era(double jd_ut1);

/* Greenwich Mean Sidereal Time, IAU 1982 (Astronomical Almanac):
 *
 *   GMST[s] = 67310.54841 + (876600*3600 + 8640184.812866) * T
 *           + 0.093104 * T^2 - 6.2e-6 * T^3
 *
 * with T = (JD(UT1) - JD_J2000) / 36525 and 240 seconds of time to the
 * degree. Same physical angle as spody_iau2006_era, measured from the
 * equinox instead of the CIO.
 *
 * The caller decides what "UT1" means, and the two callers here mean
 * different things on purpose. A frame conversion passes a true UT1
 * built from the EOP table. SGP4 passes the TLE epoch unadjusted,
 * because that is what the element set was fitted against and what the
 * published verification vectors reproduce: feeding it dUT1 would be
 * more correct as time and less correct as an answer. This is not an
 * inconsistency to tidy away.
 *
 * Returns GMST in radians, reduced to [0, 2*pi). */
double spody_gmst1982(double jd_ut1);

/* Build the polar-motion matrix W(t):
 *
 *   W = R3(-s'(t)) . R2(xp) . R1(yp)
 *
 * with TIO locator s' = -47e-6 * t arcsec (eq. 5.13). All angle
 * arguments are in radians; t_tt_cy is TT Julian centuries past
 * J2000 (used only for the s' linear term, which is below mas
 * level over a century). */
void spody_iau2006_polar_motion(double t_tt_cy, double xp_rad, double yp_rad,
                                  double W[3][3]);

/* Forward declaration of the typedef-name `ForceModelContext`.
 * Defined in full in spody_forcemodels.h. Declared here so
 * spody_bf_rotation_earth's signature uses the SAME tag-name resolution
 * the .c file does -- a `struct ForceModelContext *` parameter declared
 * before the typedef is in scope would otherwise become a fresh
 * function-local struct type under strict C, and gcc / clang flag the
 * resulting prototype-vs-definition mismatch as a hard error. MSVC is
 * lenient and resolves the name to the global struct, hiding the bug
 * on Windows. The forward typedef below makes the two compilers agree.
 */
struct ForceModelContext;
typedef struct ForceModelContext ForceModelContext;

/* Top-level rotation: ICRF (GCRS) <-> ITRF at the given ET.
 *
 * Mirrors `spody_bf_rotation_moon` in shape (matches the
 * `spody_bf_rotation_fn` typedef in spody_forcemodels.h) so the
 * spherical-harmonics force can call it generically.
 *
 * Pulls runtime inputs from the ForceModelContext:
 *   ctx->eop      : interpolated for xp, yp, dUT1, dX, dY at `et`
 *   ctx->iau2006  : X, Y, s series evaluated at TT
 *
 * Both pointers MUST be non-NULL when this callback is registered
 * as ctx->get_bf_rotation. The Moon callback ignores them, so a
 * lunar run can leave them unset (NULL). */
void spody_bf_rotation_earth(const ForceModelContext *ctx, double et,
                              double R_icrf_to_bf[3][3],
                              double R_bf_to_icrf[3][3]);

/* Rotation between TEME and ICRF (GCRS) at the given ET.
 *
 * TEME has no official definition, and the ambiguity is not academic.
 * Building it on the IAU 2006 mean equinox instead of the 1982 one
 * moves the frame by up to 67 mas between 1981 and 2025 -- 13.7 m at
 * geostationary radius. This function takes the 1982 equinox, because
 * that is the one SGP4 itself turns on: its internal sidereal angle is
 * GMST 1982, checked to 5e-09 km against 634 published reference
 * states. A TEME built on the 2006 equinox would be a frame the
 * propagator feeding it does not share.
 *
 * Cross-checked two ways. Against astropy, which composes the frame as
 * this one does, agreement is 4.3 uas over 44 years. Against the
 * classical route -- bias-precession-nutation and the equation of the
 * equinoxes, which uses no UT1 at all and so would expose any failure
 * of GMST and ERA to cancel -- agreement is 6.7 uas rms once the
 * equinox convention above is accounted for.
 * And end to end against Skyfield, which carries its own precession
 * and nutation rather than erfa: propagating a real element set and
 * printing TEME and ICRF side by side, the two disagreements match to
 * two centimetres, so whatever moves the state is the propagator and
 * not this rotation.

 * TEME -- true equator, mean equinox of date -- is the frame SGP4 is
 * defined in, and the only frame a GP element set can be propagated
 * into. It is quasi-inertial: it does not turn with the Earth, so the
 * same matrix rotates position and velocity and there is no
 * omega x r term. What little it does turn is precession, some
 * 7.7e-12 rad/s.
 *
 * The chain collapses. Writing the standard one as
 *
 *   ITRS <- GCRS = W . R3(ERA) . Q^T          (SOFA iauC2t06a)
 *
 * and the classical entry as PEF <- TEME = R3(GMST) -- GMST and not
 * GAST, because TEME's origin is the MEAN equinox -- and taking the
 * pseudo-earth-fixed frame for the terrestrial intermediate one:
 *
 *   GCRS <- TEME = Q . R3(-ERA) . W^T . W . R3(GMST)
 *                = Q . R3(GMST - ERA)
 *
 * Polar motion cancels outright: both routes pass through the same
 * intermediate frame, so W meets its own transpose. xp and yp are
 * still read from the EOP table here, and still unused, because the
 * dUT1 that comes with them is not.
 *
 * What PEF and TIRS do not share is the TIO locator s', -47 uas per
 * century, so the angle actually carried is GMST - ERA - s'. Checked
 * against astropy on eight epochs between 2019 and 2025: the residual
 * about the pole axis is 9.1 uas rms without that term and 2.7 with
 * it, and the 2.7 that remain are dUT1 interpolated differently on the
 * two sides. The residual in the other two axes is larger and is not
 * ours -- it tracks the dX/dY celestial-pole offsets, which this
 * function applies and astropy does not.
 *
 * Same contract as spody_bf_rotation_earth: ctx->eop and ctx->iau2006
 * must both be non-NULL, and both matrices come back as the identity
 * if they are not or if the EOP interpolation fails. */

 
void spody_teme2icrf_rotation(const ForceModelContext *ctx, double et,
                              double R_teme_to_icrf[3][3],
                              double R_icrf_to_teme[3][3]);


#ifdef __cplusplus
}
#endif

#endif /* SPODY_EARTH_ORIENTATION_H */
