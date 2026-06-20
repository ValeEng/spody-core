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

/* Per-thread handle. Today this is a thin wrapper around the shared
 * data pointer (the evaluators are stateless); kept as a separate
 * type for symmetry with MappedEphemeris / MappedEOP, so adding
 * per-thread caches in the future does not change the API. */
typedef struct MappedIAU2006 {
    const MappedIAU2006Data *data;
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

/* Earth Rotation Angle, Capitaine's strictly-linear function of UT1
 * (TN 36 eq. 5.15):
 *
 *   ERA(Tu) = 2*pi * (0.7790572732640 + 1.00273781191135448 * Tu)
 *
 * where Tu = JD(UT1) - 2451545.0 (days, fractional).
 *
 * Returns ERA in radians, reduced to [0, 2*pi). */
double spody_iau2006_era(double jd_ut1);

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
void spody_bf_rotation_earth(const struct ForceModelContext *ctx, double et,
                              double R_icrf_to_bf[3][3],
                              double R_bf_to_icrf[3][3]);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_EARTH_ORIENTATION_H */
