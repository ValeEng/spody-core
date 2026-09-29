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
#ifndef SPODY_ECLIPSE_H
#define SPODY_ECLIPSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "spody_const.h"
#include "spody_math.h"      /* SpodyBodyShape */

//debug
#define DEBUG_ECLIPSE 0 // 0 = no debug | 1 = debug |---> CODE TESTING

#define BUFFER_SIZE_ECL 256

/* Upper bound on the number of occulting bodies accepted in one call.
 * Sized for "central body + every third body"; spody_forcemodels.h
 * checks its own SPODY_FM_MAX_THIRD against this at compile time. */
#define SPODY_ECL_MAX_OCCULTERS 9

/*
 * Fraction of the Sun's light reaching the satellite: 1 in full
 * sunlight, 0 in total eclipse, in between in penumbra. Single entry
 * point of the eclipse machinery -- both the SRP force and the
 * eclipse event go through it.
 *
 * The geometry is Montenbruck & Gill's: the Sun and every occulting
 * body are projected onto the sky plane as discs with the angular
 * radius they subtend at the satellite, and the shadow problem
 * becomes planar circle overlap.
 *
 * Several bodies can occult at once -- the Moon crossing the Sun
 * while the satellite is entering the Earth's penumbra, the Earth
 * seen from a lunar orbit. The hidden part of the solar disc is their
 * UNION, so the overlaps must not be counted twice; the sum is
 * inclusion-exclusion,
 *
 *     hidden = sum_i g_i - sum_i<j g_ij + ...
 *
 * with g_i the fraction of the solar disc covered by body i and g_ij
 * the fraction covered by bodies i and j at the same time. Terms of
 * order three and higher are dropped (they need three bodies over the
 * same piece of the Sun simultaneously) and the result is clamped
 * into the bracket that holds whatever those terms would have been:
 *
 *     1 - sum_i g_i  <=  lit  <=  min_i (1 - g_i)
 *
 * Occulting bodies are spheres or oblate spheroids (SpodyBodyShape).
 * A spheroid enters the disc model with the angular radius of its
 * limb in the plane (satellite, body centre, Sun): the affine map that
 * stretches the polar axis by r_eq / r_pol turns it into a sphere and
 * keeps tangency, so the limb point comes in closed form. The
 * spheroid lies between its polar and its equatorial sphere, so the
 * limb is only computed where those two disagree, the seconds around
 * each contact. The same plane-of-the-Sun construction as Adhya,
 * Sibthorpe, Ziebart & Cross (J. Spacecraft Rockets 41(1), 2004),
 * extended here from a lit/shadow state to the lit fraction.
 *
 * Arguments -- all vectors in one common frame [km]; only relative
 * geometry is used, so which frame it is does not matter as long as
 * the poles in `occ` are in it too:
 *   sat2sun     satellite -> Sun
 *   sun_radius  km
 *   sat2occ     satellite -> centre of each occulting body, n_occ rows
 *   occ         shape of each occulter, one per row; rows with
 *               r_eq <= 0 are ignored
 *   n_occ       number of rows, <= SPODY_ECL_MAX_OCCULTERS
 *
 * The occulter list is whatever the caller decides can cast a shadow;
 * it must never contain the Sun itself.
 */
double spody_get_satlitfraction(const double sat2sun[3], double sun_radius,
                                const double sat2occ[][3],
                                const SpodyBodyShape occ[], int n_occ);

/*
 * Signed predicate of an eclipse event against ONE occulting body:
 * positive when the satellite is more lit than `threshold`, negative
 * when less, zero on the crossing. What a root finder needs, which
 * the lit fraction alone is not at the two ends of its range: it is
 * exactly 0 all through the umbra and exactly 1 all through sunlight,
 * so "fraction - 0" and "fraction - 1" are flat on one side of their
 * root (a bracketing solver then stops on any point of the flat side,
 * and "fraction - 1" never even changes sign). With a, b, c the Sun's
 * and the body's angular radii and their separation (as in
 * spody_get_satlitfraction):
 *
 *     threshold <= 0   c - (b - a)          umbra contact (second/third)
 *     threshold >= 1   c - (a + b)          penumbra contact (first/fourth)
 *     otherwise        lit fraction - threshold
 *
 * The contact forms are the boundaries of the disc-overlap cases that
 * spody_get_satlitfraction already branches on (Montenbruck & Gill's
 * conical shadow: full sunlight for c >= a + b, total eclipse for
 * c <= b - a), written as signed angles, linear through the root. In
 * between, the fraction crosses the threshold inside the penumbra,
 * where its slope is finite. A satellite inside the body returns -1.
 * For a spheroid, b is its limb radius, computed at every call (not
 * only near the contacts as in the force) so the residual stays
 * continuous.
 *
 * Unlike spody_get_satlitfraction there is no sunward-side screening:
 * there c is near 90 degrees or more and every form is positive.
 */
double spody_get_eclipse_residual(const double sat2sun[3], double sun_radius,
                                  const double sat2occ[3],
                                  const SpodyBodyShape *occ,
                                  double threshold);

#ifdef __cplusplus
}
#endif

#endif // SPODY_ECLIPSE_H
