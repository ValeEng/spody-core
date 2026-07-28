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
 * Arguments -- all vectors in one common frame [km]; only relative
 * geometry is used, so which frame it is does not matter:
 *   sat2sun     satellite -> Sun
 *   sun_radius  km
 *   sat2occ     satellite -> centre of each occulting body, n_occ rows
 *   occ_radius  km, one per row; rows with radius <= 0 are ignored
 *   n_occ       number of rows, <= SPODY_ECL_MAX_OCCULTERS
 *
 * The occulter list is whatever the caller decides can cast a shadow;
 * it must never contain the Sun itself.
 */
double spody_get_satlitfraction(const double sat2sun[3], double sun_radius,
                                const double sat2occ[][3],
                                const double occ_radius[], int n_occ);

#ifdef __cplusplus
}
#endif

#endif // SPODY_ECLIPSE_H
