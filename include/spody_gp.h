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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */
#ifndef SPODY_GP_H
#define SPODY_GP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spody_sgp4.h"

/* One general-perturbation element set to one osculating ICRF state at
 * its own epoch: the bridge between a catalogue and an integrator.
 *
 * The elements arrive already decoded -- there is no card parser here
 * on purpose, because the caller has them as numbers, from a catalogue
 * or from whatever read the card. A decoder living here would be a
 * second decoder in the tree with one reader.
 *
 * What this does is short and the shortness is the point: run SGP4 to
 * tsince = 0, which is what turns mean elements into an osculating
 * state, then rotate TEME to ICRF. Both halves are already checked
 * against independent implementations; this is the wire between them.
 *
 * WHAT THE ANSWER IS WORTH. The rotation is exact to centimetres and
 * the propagator is bit-identical to the reference implementation, and
 * neither of those is the accuracy of the state. A GP element set is a
 * fit, worth kilometres, and this function cannot make it better --
 * it makes it *comparable*, which is a different service. Treating the
 * result as a precise osculating state because the digits are long is
 * the mistake this comment exists to prevent.
 *
 * Arguments:
 *   el           : mean elements. epoch_mjd is UTC, per the TLE
 *                  convention.
 *   eop_file     : path to IERS finals2000A.all.
 *   iau2006_dir  : path to the directory holding tab5.2{a,b,d}.txt.
 *   et_s         : out, ET (TDB seconds past J2000) of the epoch. The
 *                  state means nothing without it; it is returned
 *                  rather than left to the caller precisely because it
 *                  is the field a caller forgets.
 *   r_icrf_km    : out, position [km].
 *   v_icrf_kms   : out, velocity [km/s].
 *
 * Returns 0 on success, non-zero on a bad element set or unreadable
 * tables. */
int spody_gp_to_state_icrf(const spody_sgp4_elements *el,
                           const char *eop_file,
                           const char *iau2006_dir,
                           double *et_s,
                           double r_icrf_km[3],
                           double v_icrf_kms[3]);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_GP_H */
