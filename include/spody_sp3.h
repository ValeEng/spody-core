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
 * IGS .sp3 -> SpOdy SPDYOUT_ trajectory binary converter.
 *
 * Reads a precise-orbit ephemeris file in the SP3 standard format (the
 * IGS distribution for GPS / GLONASS / Galileo / Beidou precise orbits;
 * see Hilla 2010, "The Extended Standard Product 3 Orbit Format (SP3-c)"
 * and the SP3-d update 2016) and writes a SpOdy SPDYOUT_ binary holding
 * one ICRF position record per SP3 epoch for the requested satellite.
 *
 * Purpose: produce reference trajectories that the val test rigs can
 * diff against a SpOdy propagation -- "did we predict GPS PRN 11 to
 * within X metres over 24 hours?" The val rig just reads two SPDYOUT_
 * binaries (one from us, one from IGS) and reports max/RMS distance.
 *
 * Format notes:
 *   - SP3 positions are in km in an Earth-fixed (ITRF) frame, typically
 *     ITRF2014 or ITRF2020 depending on the analysis centre. The
 *     differences between ITRF realisations are sub-cm and well below
 *     orbital-propagation accuracy targets, so we treat 'ITRF' as one
 *     frame for the rotation. The header line marks the realisation
 *     for the reader's information; we do not transform between
 *     realisations.
 *   - SP3 epochs are in GPS time. GPS time runs TAI - 19 s exactly
 *     (no leap seconds since the GPS epoch 1980-01-06), so TT = GPS +
 *     51.184 s and ET_TDB ~ TT (microsecond-level diff, negligible).
 *   - The output is a SpOdy SPDYOUT_ binary with state_dim=6: we keep
 *     the standard layout (t, x, y, z, vx, vy, vz) so the existing
 *     reader works without changes, but write vx=vy=vz=0 since SP3
 *     normally carries position only (the 'EV' velocity records are
 *     optional and ignored here). Diff tools must compare positions
 *     only against an SP3-derived reference.
 *
 * Lifecycle: this is a one-shot file conversion. The function opens
 * the EOP + IAU 2006 data internally, performs the conversion, frees
 * everything before returning. Not for hot paths -- intended for offline
 * benchmark preparation.
 */
#ifndef SPODY_SP3_H
#define SPODY_SP3_H

#ifdef __cplusplus
extern "C" {
#endif

/* Convert one satellite's track from *input_sp3* into *output_bin*.
 *
 * Returns 0 on success, non-zero on any IO / parse / write failure
 * (error message printed to stderr).
 *
 * Arguments:
 *   input_sp3    : path to the IGS .sp3 file (text).
 *   output_bin   : path for the SPDYOUT_ binary to be written.
 *   sat_id       : 3-character SP3 satellite identifier, e.g. "G11"
 *                  for GPS PRN 11, "E03" for Galileo E03, "R23" for
 *                  GLONASS slot 23. Case-sensitive, matched to the
 *                  first 3 characters of every "P<id> x y z ..." line.
 *   eop_file     : path to IERS finals2000A.all (xp, yp, dUT1, dX, dY).
 *   iau2006_dir  : path to the directory containing tab5.2{a,b,d}.txt
 *                  (IAU 2006 X / Y / s+XY/2 series).
 *
 * Side effects: informational lines printed to stderr (satellite id,
 * number of epochs read, time range, output record count).
 */
int spody_convert_sp3_to_state_icrf(const char *input_sp3,
                                    const char *output_bin,
                                    const char *sat_id,
                                    const char *eop_file,
                                    const char *iau2006_dir);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_SP3_H */
