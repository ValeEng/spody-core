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
 * RINEX-NAV (GLONASS) -> SpOdy SPDYOUT_ state-vector binary converter.
 *
 * GLONASS is unique among GNSS constellations: its broadcast nav
 * messages carry the satellite state vector (position, velocity, and
 * luni-solar acceleration) directly, in PZ-90 every 30 min. GPS /
 * Galileo / BeiDou instead broadcast Keplerian elements, requiring an
 * extra Keplerian-propagation step to recover (r, v).
 *
 * This converter reads a RINEX 3.x navigation file, extracts every
 * (TOC, r, v) record for the requested GLONASS satellite, rotates
 * each into ICRF, and writes a SpOdy SPDYOUT_ binary suitable for
 * direct diff against a propagator output -- without the orbit-
 * determination workaround that GPS-only-from-sp3 would require.
 *
 * Frame transformation. The broadcast (r, v) is in PZ-90, which is
 * the Russian terrestrial frame (an ITRF realisation; PZ-90.11 differs
 * from ITRF2014 by ~few cm in position and arc-seconds in rotation
 * -- negligible for orbital propagation, so we treat PZ-90 ~ ITRF):
 *
 *   r_ICRF = R_bf2i(et) . r_ECEF
 *   v_ICRF = R_bf2i(et) . v_ECEF  +  omega_earth x r_ICRF
 *
 * The omega_earth term is mandatory: in the rotating PZ-90 frame the
 * Earth-rotation contribution to velocity is removed, so we add it
 * back to get inertial velocity. omega_earth ~ 7.2921151467e-5 rad/s
 * along the celestial pole (we use the nominal IERS value; polar-
 * motion corrections are <1e-7 in magnitude and well below broadcast
 * precision).
 *
 * Time scale. GLONASS broadcast TOC is UTC per RINEX 3 spec sect. 6.
 * We chain UTC -> TAI -> TT -> TDB using the leap-second table and
 * the deltet periodic term in spody_time.c (exact at any post-1972
 * epoch) plus TT - TAI = 32.184 s.
 *
 * Wire format we emit (must match spody/src/sim_run.c and the sp3
 * converter):
 *   8 bytes  "SPDYOUT_"
 *   4 bytes  uint32 LE version=1
 *   4 bytes  uint32 LE state_dim=6
 *   8 bytes  reserved
 *   then per record: 7 little-endian doubles
 *     (t [s past J2000 TDB], x, y, z [km, ICRF], vx, vy, vz [km/s, ICRF])
 *
 * Unlike the sp3 converter (which writes vx=vy=vz=0), the GLONASS
 * binary carries real velocity in the v slots. Diff tools can compare
 * both position and velocity against a propagation.
 */
#ifndef SPODY_GLONASS_H
#define SPODY_GLONASS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Convert one satellite's broadcast track from one or more RINEX
 * input files into a single concatenated *output_bin*. Returns 0 on
 * success, non-zero on any IO / parse / write failure (error message
 * printed to stderr).
 *
 * Multi-file mode is the workhorse for week-long (or longer) GNSS
 * validations: IGS-BKG ships GLONASS broadcast nav as one RINEX file
 * per UTC day, so a 7-day window is 7 daily files passed in
 * chronological order. The converter concatenates them into one
 * SPDYOUT_ binary with a single header and a continuous, 0-anchored
 * time axis (every record's t = et - et_of_first_record_overall, no
 * gaps at day boundaries). Calling with n_inputs==1 reproduces the
 * single-file behaviour bit-for-bit.
 *
 * Arguments:
 *   n_inputs          : number of RINEX paths (>= 1).
 *   input_rnx_paths   : array of length n_inputs, paths to RINEX 3.x
 *                       nav files in chronological order.
 *   output_bin        : path for the SPDYOUT_ binary to be written.
 *   sat_id            : 3-char GLONASS slot id, e.g. "R03". Case-sensitive.
 *   eop_file          : path to IERS finals2000A.all.
 *   iau2006_dir       : path to the directory containing tab5.2{a,b,d}.txt.
 *
 * Side effects: per-file and aggregate informational summaries printed
 * to stderr (sat id, record count, time span).
 */
int spody_convert_glonass_to_state_icrf(int n_inputs,
                                        const char *const *input_rnx_paths,
                                        const char *output_bin,
                                        const char *sat_id,
                                        const char *eop_file,
                                        const char *iau2006_dir);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_GLONASS_H */
