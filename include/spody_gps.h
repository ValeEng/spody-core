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
 * RINEX-NAV (GPS) -> SpOdy SPDYOUT_ state-vector binary converter.
 *
 * GPS broadcast nav messages carry a *Kepler-with-corrections* set
 * of 16 orbital parameters (a, e, i0, Omega0, omega, M0, plus six
 * second-order harmonic correction coefficients Cuc/Cus/Crc/Crs/Cic/Cis
 * and three time-rate terms Delta_n, OmegaDot, iDot) per record,
 * valid for the ~4 h fit interval after Toe (time of ephemeris). The
 * propagation algorithm is the canonical one in IS-GPS-200 sect.
 * 20.3.3.4.3 -- solve Kepler, apply corrections, transform from the
 * orbital plane to ECEF (WGS-84). For velocity we add the analytic
 * derivatives (Remondi 2004, "Computing satellite velocity using
 * the broadcast ephemeris").
 *
 * The result is in WGS-84 ECEF; we rotate to the central-body
 * inertial frame the same way spody_glonass does:
 *
 *   r_ICRF = R_bf2i(et) . r_ECEF
 *   v_ICRF = R_bf2i(et) . v_ECEF  +  omega_earth x r_ICRF
 *
 * where omega_earth uses the full ICRF axis (third column of R_bf2i)
 * rather than nominal (0, 0, omega) -- at GPS altitude the nominal
 * z-hat approximation misses ~few m/s of velocity due to the J2024
 * precession-driven ITRS z-axis tilt.
 *
 * Time scale. GPS broadcast TOC is GPS time per RINEX 3.05 sect.
 * 6.10.1: TT = GPST + 51.184 s exactly (= GPS->TAI of 19 s plus
 * TAI->TT of 32.184 s). No leap-second table required, unlike
 * GLONASS UTC. TT -> TDB adds the deltet periodic term
 * (spody_tdb_minus_tt, +/-1.657 ms).
 *
 * Wire format we emit (must match spody/src/sim_run.c, spody_sp3.c
 * and spody_glonass.c):
 *   8 bytes  "SPDYOUT_"
 *   4 bytes  uint32 LE version=1
 *   4 bytes  uint32 LE state_dim=6
 *   8 bytes  reserved
 *   then per record: 7 little-endian doubles
 *     (t [s past 1st record], x, y, z [km, ICRF], vx, vy, vz [km/s, ICRF])
 *
 * Like the glonass converter (and unlike the sp3 converter), the GPS
 * broadcast binary carries real velocity in the v slots. Diff tools
 * can compare both position and velocity against a propagation.
 *
 * Multi-file: IGS-BKG ships GPS broadcast as one RINEX-NAV file per
 * UTC day. For week-or-more-long validation, pass the files in
 * chronological order and the converter concatenates them into one
 * SPDYOUT_ binary with a continuous 0-anchored time axis.
 */
#ifndef SPODY_GPS_H
#define SPODY_GPS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Convert one GPS satellite's broadcast track from one or more
 * RINEX inputs into a single concatenated *output_bin*. Returns 0
 * on success, non-zero on any IO / parse / write failure (error
 * message printed to stderr).
 *
 * Arguments:
 *   n_inputs          : number of RINEX paths (>= 1).
 *   input_rnx_paths   : array of length n_inputs, paths to RINEX 3.x
 *                       GPS-or-mixed nav files in chronological order.
 *   output_bin        : path for the SPDYOUT_ binary to be written.
 *   sat_id            : 3-char GPS PRN, e.g. "G11". Case-sensitive.
 *   eop_file          : path to IERS finals2000A.all.
 *   iau2006_dir       : path to the directory containing tab5.2{a,b,d}.txt.
 *
 * Side effects: per-file and aggregate informational summaries printed
 * to stderr (sat id, record count, time span).
 */
int spody_convert_gps_to_state_icrf(int n_inputs,
                                    const char *const *input_rnx_paths,
                                    const char *output_bin,
                                    const char *sat_id,
                                    const char *eop_file,
                                    const char *iau2006_dir);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_GPS_H */
