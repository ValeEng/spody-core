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
 * CCSDS OEM (Orbit Ephemeris Message) -> SpOdy SPDYOUT_ state-vector
 * binary converter.
 *
 * Reads one or more CCSDS 502.0-B OEM text files (the format used by
 * e.g. the NASA/JSC public ISS trajectory ephemerides) and writes a
 * SpOdy SPDYOUT_ binary with full (r, v) records, suitable as a
 * reference trajectory for the diff-validation workflow and for the
 * `spody calibrate` density-scale fit.
 *
 * Supported subset (v1, rejected with a clear message otherwise):
 *   - REF_FRAME   : ICRF, EME2000 or J2000. Per the SPICE convention
 *                   J2000 == ICRF (the ~mas frame bias is ignored), so
 *                   all three map to the engine's ICRF axes with no
 *                   rotation and the converter needs no EOP data.
 *   - TIME_SYSTEM : UTC (leap-second chain + TDB periodic term) or
 *                   TDB (taken as ET directly).
 *   - Epochs in calendar ISO form YYYY-MM-DDThh:mm:ss[.sss]; the
 *     day-of-year OEM variant is not supported.
 *   - Ephemeris rows: epoch + x y z vx vy vz in km / km/s (the CCSDS
 *     mandatory units); optional trailing acceleration columns are
 *     ignored. COMMENT lines and covariance blocks are skipped.
 *
 * Output wire format (must match spody/src/sim_run.c):
 *   8 bytes  "SPDYOUT_"
 *   4 bytes  uint32 LE version = 1
 *   4 bytes  uint32 LE state_dim = 6
 *   8 bytes  reserved (two uint32 zero)
 *   then per record 7 LE doubles:
 *     (t [s, 0-anchored on the first record], x, y, z [km, ICRF],
 *      vx, vy, vz [km/s])
 *
 * Multi-file mode concatenates every input into one binary with a
 * continuous 0-anchored time axis (pass the files in chronological
 * order). Records whose epoch does not advance past the previously
 * written one are skipped -- consecutive daily OEM releases overlap,
 * and the earlier file's tail would otherwise duplicate epochs.
 *
 * Arguments:
 *   n_inputs        : number of input OEM paths (>= 1).
 *   input_oem_paths : the OEM text files, chronological order.
 *   output_bin      : path for the SPDYOUT_ binary to be written.
 *
 * Returns 0 on success, non-zero on I/O or parse failure.
 */
#ifndef SPODY_OEM_H
#define SPODY_OEM_H

#ifdef __cplusplus
extern "C" {
#endif

int spody_convert_oem_to_state_icrf(int n_inputs,
                                    const char *const *input_oem_paths,
                                    const char *output_bin);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_OEM_H */
