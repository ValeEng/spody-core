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
 * ICGEM .gfc -> SpOdy GRGM-style .tab converter.
 *
 * Reads a spherical-harmonic gravity-field file in the ICGEM .gfc text
 * format (the de-facto standard for modern Earth/lunar/Mars models:
 * EGM2008, EGM2020, EIGEN-6C4, GOCE-DIR, GRACE GGM05C, ...) and writes
 * a .tab file in the exact CSV layout consumed by
 * spody_load_HarmonicGravityData (the loader currently used for the
 * GRGM1200B lunar model). Adding a new model is then a one-liner:
 * convert once, drop the .tab next to the others, point the TOML at it.
 *
 * Format references:
 *   - ICGEM .gfc: F. Barthelmes & C. Foerste, "The ICGEM-format",
 *     GFZ Potsdam, doc. v1.0, 2011 (kept current on the ICGEM portal
 *     http://icgem.gfz-potsdam.de/).
 *   - EGM2008: N.K. Pavlis, S.A. Holmes, S.C. Kenyon, J.K. Factor,
 *     "The development and evaluation of the Earth Gravitational Model
 *     2008 (EGM2008)", J. Geophys. Res. 117, B04406, 2012.
 *   - SpOdy .tab header line: see
 *     spody-core/src/spody_harmonics.c::read_spherical_harmonics_file.
 *
 * Lifecycle: this is a one-shot file conversion (no live state to keep).
 * The function allocates ~75 MB of working memory for a full N=2190
 * Earth model and frees it before returning. Designed to be called from
 * the SpOdy app's `convert harmonics_icgem` subcommand.
 */
#ifndef SPODY_ICGEM_H
#define SPODY_ICGEM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Convert *input_gfc* to *output_tab*. Returns 0 on success, non-zero on
 * any IO / parse / write failure (error message printed to stderr).
 *
 * max_degree:
 *   > 0  truncate the written file to this degree (must be <= the
 *        intrinsic max_degree declared in the .gfc header).
 *   <= 0 keep the full model. Recommended: the SpOdy harmonics loader
 *        truncates again at read time according to the user's TOML, so
 *        storing the full coefficient set keeps every future run free
 *        to pick its own truncation without a re-conversion.
 *
 * Side effects: informational lines printed to stderr (model name, GM,
 * radius, declared max_degree, tide system, count of skipped time-
 * variable rows) so the operator can confirm the file's identity.
 */
int spody_convert_icgem_to_tab(const char *input_gfc,
                               const char *output_tab,
                               int max_degree);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_ICGEM_H */
