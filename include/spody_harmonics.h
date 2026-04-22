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
// Normalized Derived Associated Lengendre Polynomials (of the 1st kind)
// per the method of Ref.[1]. "Fully" normalized for compatibility with
// the EGM96 coefficients per Ref[2]
//
// REFERENCES:
// [1] Lundberg, J.B., and Schutz, B.E., "Recursion Formulas of Legendre
//     Functions for Use with Nonsingular Geopotential Models", Journal
//     of Guidance, Dynamics, and Control, Vol. 11, No.1, Jan.-Feb. 1988.
//
// [2] Heiskanen, W.A., and Moritz, H., "Physical Geodesy", W.H. Freeman
//     and Company, San Francisco, 1967.
//
// [3] Pines, S., "Uniform Representation of the Gravitational Potential
//     and its Derivatives", AIAA Journal, Vol. 11, No. 11, 1973.

#ifndef SPODY_HARMONICS_H
#define SPODY_HARMONICS_H

#ifdef __cplusplus
extern "C" {
#endif

#define DEBUG_HARMONICS 0
#define BUFFER_SIZE_HG 256

typedef struct {
    int N;
    double R_ref;
    double GM;
    double *C; // normalized coefficients (GRAIL) 
    double *S; // normalized coefficients (GRAIL)
    double *recurr_a; // common recurrence factors
    double *recurr_b; // common recurrence factors
} HarmonicGravityData;

typedef struct {
    HarmonicGravityData *hgd;
    double *A_row0; // Row n+1 (future)
    double *A_row1; // Row n (current)
    double *A_row2; // Row n-1 (previus)
    double *real;
    double *imag;
} HarmonicGravity;

void spody_get_hgaccbodyfixed(HarmonicGravity *hg, double pos[3], double acc_out[3]);
int spody_load_HarmonicGravityData(HarmonicGravityData *hgd, const char *filename, int degree);
int spody_setup_HarmonicGravity(HarmonicGravity *hg, HarmonicGravityData *hgd);
int spody_free_HarmonicGravity(HarmonicGravity *hg);
int spody_free_HarmonicGravityData(HarmonicGravityData *hgd);

#ifdef __cplusplus
}
#endif

#endif // SPODY_HARMONICS_H