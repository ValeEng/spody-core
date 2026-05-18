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

    //debug
#define DEBUG_HARMONICS 0 // 0 = no debug | 1 = debug |---> CODE TESTING

#define BUFFER_SIZE_HG 256

typedef struct {
    double *C; // normalized coefficients (GRAIL)
    double *S; // normalized coefficients (GRAIL)
    double *recurr_a; // common recurrence factors
    double *recurr_b; // common recurrence factors
    double *recurr_s3; // sqrt((n-m)(n+m+1)(2-delta)/2)
    double *recurr_s4; // sqrt((n+m+2)(n+m+1)(2n+1)(2-delta) / (2(2n+3)))
    double R_ref;
    double GM;
    int N;
} HarmonicGravityData;

typedef struct {
    const HarmonicGravityData *hgd;
    double *A_row0; // Row n+1 (future)
    double *A_row1; // Row n (current)
    double *A_row2; // Row n-1 (previus)
    double *real;
    double *imag;
} HarmonicGravity;

void spody_get_hgaccbodyfixed(HarmonicGravity *hg, double pos[3], double acc_out[3]);

/* High-performance variant of spody_get_hgaccbodyfixed.
 * Same numerical contract as the reference but optimized for throughput:
 *   - column loop is branch-free (relies on recurr_b being 0 where needed)
 *   - accumulation loop is peeled at m=0 and m=n so the central body is
 *     branch-free and amenable to SIMD reduction
 *   - hot pointers tagged restrict so the auto-vectorizer can reason
 *   - inner loops carry #pragma omp simd hints for GCC/Clang (active when
 *     the project is built with SPODY_ENABLE_OMP_SIMD=ON)
 *
 * Use this in production hot paths (per-step RHS in propagation). The
 * unsuffixed function above remains as the algorithmic reference for
 * regression tests and audit. */
void spody_get_hgaccbodyfixed_hpc(HarmonicGravity *hg, double pos[3], double acc_out[3]);

int spody_load_HarmonicGravityData(HarmonicGravityData *hgd, const char *filename, int degree);
int spody_setup_HarmonicGravity(HarmonicGravity *hg, const HarmonicGravityData *hgd);
int spody_free_HarmonicGravity(HarmonicGravity *hg);
int spody_free_HarmonicGravityData(HarmonicGravityData *hgd);

#ifdef __cplusplus
}
#endif

#endif // SPODY_HARMONICS_H