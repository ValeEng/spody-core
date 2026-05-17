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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spody_harmonics.h"


static int read_spherical_harmonics_file(FILE *file, HarmonicGravityData *pm, int degree) {
    char line[BUFFER_SIZE_HG];

    if (!fgets(line, sizeof(line), file)) return -1;

    char *token = strtok(line, ",");
    if (!token) return -1;
    //printf("Token N read: %s\n", token);

    pm->R_ref = strtod(token, NULL);

    token = strtok(NULL, ",");
    if (!token) return -1;
    //printf("Token N read: %s\n", token);

    pm->GM = strtod(token, NULL);

    // skip token
    strtok(NULL, ","); 
    
    token = strtok(NULL, ",");
    //printf("Token N read: %s\n", token);
    if (!token) return -1;

    // we read the maximum degree from file
    // but we limit to 'degree' requested by user.
    double file_N = strtod(token, NULL);
    if (degree > (int)file_N)  {
        printf(" Degree %d out of bounds. File maximum degree is %d\n", degree, (int)file_N);
        return -1;
    }
    pm->N = (degree <= (int)file_N) ? degree : (int)file_N;

    #if DEBUG_HARMONICS == 1
    printf("\n\nR_ref = %f | GM = %.10e | N_used = %d\n\n", pm->R_ref, pm->GM, pm->N);
    #endif

    int size = (pm->N + 2) * (pm->N + 3) / 2;
    pm->C = calloc(size, sizeof(double)); // calloc init to 0
    pm->S = calloc(size, sizeof(double));

    int loaded_count = 0;
    while (fgets(line, sizeof(line), file)) {

        token = strtok(line, ","); 
        //printf("Token N read: %s\n", token);

        if (!token) continue;
        
        int n = (int)strtod(token, NULL);
        if (n > pm->N) break; // Stop if beyond desired degree

        token = strtok(NULL, ",");
        int m = (int)strtod(token, NULL);
        //printf("Token N read: %s\n", token);

        token = strtok(NULL, ",");
        double C_val = strtod(token, NULL);
        //printf("Token N read: %s\n", token);
        
        token = strtok(NULL, ",");
        double S_val = strtod(token, NULL);
        //printf("Token N read: %s\n", token);
        //printf("Read coeffs n=%d m=%d | C=%.12e | S=%.12e\n", n, m, C_val, S_val);

        int index = (n * (n + 1) / 2) + m; //tringular indexing
        //printf("index: %d\n", index);

        if (index < size) {
            pm->C[index] = C_val;
            pm->S[index] = S_val;
        }
    }

    return 0;
}

void spody_get_hgaccbodyfixed(HarmonicGravity *hg, double pos[3], double acc_out[3]) {

    //https://hanspeterschaub.info/Papers/Martin2020a.pdf
    
    double r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    double r = sqrt(r2);
    double ir = 1.0 / r; 

    double s = pos[0] * ir; //cos phi * cos lambda
    double t = pos[1] * ir; //cos phi * sin lambda
    double u = pos[2] * ir; //sin phi | phi = latitude
    
    double rho = hg->hgd->R_ref * ir; // R_ref/r
    int N_max = hg->hgd->N; //max degree

    // buffers pointers
    double *real = hg->real;
    double *imag = hg->imag; 
    
    // Setup Rolling Buffers
    // r_n  = current row (n)
    // r_np1 = row +1 (n+1)
    // r_nm1 = row -1 (n-1)
    double *r_np1 = hg->A_row0; 
    double *r_n   = hg->A_row1; 
    double *r_nm1 = hg->A_row2; 
    double *tmp;

    real[0] = 1.0; //real(s + it)^m
    imag[0] = 0.0; //imag(s + it)^m
    if (N_max >= 1) {
        real[1] = s; imag[1] = t;
    }
    for (int m = 2; m <= N_max; m++) { // TBD understand if we can put it in the next loop
        real[m] = s * real[m-1] - t * imag[m-1];
        imag[m] = t * real[m-1] + s * imag[m-1];
    }

    // setup first row of A
    r_n[0] = 1.0; // A[0][0]
    double a1 = 0, a2 = 0, a3 = 0, a4 = 0; // accumulator for acceleration components

    // setup rho_n.
    // n=0 -> rho^2 | n=1 -> rho^3 | n=2 -> rho^4.
    // Sign convention: this prefactor yields acc_out = -grad(V_pert) so
    // callers can simply do dvdt += acc_out (i.e. acc_out IS the acceleration).
    double rho_n = (hg->hgd->GM / (hg->hgd->R_ref * hg->hgd->R_ref)) * (rho * rho * rho);

    // MAIN LOOP (starting from 0 for row 1 generation)
    for (int n = 0; n <= N_max; n++) {
        // calculate row[n+1](r_n1) with row[n] (r_n) 
        // if(n=N_max), calculate a row[N_max+1] for a3 evaluation
          
        int np1 = n + 1; //Warning: don't exceed allocated memery! --> np1 <= N_target + 2 (the current allocation is N+3)
        int idx_diag = (np1 * (np1 + 1) / 2) + np1;
        
        // diagonal
        r_np1[np1] = hg->hgd->recurr_a[idx_diag] * r_n[n]; // A[n+1][n+1] = factorA(n+1) * A[n][n]

        // column (m < n+1)
        for (int m = 0; m <= n; m++) {
            
            int idx = (np1 * (np1 + 1) / 2) + m;
            
            if (np1 - m >= 2){ 
                // A[n+1][m] = factorA(n+1) * u * A[n][m] - factorB(n+1) * A[n-1][m]
                r_np1[m] = hg->hgd->recurr_a[idx] * u * r_n[m] - hg->hgd->recurr_b[idx] * r_nm1[m]; 
            }else{
                r_np1[m] = hg->hgd->recurr_a[idx] * u * r_n[m]; 
            }
        }

        // gravity accumaltion
        // from n >= 2, we need only the disturbing part
        // n is relative to r_n position 
        
        if (n >= 2) {

            double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
            int row_start = n * (n + 1) / 2;

            for (int m = 0; m <= n; m++) {

                int k = row_start + m;

                double C = hg->hgd->C[k];
                double S = hg->hgd->S[k];

                double D = C * real[m] + S * imag[m];

                if (m > 0) {
                    double E = C * real[m-1] + S * imag[m-1];
                    double F = S * real[m-1] - C * imag[m-1];
                    s1 += m * r_n[m] * E;
                    s2 += m * r_n[m] * F;
                }

                if (m != n){ //because r_n[n+1] is not defined and recurr_s3 has m-n=0 at numerator
                    s3 += hg->hgd->recurr_s3[k] * r_n[m+1] * D;
                }

                s4 += hg->hgd->recurr_s4[k] * r_np1[m+1] * D;

                #if DEBUG_HARMONICS == 1
                printf("s1 += %.21e | s2 += %.21e | s3 += %.21e | s4 += %.21e (n=%03d m=%03d)\n", m * r_n[m] * D, m * r_n[m] * (S * real[m] - C * imag[m]), r_np1[m] * D, hg->hgd->recurr_s4[k]* r_np1[m+1] * D, n, m);
                #endif
            }
            //

            rho_n *= rho; 
            a1 += rho_n * s1;
            a2 += rho_n * s2;         
            a3 += rho_n * s3;
            a4 += rho_n * s4;
            #if DEBUG_HARMONICS == 1
            printf("After n=%03d: rho_n=%.21f | a1=%.21e | a2=%.21e | a3=%.21e | a4=%.21e\n", n, rho_n, a1, a2, a3, a4); 
            #endif
        }

        // buffer rolling 
        tmp = r_nm1;
        r_nm1 = r_n;    // n   go to --> n-1
        r_n = r_np1;    // n+1 go to --> n 
        r_np1 = tmp;    // n-1 go to --> n+2
    }
 
    // from Pines coordinates to body fixed
    acc_out[0] = a1 - s * a4;
    acc_out[1] = a2 - t * a4;
    acc_out[2] = a3 - u * a4;

}

/* ============================================================
 * High-performance variant: branch-free + peeled + SIMD hints.
 * Same numerical contract as spody_get_hgaccbodyfixed.
 *
 * To take effect, the project must be built with SPODY_ENABLE_OMP_SIMD=ON
 * (which adds /openmp:experimental on MSVC or -fopenmp-simd on GCC/Clang).
 *
 * MSVC notes (and reason for compiler-gated macros below):
 *   - `/openmp:experimental` ignores `reduction` on `simd` (warning C4849),
 *     so vectorizing the accumulation loop without proper privatization
 *     produces wrong results.
 *   - Even on the column loop (no reduction), MSVC's SIMD output drifts
 *     numerically beyond the test tolerance, presumably from FMA fusion
 *     or fast-math behaviour enabled by `/openmp:experimental`.
 *   - On MSVC we therefore skip the pragmas entirely and let `/O2`
 *     auto-vectorize what it can. On GCC/Clang the explicit pragmas
 *     light up SIMD with the correct reduction semantics.
 * ============================================================ */
#if defined(__clang__) || defined(__GNUC__)
    #define SPODY_SIMD                  _Pragma("omp simd")
    #define SPODY_SIMD_REDUCTION_S1234  _Pragma("omp simd reduction(+:s1,s2,s3,s4)")
#else
    #define SPODY_SIMD                  /* no pragma on MSVC: numerical drift */
    #define SPODY_SIMD_REDUCTION_S1234  /* no pragma on MSVC: reduction unsupported */
#endif

/* Tell the compiler the pointers it tags don't alias.
 * Required for MSVC's auto-vectorizer to clear reason 1502
 * ("loop body too complex"). C99 has `restrict`; MSVC has `__restrict`. */
#if defined(_MSC_VER) && !defined(__clang__)
    #define SPODY_RESTRICT __restrict
#else
    #define SPODY_RESTRICT restrict
#endif
void spody_get_hgaccbodyfixed_hpc(HarmonicGravity *hg, double pos[3], double acc_out[3]) {

    double r2 = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    double r = sqrt(r2);
    double ir = 1.0 / r;

    double s = pos[0] * ir;
    double t = pos[1] * ir;
    double u = pos[2] * ir;

    double rho = hg->hgd->R_ref * ir;
    int N_max = hg->hgd->N;

    double * SPODY_RESTRICT real = hg->real;
    double * SPODY_RESTRICT imag = hg->imag;

    double * SPODY_RESTRICT r_np1 = hg->A_row0;
    double * SPODY_RESTRICT r_n   = hg->A_row1;
    double * SPODY_RESTRICT r_nm1 = hg->A_row2;
    double *tmp;

    real[0] = 1.0;
    imag[0] = 0.0;
    if (N_max >= 1) {
        real[1] = s; imag[1] = t;
    }
    // serial recurrence: cannot SIMD (each iter depends on the previous)
    for (int m = 2; m <= N_max; m++) {
        real[m] = s * real[m-1] - t * imag[m-1];
        imag[m] = t * real[m-1] + s * imag[m-1];
    }

    r_n[0] = 1.0;
    double a1 = 0, a2 = 0, a3 = 0, a4 = 0;

    // Sign convention: see the reference kernel above. acc_out IS the
    // acceleration (-grad V_pert), so callers do dvdt += acc_out.
    double rho_n = (hg->hgd->GM / (hg->hgd->R_ref * hg->hgd->R_ref)) * (rho * rho * rho);

    // pull pointers out so the SIMD loops see plain non-aliasing arrays
    const double * SPODY_RESTRICT recurr_a  = hg->hgd->recurr_a;
    const double * SPODY_RESTRICT recurr_b  = hg->hgd->recurr_b;
    const double * SPODY_RESTRICT recurr_s3 = hg->hgd->recurr_s3;
    const double * SPODY_RESTRICT recurr_s4 = hg->hgd->recurr_s4;
    const double * SPODY_RESTRICT Cv        = hg->hgd->C;
    const double * SPODY_RESTRICT Sv        = hg->hgd->S;

    for (int n = 0; n <= N_max; n++) {

        int np1 = n + 1;
        int row_np1_start = (np1 * (np1 + 1)) / 2;
        int idx_diag = row_np1_start + np1;

        // diagonal
        r_np1[np1] = recurr_a[idx_diag] * r_n[n];

        // ---- column (Phase 2): branch-free + SIMD ----
        // m hoisted outside the for() because MSVC's OpenMP requires the
        // canonical loop form (init-expr cannot be a declaration).
        int m;
        SPODY_SIMD
        for (m = 0; m <= n; m++) {
            int idx = row_np1_start + m;
            r_np1[m] = recurr_a[idx] * u * r_n[m] - recurr_b[idx] * r_nm1[m];
        }

        if (n >= 2) {

            double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
            int row_start = (n * (n + 1)) / 2;

            // ---- m = 0: only s3 and s4 contribute ----
            {
                int k = row_start;
                double C = Cv[k];
                double S = Sv[k];
                double D = C * real[0] + S * imag[0];
                s3 += recurr_s3[k] * r_n[1]   * D;
                s4 += recurr_s4[k] * r_np1[1] * D;
            }

            // ---- m = 1..n-1: branch-free + SIMD reduction ----
            SPODY_SIMD_REDUCTION_S1234
            for (m = 1; m < n; m++) {
                int k = row_start + m;
                double C = Cv[k];
                double S = Sv[k];

                double D = C * real[m]   + S * imag[m];
                double E = C * real[m-1] + S * imag[m-1];
                double F = S * real[m-1] - C * imag[m-1];

                s1 += m * r_n[m] * E;
                s2 += m * r_n[m] * F;
                s3 += recurr_s3[k] * r_n[m+1]   * D;
                s4 += recurr_s4[k] * r_np1[m+1] * D;
            }

            // ---- m = n: s1, s2, s4 (recurr_s3 is 0 there) ----
            {
                int k = row_start + n;
                double C = Cv[k];
                double S = Sv[k];
                double D = C * real[n]   + S * imag[n];
                double E = C * real[n-1] + S * imag[n-1];
                double F = S * real[n-1] - C * imag[n-1];
                s1 += n * r_n[n] * E;
                s2 += n * r_n[n] * F;
                s4 += recurr_s4[k] * r_np1[n+1] * D;
            }

            rho_n *= rho;
            a1 += rho_n * s1;
            a2 += rho_n * s2;
            a3 += rho_n * s3;
            a4 += rho_n * s4;
        }

        tmp = r_nm1;
        r_nm1 = r_n;
        r_n   = r_np1;
        r_np1 = tmp;
    }

    acc_out[0] = a1 - s * a4;
    acc_out[1] = a2 - t * a4;
    acc_out[2] = a3 - u * a4;
}

int spody_load_HarmonicGravityData(HarmonicGravityData *hgd, const char *filename, int degree){
    
    FILE *file = fopen(filename,"r");
    if (!file) {
        perror("Errore apertura file");
        return -1;
    }

    // 2. Leggere il file
    if (read_spherical_harmonics_file(file, hgd, degree) != 0) {
        fclose(file);
        return -1;
    }

    fclose(file);

    int file_degree = hgd->N; //to be sure we use the correct degree in all functions
    int size = (file_degree + 2) * (file_degree + 3) / 2;
    hgd->recurr_a = malloc(size * sizeof(double));
    hgd->recurr_b = malloc(size * sizeof(double));
    hgd->recurr_s3 = malloc(size * sizeof(double));
    hgd->recurr_s4 = malloc(size * sizeof(double));
    memset(hgd->recurr_a, 0, size * sizeof(double));
    memset(hgd->recurr_b, 0, size * sizeof(double));
    memset(hgd->recurr_s3, 0, size * sizeof(double));
    memset(hgd->recurr_s4, 0, size * sizeof(double));

    for (int n = 0; n <= file_degree + 1; n++) {

        for (int m = 0; m <= n; m++) {
            int i = (n * (n + 1) / 2) + m;

            if (n == m) {
                //printf("Calculating recurr_a for n=%d m=%d at index %d\n", n, m, i);
                // diagonal elements (sectorial)
                if (n == 0) {
                    hgd->recurr_a[i] = 1.0; //start
                    //printf("Set recurr_a[%d] = %.12f\n", i, hg->recurr_a[i]);
                }else if (n == 1) {
                    hgd->recurr_a[i] = sqrt(3.0); // from Pines 19.. DeMars 2008 and GMAT
                }else {
                    // diagonal elements from A(n-1, n-1) to A(n, n)
                    hgd->recurr_a[i] = sqrt((double)(2.0 * n + 1.0) / (2.0 * n)); //the same : sqrt((double)(1.0 + (1.0 / (2.0 * n))));
                }
            }else {
                //printf("Calculating recurr_a and recurr_b for n=%d m=%d at index %d\n", n, m, i);
                // column (Tesseral and zonal)
                hgd->recurr_a[i] = sqrt((double)((2.0 * n + 1.0) * (2.0 * n - 1.0)) / ((n - m) * (n + m)));
                // multiplay : u * A[n-1][m]

                if (n - m >= 2) {
                    hgd->recurr_b[i] = sqrt((double)((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0)) /
                                            ((2.0 * n - 3.0) * (n + m) * (n - m)));
                    // multiplay : A[n-2][m]
                } else {
                    // we dont need this brench thanks to the memset
                    //hg->recurr_b[i] = 0.0;
                }
            }

            // gravity-accumulation factors (precomputed, only n >= 2 are used)
            if (n >= 2) {
                double delta = (m == 0) ? 1.0 : 0.0;
                if (m != n) {
                    hgd->recurr_s3[i] = sqrt((n - m) * (n + m + 1.0) * (2.0 - delta) / 2.0);
                }
                hgd->recurr_s4[i] = sqrt((n + m + 2.0) * (n + m + 1.0) * (2.0 * n + 1.0) * (2.0 - delta) /
                                        (2.0 * (2.0 * n + 3.0)));
            }
        }
    }


    return 0;

}

int spody_setup_HarmonicGravity(HarmonicGravity *hg, HarmonicGravityData *hgd){
    //TBD clean sizes to improve memory usage e performance
    // we have allocated more memory than necessary for simplicity (N+2)
    
    hg->hgd = hgd;
    int degree = hg->hgd->N;
        
    int r_size = (degree + 3) * sizeof(double);
    hg->A_row0 = malloc(r_size);
    hg->A_row1 = malloc(r_size);
    hg->A_row2 = malloc(r_size);
    memset(hg->A_row0, 0, r_size);
    memset(hg->A_row1, 0, r_size);
    memset(hg->A_row2, 0, r_size);
    //TBD better calloc ? 

    int dim = degree + 2; 
    hg->real = malloc(dim * sizeof(double));
    hg->imag = malloc(dim * sizeof(double));
    
    return 0;
}

int spody_free_HarmonicGravity(HarmonicGravity *hg) {
    if (hg->A_row0) free(hg->A_row0);
    if (hg->A_row1) free(hg->A_row1);
    if (hg->A_row2) free(hg->A_row2);
    if (hg->imag) free(hg->imag);
    if (hg->real) free(hg->real);

    return 0;
}

int spody_free_HarmonicGravityData(HarmonicGravityData *hgd) {
    if (hgd->C) free(hgd->C);
    if (hgd->S) free(hgd->S);
    if (hgd->recurr_a) free(hgd->recurr_a);
    if (hgd->recurr_b) free(hgd->recurr_b);
    if (hgd->recurr_s3) free(hgd->recurr_s3);
    if (hgd->recurr_s4) free(hgd->recurr_s4);
    return 0;
}