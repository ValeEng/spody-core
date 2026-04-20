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
#ifndef SPODY_MATH_H
#define SPODY_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

void spody_transpose_matrix(double in[3][3], double out[3][3]); 
void spody_getrotmatrix_XZX(double alpha, double beta, double gamma, double R[3][3]);
void spody_rotate_vector(const double C[3][3], const double v_in[3], double v_out[3]);


#ifdef __cplusplus
}
#endif

#endif // SPODY_MATH_H