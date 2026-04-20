#ifndef SPODY_MATH_H
#define SPODY_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

void spody_transpose_matrix(double in[3][3], double out[3][3]); //TBD dont work 
void spody_getrotmatrix_XZX(double alpha, double beta, double gamma, double R[3][3]);
void spody_rotate_vector(const double C[3][3], const double v_in[3], double v_out[3]);


#ifdef __cplusplus
}
#endif

#endif // SPODY_MATH_H