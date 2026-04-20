#include "spody_math.h"
#include <math.h>

void spody_transpose_matrix(double in[3][3], double out[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[j][i] = in[i][j];
        }
    }
}

void spody_getrotmatrix_XZX(double alpha, double beta, double gamma, double R[3][3]) {
    // Rotation matrix for XZX Euler angles
    double ca = cos(alpha);
    double sa = sin(alpha);
    double cb = cos(beta);
    double sb = sin(beta);
    double cg = cos(gamma);
    double sg = sin(gamma);

    R[0][0] = ca * cg - sa * cb * sg;
    R[0][1] = -ca * sg - sa * cb * cg;
    R[0][2] = sa * sb;

    R[1][0] = sa * cg + ca * cb * sg;
    R[1][1] = -sa * sg + ca * cb * cg;
    R[1][2] = -ca * sb;

    R[2][0] = sb * sg;
    R[2][1] = sb * cg;
    R[2][2] = cb;
}

void spody_rotate_vector(const double C[3][3], const double v_in[3], double v_out[3]){
  
    v_out[0] = C[0][0]*v_in[0] + C[0][1]*v_in[1] + C[0][2]*v_in[2];

    v_out[1] = C[1][0]*v_in[0] + C[1][1]*v_in[1] + C[1][2]*v_in[2];

    v_out[2] = C[2][0]*v_in[0] + C[2][1]*v_in[1] + C[2][2]*v_in[2];
    
}