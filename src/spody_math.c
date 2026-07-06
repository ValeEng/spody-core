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
#include "spody_math.h"
#include <math.h>
#include "spody_const.h"   /* PI */

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

void spody_bf_to_geodetic(const double r_bf_km[3], double a_km,
                            double inv_f, double *lat_rad,
                            double *lon_rad, double *alt_km) {
    double x = r_bf_km[0], y = r_bf_km[1], z = r_bf_km[2];
    double f   = 1.0 / inv_f;
    double e2  = f * (2.0 - f);          /* first eccentricity^2  */
    double b   = a_km * (1.0 - f);       /* semi-minor axis       */
    double ep2 = e2 / (1.0 - e2);        /* second eccentricity^2 */
    double p   = sqrt(x * x + y * y);

    if (lon_rad) *lon_rad = atan2(y, x);

    /* Polar guard: within ~1 m of the spin axis the Bowring quotient
     * degenerates; the exact answer there is trivial. */
    if (p < 1e-3) {
        if (lat_rad) *lat_rad = (z >= 0.0) ? PI / 2.0 : -PI / 2.0;
        if (alt_km)  *alt_km  = fabs(z) - b;
        return;
    }

    /* Bowring (1976) with a fixed re-anchoring of the parametric
     * latitude: the classic single evaluation is ~1e-9 rad near the
     * surface but degrades to ~5e-9 rad / 20 cm at GEO heights; two
     * extra passes converge to machine precision for any Earth-orbit
     * altitude at the cost of a handful of trig calls. Fixed count,
     * branch-free. */
    {
        double theta = atan2(z * a_km, p * b);
        double lat = 0.0;
        double sl, cl, N;
        int it;
        for (it = 0; it < 3; ++it) {
            double st = sin(theta), ct = cos(theta);
            lat = atan2(z + ep2 * b * st * st * st,
                        p - e2 * a_km * ct * ct * ct);
            theta = atan2((1.0 - f) * sin(lat), cos(lat));
        }
        sl = sin(lat);
        cl = cos(lat);
        N = a_km / sqrt(1.0 - e2 * sl * sl);
        if (lat_rad) *lat_rad = lat;
        if (alt_km) {
            /* Away from the poles p/cos is well conditioned; near
             * them (|lat| > ~80 deg) the z-form is the stable one. */
            if (fabs(cl) > 0.17)
                *alt_km = p / cl - N;
            else
                *alt_km = z / sl - N * (1.0 - e2);
        }
    }
}
size_t spody_bracket_index(const double *xs, size_t n, double x) {
    if (x <= xs[0]) return 0;
    if (x >= xs[n - 1]) return n - 2;
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (xs[mid] <= x) lo = mid;
        else hi = mid;
    }
    return lo;
}

double spody_interp_linear(const double *xs, const double *ys,
                           size_t n, double x) {
    if (n == 1 || x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    size_t i = spody_bracket_index(xs, n, x);
    double t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    return ys[i] + t * (ys[i + 1] - ys[i]);
}

double spody_dot3(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void spody_cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}
