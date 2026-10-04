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
#include <stdlib.h>
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
double spody_dot3(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void spody_cross3(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

void spody_iau_pole(double ra0_deg, double ra1_deg, double dec0_deg,
                    double dec1_deg, double et, double pole[3]) {
    double T   = et / (SECONDSxDAY * DAYS_PER_JULIAN_CY);
    double ra  = (ra0_deg  + ra1_deg  * T) * DEG2RAD;
    double dec = (dec0_deg + dec1_deg * T) * DEG2RAD;
    pole[0] = cos(dec) * cos(ra);
    pole[1] = cos(dec) * sin(ra);
    pole[2] = sin(dec);
}

double spody_body_shape_distance(const SpodyBodyShape *shape,
                                 const double r_rel[3]) {
    double d2 = r_rel[0] * r_rel[0] + r_rel[1] * r_rel[1] + r_rel[2] * r_rel[2];
    if (!(shape->r_pol > 0.0) || shape->r_pol == shape->r_eq) return sqrt(d2);

    /* The spheroid is symmetric about its axis: the point's height
     * along the pole and its distance from the axis are all
     * spody_bf_to_geodetic needs. */
    double z = spody_dot3(r_rel, shape->pole);
    double r_axis[3] = { sqrt(fmax(0.0, d2 - z * z)), 0.0, z };
    double alt_km = 0.0;
    spody_bf_to_geodetic(r_axis, shape->r_eq,
                         shape->r_eq / (shape->r_eq - shape->r_pol),
                         NULL, NULL, &alt_km);
    return shape->r_eq + alt_km;
}

/* Degeneracy thresholds of the RIC axes (copied verbatim by the
 * Python twin spopy.rotations.ric_to_icrf). */
static const double ric_min_r_km = 1.0e-9;
static const double ric_min_h    = 1.0e-12;

int spody_getrotmatrix_ric2icrf(const double r[3], const double v[3],
                                double R[3][3]) {
    double rn = sqrt(spody_dot3(r, r));
    if (rn < ric_min_r_km) return -1;
    double h[3];
    spody_cross3(r, v, h);
    double hn = sqrt(spody_dot3(h, h));
    if (hn < ric_min_h) return -1;
    double rh[3] = { r[0] / rn, r[1] / rn, r[2] / rn };
    double ch[3] = { h[0] / hn, h[1] / hn, h[2] / hn };
    double ih[3];
    spody_cross3(ch, rh, ih);
    for (int k = 0; k < 3; ++k) {
        R[k][0] = rh[k];
        R[k][1] = ih[k];
        R[k][2] = ch[k];
    }
    return 0;
}

int spody_getrotmatrix_icrf2ric(const double r[3], const double v[3],
                                double R[3][3]) {
    double Rt[3][3];
    if (spody_getrotmatrix_ric2icrf(r, v, Rt) != 0) return -1;
    spody_transpose_matrix(Rt, R);
    return 0;
}

int spody_symmat_cholesky(int n, const double *a, double *l) {
    for (int i = 0; i < n * n; ++i) l[i] = 0.0;
    for (int j = 0; j < n; ++j) {
        double d = a[j * n + j];
        for (int k = 0; k < j; ++k) d -= l[j * n + k] * l[j * n + k];
        if (!(d > 0.0)) return j + 1;          /* also catches NaN */
        double ljj = sqrt(d);
        l[j * n + j] = ljj;
        for (int i = j + 1; i < n; ++i) {
            double s = a[i * n + j];
            for (int k = 0; k < j; ++k) s -= l[i * n + k] * l[j * n + k];
            l[i * n + j] = s / ljj;
        }
    }
    return 0;
}

/* Cyclic Jacobi: each rotation in the (p, q) plane zeroes a[p][q]
 * (Golub & Van Loan, Matrix Computations, 4th ed., algorithms 8.5.1 and
 * 8.5.3). The angle uses the smaller root t of t^2 + 2 theta t - 1 = 0,
 * theta = (a_qq - a_pp) / (2 a_pq), which keeps |rotation| <= pi/4. */
int spody_symmat_eigen_jacobi(int n, const double *a, double *w, double *v) {
    double *m = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    double *e = v ? v : (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!m || !e) {
        free(m);
        if (!v) free(e);
        return -2;
    }
    for (int i = 0; i < n; ++i)               /* symmetric copy of the lower triangle */
        for (int j = 0; j <= i; ++j)
            m[i * n + j] = m[j * n + i] = a[i * n + j];
    for (int i = 0; i < n * n; ++i) e[i] = 0.0;
    for (int i = 0; i < n; ++i) e[i * n + i] = 1.0;

    int rc = -1;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) off += fabs(m[p * n + q]);
        if (off == 0.0) { rc = 0; break; }
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double apq = m[p * n + q];
                double g = 100.0 * fabs(apq);
                double app = m[p * n + p], aqq = m[q * n + q];
                if (sweep > 3 && fabs(app) + g == fabs(app)
                              && fabs(aqq) + g == fabs(aqq)) {
                    m[p * n + q] = m[q * n + p] = 0.0;
                    continue;
                }
                if (apq == 0.0) continue;
                double theta = (aqq - app) / (2.0 * apq);
                double t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                if (theta < 0.0) t = -t;
                double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
                for (int k = 0; k < n; ++k) {  /* columns p, q of m and e */
                    double mkp = m[k * n + p], mkq = m[k * n + q];
                    m[k * n + p] = c * mkp - s * mkq;
                    m[k * n + q] = s * mkp + c * mkq;
                    double ekp = e[k * n + p], ekq = e[k * n + q];
                    e[k * n + p] = c * ekp - s * ekq;
                    e[k * n + q] = s * ekp + c * ekq;
                }
                for (int k = 0; k < n; ++k) {  /* rows p, q of m */
                    double mpk = m[p * n + k], mqk = m[q * n + k];
                    m[p * n + k] = c * mpk - s * mqk;
                    m[q * n + k] = s * mpk + c * mqk;
                }
                m[p * n + q] = m[q * n + p] = 0.0;
            }
        }
    }

    for (int i = 0; i < n; ++i) w[i] = m[i * n + i];
    for (int i = 1; i < n; ++i) {             /* insertion sort, ascending */
        for (int j = i; j > 0 && w[j - 1] > w[j]; --j) {
            double tw = w[j]; w[j] = w[j - 1]; w[j - 1] = tw;
            for (int k = 0; k < n; ++k) {
                double te = e[k * n + j];
                e[k * n + j] = e[k * n + j - 1];
                e[k * n + j - 1] = te;
            }
        }
    }
    free(m);
    if (!v) free(e);
    return rc;
}
