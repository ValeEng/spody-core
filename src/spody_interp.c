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
#include "spody_interp.h"

/*
 * Cubic Hermite basis on s in [0, 1]:
 *
 *   H00(s) =  2 s^3 - 3 s^2 + 1     value at endpoint a
 *   H10(s) =    s^3 - 2 s^2 + s     slope at endpoint a (scaled by h)
 *   H01(s) = -2 s^3 + 3 s^2         value at endpoint b
 *   H11(s) =    s^3 -   s^2         slope at endpoint b (scaled by h)
 *
 * y(t) = H00 y_a + H10 h y'_a + H01 y_b + H11 h y'_b
 *
 * Derivative is exact for polynomials of degree <= 3 given matching
 * endpoint values and slopes (this is what makes it the natural
 * dense-output companion to 5(4) and 8(7) RK pairs).
 */
static void hermite_basis(double s, double *H00, double *H10,
                          double *H01, double *H11) {
    double s2 = s * s;
    double s3 = s2 * s;
    *H00 =  2.0 * s3 - 3.0 * s2 + 1.0;
    *H10 =        s3 - 2.0 * s2 + s;
    *H01 = -2.0 * s3 + 3.0 * s2;
    *H11 =        s3 -       s2;
}

/* Derivative of the basis with respect to s.
 * dH00/ds, dH10/ds, dH01/ds, dH11/ds. */
static void hermite_basis_d(double s, double *dH00, double *dH10,
                            double *dH01, double *dH11) {
    double s2 = s * s;
    *dH00 =  6.0 * s2 - 6.0 * s;
    *dH10 =  3.0 * s2 - 4.0 * s + 1.0;
    *dH01 = -6.0 * s2 + 6.0 * s;
    *dH11 =  3.0 * s2 - 2.0 * s;
}

/* Clamp s into [0, 1] so that queries slightly outside the interval
 * (typical floating-point edge case at grid boundaries) return the
 * endpoint instead of extrapolating. */
static double clamp_unit(double s) {
    if (s < 0.0) return 0.0;
    if (s > 1.0) return 1.0;
    return s;
}

double spody_hermite_cubic_1d(double t,
                              double t_a, double t_b,
                              double y_a, double dy_a,
                              double y_b, double dy_b) {
    double h = t_b - t_a;
    double s = clamp_unit((t - t_a) / h);
    double H00, H10, H01, H11;
    hermite_basis(s, &H00, &H10, &H01, &H11);
    return H00 * y_a + H10 * h * dy_a
         + H01 * y_b + H11 * h * dy_b;
}

void spody_hermite_cubic(double t,
                         double t_a, double t_b,
                         const double *y_a, const double *dy_a,
                         const double *y_b, const double *dy_b,
                         int dim,
                         double *y_out, double *dy_out) {
    double h = t_b - t_a;
    double s = clamp_unit((t - t_a) / h);

    double H00, H10, H01, H11;
    hermite_basis(s, &H00, &H10, &H01, &H11);
    for (int k = 0; k < dim; ++k) {
        y_out[k] = H00 * y_a[k] + H10 * h * dy_a[k]
                 + H01 * y_b[k] + H11 * h * dy_b[k];
    }

    if (dy_out) {
        double dH00, dH10, dH01, dH11;
        hermite_basis_d(s, &dH00, &dH10, &dH01, &dH11);
        /* dy/dt = (1/h) * dy/ds. The basis derivatives H00..H11 carry
         * a factor of h on the slope-at-endpoint terms; dividing by h
         * the value terms gives the consistent dy/dt expression. */
        for (int k = 0; k < dim; ++k) {
            dy_out[k] = (dH00 * y_a[k] + dH01 * y_b[k]) / h
                      +  dH10 * dy_a[k] + dH11 * dy_b[k];
        }
    }
}

void spody_hermite_dense_rv6(double t,
                             double t_a, const double y_a[6],
                             double t_b, const double y_b[6],
                             double y_out[6]) {
    /* Position is interpolated via cubic Hermite using the velocity
     * portion of y as the slope; velocity is the analytic derivative
     * of that cubic, so r and v stay C^1 consistent. */
    spody_hermite_cubic(t, t_a, t_b,
                        y_a,     y_a + 3,
                        y_b,     y_b + 3,
                        3,
                        y_out,   y_out + 3);
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
