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
#include "spody_eclipse.h"

/* ==========================================================================
 * Eclipse geometry -- Montenbruck & Gill, GMAT-like.
 *
 * Everything below lives on the sky plane seen from the satellite: the
 * Sun is a disc of angular radius a, every occulting body a disc of
 * angular radius b_i, and shadowing is circle overlap. The functions
 * work in FRACTIONS of the solar disc rather than absolute areas --
 * every term of the inclusion-exclusion sum carries the same PI*a^2
 * denominator, so the factor cancels and never enters the arithmetic.
 * ========================================================================== */

/* Pull an argument back into the acos() domain. Dot products and
 * ratios land a few ULP outside [-1, 1] on nearly degenerate
 * geometry. */
static double clamp1(double x) {
    if (x >  1.0) return  1.0;
    if (x < -1.0) return -1.0;
    return x;
}

/* Overlap area of two discs of radii r1, r2 whose centres are d
 * apart. Radii and distance are angles [rad], the result is an
 * angular area [rad^2]. */
static double disc_overlap_area(double r1, double r2, double d) {
    if (d >= r1 + r2) return 0.0;                 /* disjoint          */
    if (d <= r2 - r1) return PI * r1 * r1;        /* disc 1 inside 2   */
    if (d <= r1 - r2) return PI * r2 * r2;        /* disc 2 inside 1   */

    double aa = r1 * r1;
    double bb = r2 * r2;
    double x  = (d * d + aa - bb) / (2 * d);
    double y  = sqrt(fmax(0.0, aa - x * x));
    return (aa * acos(clamp1(x / r1)))
         + (bb * acos(clamp1((d - x) / r2)))
         - (d * y);
}

/* Fraction of the solar disc (angular radius a) hidden by ONE
 * occulting disc (angular radius b) whose centre sits at angular
 * distance c from the Sun's centre. The three degenerate
 * configurations return exact values rather than going through the
 * area formula, so no spurious rounding is introduced. */
static double disc_overlap_fraction(double a, double b, double c) {
    if (a + b <= c) return 0.0;                   /* full sunlight     */
    if (c <= b - a) return 1.0;                   /* total eclipse     */
    if (c <= a - b) return (b * b) / (a * a);     /* annular           */
    return disc_overlap_area(a, b, c) / (PI * (a * a));
}

/* Common area of three discs given their radii and centres on the
 * plane. Used for the pairwise correction of the inclusion-exclusion
 * sum, where disc 0 is always the Sun.
 *
 * The general configuration -- three circles crossing pairwise -- has
 * a boundary made of one arc per circle meeting at three vertices.
 * The area then follows from Green's theorem, 1/2 * contour integral
 * of (x dy - y dx), evaluated arc by arc. Unlike the classic
 * "triangle plus three circular segments" recipe this stays correct
 * when an arc is longer than half its circle, which happens as soon
 * as two of the discs nearly coincide. */
static double triple_overlap_area(const double r[3], const double O[3][2]) {
    /* Pairwise distances. One empty pair empties the whole thing. */
    double d[3][3] = { { 0.0 } };
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double dx = O[j][0] - O[i][0];
            double dy = O[j][1] - O[i][1];
            d[i][j] = d[j][i] = sqrt(dx * dx + dy * dy);
            if (d[i][j] >= r[i] + r[j]) return 0.0;
        }
    }

    /* A disc swallowed by another one reduces the triple overlap to
     * the two-disc lens of the survivors. */
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            if (d[i][j] + r[i] <= r[j]) {
                /* Lowest index first so that when the survivor pair
                 * contains the Sun (index 0) the call is the very same
                 * expression disc_overlap_fraction evaluates, and two
                 * coincident bodies cancel to the last bit. */
                int k  = 3 - i - j;
                int lo = (i < k) ? i : k;
                int hi = (i < k) ? k : i;
                return disc_overlap_area(r[lo], r[hi], d[lo][hi]);
            }
        }
    }

    /* Vertices: the pairwise intersection points that fall inside the
     * third disc. A convex region bounded by three circles has at
     * most three of them. */
    double vx[4], vy[4];
    int    va[4], vb[4];               /* the two circles each vertex is on */
    int    nv = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            int    k  = 3 - i - j;
            double dd = d[i][j];
            double ex = (O[j][0] - O[i][0]) / dd;
            double ey = (O[j][1] - O[i][1]) / dd;
            double t  = (dd * dd + r[i] * r[i] - r[j] * r[j]) / (2.0 * dd);
            double h  = sqrt(fmax(0.0, r[i] * r[i] - t * t));
            for (int s = -1; s <= 1; s += 2) {
                double px = O[i][0] + t * ex - s * h * ey;
                double py = O[i][1] + t * ey + s * h * ex;
                double qx = px - O[k][0];
                double qy = py - O[k][1];
                if (qx * qx + qy * qy > r[k] * r[k]) continue;
                if (nv >= 4) return 0.0;        /* degenerate: subtract nothing */
                vx[nv] = px; vy[nv] = py;
                va[nv] = i;  vb[nv] = j;
                nv++;
            }
        }
    }

    if (nv == 0) return 0.0;
    if (nv == 2) {
        /* One lens sits entirely inside the third disc: the triple
         * overlap IS that lens. */
        return disc_overlap_area(r[va[0]], r[vb[0]], d[va[0]][vb[0]]);
    }
    if (nv != 3) return 0.0;           /* numerically degenerate */

    /* Counter-clockwise order around the vertex centroid, which lies
     * inside the (convex) region. */
    double cx  = (vx[0] + vx[1] + vx[2]) / 3.0;
    double cy  = (vy[0] + vy[1] + vy[2]) / 3.0;
    double ang[3];
    int    idx[3] = { 0, 1, 2 };
    for (int i = 0; i < 3; ++i) ang[i] = atan2(vy[i] - cy, vx[i] - cx);
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (ang[idx[j]] < ang[idx[i]]) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }
        }
    }

    /* Walk the boundary. Consecutive vertices share exactly one
     * circle -- that is the arc between them -- and the region being
     * convex, every arc runs counter-clockwise about its own centre,
     * so the CCW angular gap is the one to take. */
    double area = 0.0;
    for (int s = 0; s < 3; ++s) {
        int p = idx[s];
        int q = idx[(s + 1) % 3];
        int k = (va[p] == va[q] || va[p] == vb[q]) ? va[p] : vb[p];
        double t1 = atan2(vy[p] - O[k][1], vx[p] - O[k][0]);
        double t2 = atan2(vy[q] - O[k][1], vx[q] - O[k][0]);
        double dt = t2 - t1;
        while (dt < 0.0) dt += 2.0 * PI;
        area += 0.5 * (r[k] * r[k] * dt
                       + O[k][0] * r[k] * (sin(t2) - sin(t1))
                       - O[k][1] * r[k] * (cos(t2) - cos(t1)));
    }
    return fmax(0.0, area);
}

/* Fraction of the solar disc hidden by TWO occulting discs at the
 * same time. Inputs are the three angular radii and the three angular
 * separations already computed by the caller; the sky-plane layout is
 * reconstructed from them with the Sun at the origin and body 1 on
 * the x axis. */
static double triple_overlap_fraction(double a, double b1, double b2,
                                      double c1, double c2, double c12) {
    double r[3] = { a, b1, b2 };
    double O[3][2] = { { 0.0, 0.0 }, { c1, 0.0 }, { 0.0, 0.0 } };

    if (c1 > 0.0) {
        double x3 = (c1 * c1 + c2 * c2 - c12 * c12) / (2.0 * c1);
        O[2][0] = x3;
        O[2][1] = sqrt(fmax(0.0, c2 * c2 - x3 * x3));
    } else {
        /* Body 1 dead centre on the Sun: the triangle collapses onto
         * a line and the construction above would divide by zero. */
        O[2][0] = c2;
    }
    return triple_overlap_area(r, O) / (PI * (a * a));
}

double spody_get_satlitfraction(const double sat2sun[3], double sun_radius,
                                const double sat2occ[][3],
                                const double occ_radius[], int n_occ) {
    if (!sat2sun || !sat2occ || !occ_radius || n_occ <= 0) return 1.0;
    if (sun_radius <= 0.0) return 1.0;
    if (n_occ > SPODY_ECL_MAX_OCCULTERS) n_occ = SPODY_ECL_MAX_OCCULTERS;

    double d_sun = sqrt(sat2sun[0] * sat2sun[0]
                      + sat2sun[1] * sat2sun[1]
                      + sat2sun[2] * sat2sun[2]);
    if (!(d_sun > 0.0)) return 1.0;

    double usun[3];
    usun[0] = sat2sun[0] / d_sun;
    usun[1] = sat2sun[1] / d_sun;
    usun[2] = sat2sun[2] / d_sun;
    double a = asin(clamp1(sun_radius / d_sun));

    /* Active list: the bodies actually covering part of the Sun right
     * now. Everything downstream -- the sum, the pairs -- runs on this
     * list, not on the full one. */
    double g[SPODY_ECL_MAX_OCCULTERS];      /* hidden fraction   */
    double b[SPODY_ECL_MAX_OCCULTERS];      /* angular radius    */
    double c[SPODY_ECL_MAX_OCCULTERS];      /* separation to Sun */
    double u[SPODY_ECL_MAX_OCCULTERS][3];   /* unit sat -> body  */
    int    k = 0;

    for (int i = 0; i < n_occ; ++i) {
        if (occ_radius[i] <= 0.0) continue;
        const double *p = sat2occ[i];

        /* Screening: is the satellite on the sunward side of this
         * body? Same decision M&G open with, minus the normalisation
         * a sign test does not need. Exact -- the penumbra reaches
         * the plane through the body centre only below an altitude of
         * R*(1/cos(a) - 1), some 70 m for the Earth. */
        double d2  = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        double dps = p[0] * sat2sun[0] + p[1] * sat2sun[1] + p[2] * sat2sun[2];
        if (d2 > dps) continue;

        double d = sqrt(d2);
        if (occ_radius[i] >= d) return 0.0;     /* satellite inside the body */

        double b_i = asin(clamp1(occ_radius[i] / d));
        double u_i[3];
        u_i[0] = p[0] / d;
        u_i[1] = p[1] / d;
        u_i[2] = p[2] / d;
        double c_i = acos(clamp1((p[0] * usun[0] + p[1] * usun[1]
                                                 + p[2] * usun[2]) / d));
        double g_i = disc_overlap_fraction(a, b_i, c_i);

        #if DEBUG_ECLIPSE == 1
        printf("[eclipse] body %d | a: %.12f | b: %.12f | c: %.12f | hidden: %.12f\n",
               i, a, b_i, c_i, g_i);
        #endif

        if (g_i <= 0.0) continue;               /* discs do not overlap */
        if (g_i >= 1.0) return 0.0;             /* total eclipse        */

        g[k] = g_i; b[k] = b_i; c[k] = c_i;
        u[k][0] = u_i[0]; u[k][1] = u_i[1]; u[k][2] = u_i[2];
        k++;
    }

    if (k == 0) return 1.0;
    if (k == 1) return 1.0 - g[0];

    /* Two or more bodies over the Sun at once: their shadows on the
     * solar disc may overlap, so the plain sum would count the shared
     * part twice. Subtract it pair by pair. */
    double hidden = 0.0;
    double gmax   = 0.0;
    for (int i = 0; i < k; ++i) {
        hidden += g[i];
        if (g[i] > gmax) gmax = g[i];
    }
    for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
            /* Separation from the CHORD, not from acos(dot): the two
             * bodies are nearly aligned exactly when the correction
             * matters, and acos loses half its digits there. The
             * reconstruction of the sky-plane triangle downstream is
             * sensitive enough that the difference decides whether two
             * superposed discs are seen as one. */
            double dx   = u[i][0] - u[j][0];
            double dy   = u[i][1] - u[j][1];
            double dz   = u[i][2] - u[j][2];
            double c_ij = 2.0 * asin(clamp1(sqrt(dx*dx + dy*dy + dz*dz) / 2.0));
            if (c_ij >= b[i] + b[j]) continue;  /* the two discs miss */
            hidden -= triple_overlap_fraction(a, b[i], b[j],
                                              c[i], c[j], c_ij);
        }
    }

    /* Keep the truncated series inside its rigorous bracket. */
    double lit = 1.0 - hidden;
    if (lit > 1.0 - gmax) lit = 1.0 - gmax;
    if (lit < 0.0)        lit = 0.0;
    return lit;
}
