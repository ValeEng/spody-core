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
#include "spody_kepler.h"
#include "spody_const.h"

#include <math.h>

/* Wrap any real angle into [-pi, pi]. Used to keep the Newton seed
 * close to the true root so the iteration converges in a handful of
 * steps regardless of how the caller phased M. */
static double _wrap_pi(double a) {
    const double TWO_PI = 2.0 * PI;
    a = fmod(a, TWO_PI);
    if (a >  PI) a -= TWO_PI;
    if (a < -PI) a += TWO_PI;
    return a;
}

double spody_kepler_solve_E(double mean_anom_rad, double ecc) {
    double M = _wrap_pi(mean_anom_rad);
    /* Standard initial guess: E0 = M + e*sin(M) (Danby seed). For
     * e < ~0.6 this gets us within machine epsilon in 3-5 Newton
     * steps. The 30-iteration cap is a hard safety; we never hit it
     * in the elliptical regime this module documents. */
    double E = M + ecc * sin(M);
    for (int k = 0; k < 30; ++k) {
        double f  = E - ecc * sin(E) - M;
        double fp = 1.0 - ecc * cos(E);
        double dE = f / fp;
        E -= dE;
        if (fabs(dE) < 1.0e-14) break;
    }
    return E;
}

double spody_kepler_mean_to_true_anom(double mean_anom_rad, double ecc) {
    double E = spody_kepler_solve_E(mean_anom_rad, ecc);
    /* atan2 form avoids quadrant ambiguity. The standard half-angle
     * relation tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2) is equivalent
     * to nu = atan2(sqrt(1-e^2) * sin(E), cos(E) - e). */
    double sqrt_1me2 = sqrt(1.0 - ecc * ecc);
    return atan2(sqrt_1me2 * sin(E), cos(E) - ecc);
}

void spody_keplerian_to_cartesian(double sma_km, double ecc,
                                  double inc_rad, double raan_rad,
                                  double argp_rad, double true_anom_rad,
                                  double mu_km3_s2,
                                  double r_out[3], double v_out[3]) {
    /* Perifocal coordinates: orbit in the (xp, yp) plane, periapsis
     * along +xp. Closed-form from sma, e, nu and mu (Vallado, ch.2). */
    double p     = sma_km * (1.0 - ecc * ecc);     /* semi-latus rectum */
    double cnu   = cos(true_anom_rad);
    double snu   = sin(true_anom_rad);
    double r_mag = p / (1.0 + ecc * cnu);
    double r_pf[3] = { r_mag * cnu, r_mag * snu, 0.0 };

    double mu_p  = sqrt(mu_km3_s2 / p);
    double v_pf[3] = { -mu_p * snu, mu_p * (ecc + cnu), 0.0 };

    /* Rotation from perifocal to the reference inertial frame via the
     * 3-1-3 Euler sequence (raan about z, inc about x', argp about
     * z''). The combined rotation R = R3(raan) * R1(inc) * R3(argp)
     * is expanded below to avoid 3 sequential matvecs.
     *
     * Identical to Vallado eqs. (2-118 / 2-145) and the variant of the
     * formulas used in every aerospace toolbox; sign convention is the
     * mathematical one (counter-clockwise positive). */
    double co = cos(raan_rad), so = sin(raan_rad);
    double ci = cos(inc_rad),  si = sin(inc_rad);
    double cw = cos(argp_rad), sw = sin(argp_rad);

    double R[3][3];
    R[0][0] =  co * cw - so * sw * ci;
    R[0][1] = -co * sw - so * cw * ci;
    R[0][2] =  so * si;
    R[1][0] =  so * cw + co * sw * ci;
    R[1][1] = -so * sw + co * cw * ci;
    R[1][2] = -co * si;
    R[2][0] =  sw * si;
    R[2][1] =  cw * si;
    R[2][2] =  ci;

    for (int i = 0; i < 3; ++i) {
        r_out[i] = R[i][0] * r_pf[0] + R[i][1] * r_pf[1] + R[i][2] * r_pf[2];
        v_out[i] = R[i][0] * v_pf[0] + R[i][1] * v_pf[1] + R[i][2] * v_pf[2];
    }
}
