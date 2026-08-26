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
/*
 * SGP4/SDP4.
 *
 * This file exists so that the conformance harness
 * (tvb/validations/val_sgp4_vectors.c) compiles, runs, and reports how
 * far every published test case is from its reference value. Right now
 * that distance is "the whole orbit", which is the correct starting
 * measurement: the oracle is in place before the model, so every
 * equation added from here on has an immediate verdict.
 *
 * Order of work, from the equations in Spacetrack Report No. 3:
 *
 *   1. Recover the original mean motion and semimajor axis from the
 *      Brouwer elements (a1, delta1, a0, delta0, n0'', a0''). This is
 *      the Kozai correction, and it is the only piece that can be
 *      checked without the rest of the model working.
 *   2. Perigee-dependent drag setup: s and (q0-s)^4 modified below
 *      156 km and again below 98 km, and the simplified branch below
 *      220 km. Thresholds go into spody_const.h *with* this code, not
 *      before it -- a constant no code reads is a constant nobody has
 *      verified.
 *   3. Branch selection: period >= 225 min is deep space and a
 *      different model, not this one with extra terms.
 *   4. Near-Earth secular + periodic terms, then Kepler.
 *   5. Deep space: lunisolar terms, then the 24 h and 12 h resonances.
 *
 * Until (5) exists, a deep-space element set must fail loudly rather
 * than be propagated by the near-Earth branch: silently returning
 * near-Earth numbers for a GEO or a GNSS satellite is wrong by
 * thousands of kilometres while looking entirely healthy.
 */
#include <stddef.h>
#include <math.h>

#include "spody_sgp4.h"
#include "spody_const.h"

int spody_sgp4_init(const spody_sgp4_elements *el, spody_sgp4_sat *sat)
{
    /*
    * Reading the names. They are the ones Spacetrack Report No. 3 uses in
    * its FORTRAN listing (pp. 19-21), kept unchanged more less so that the code can
    * be read side by side with the report. FORTRAN IV typed any name starting 
    * with I-N as INTEGER, so a real variable called MDOT had to be spelled XMDOT. 
    * The leading X carries no meaning -- strip it. A trailing O means "original", the
    * value at epoch. SQ is squared, COF is "coefficient of".
    *
    *   theta  th2 th4   theta = cos i0, and its powers          p.13
    *   sinio                  sin i0
    *   betao  betao2          beta0 = sqrt(1-e0^2), and its square   p.14
    *   x3thm1                 3 theta^2 - 1
    *   x1mth2                 1 - theta^2
    *   x7thm1                 7 theta^2 - 1
    *   x1m5th                 1 - 5 theta^2
    *   xi                     xi = 1/(a0''-s): the reciprocal of how
    *                          far the orbit sits above the density
    *                          anchor s                               p.14
    *   eta  etasq  eeta       eta = a0'' e0 xi; eeta = e0 eta        p.14
    *   psisq                  psi^2 = |1 - eta^2|
    *   qoms24                 (q0 - s*)^4                            p.13
    *   coef   coef1           (q0-s)^4 xi^4, then times (1-eta^2)^-7/2
    *   pinvsq                 1/p^2, p = a0'' beta0^2 the semi-latus
    *                          rectum
    *   a3ovk2                 A(3,0)/k2 = -J3/k2
    */

    if (el == NULL || sat == NULL) {
        return SPODY_SGP4_ERR_ELEMENTS;
    }

    double theta  = cos(el->i0_deg * DEG2RAD);
    if (1.0 + theta < SGP4_MIN_1P_COSI) {
        return SPODY_SGP4_ERR_INCLINATION;
    }

    double theta2 = theta * theta;
    double x3thm1 = 3.0 * theta2 - 1.0; // "x3 theta minus 1"
    double betao2 = 1.0 - el->e0 * el->e0;
    double betao  = sqrt(betao2);
    double tval   = 1.5 * SGP4_WGS72_K2 * x3thm1 / (betao * betao2);

    double a1 = pow(SGP4_WGS72_KE/(el->n0_rev_day * REVDAY2RADMIN), 2.0/3.0);
    double d1 = tval / (a1 * a1);
    double a0 = a1 * (1.0 - d1 * (1.0/3.0 + d1 * (1.0 + 134.0/81.0 * d1)));
    double d0 = tval / (a0 * a0);

    sat->el = *el;
    sat->n0dp = el->n0_rev_day * REVDAY2RADMIN / (1.0 + d0);
    sat->a0dp = a0 / (1.0 - d0);

    sat->deep_space = ((2.0 * PI / sat->n0dp / MINUTESxDAY) >= SGP4_DEEP_SPACE_PERIOD_DAY);
    /* The deep-space half of the model does not exist yet. Refuse the
     * element set rather than let the near-Earth branch answer for it:
     * for a GEO or a GNSS satellite that answer would be wrong by
     * thousands of kilometres while looking perfectly healthy. */
    if (sat->deep_space) {
        return SPODY_SGP4_ERR_DEEP_SPACE;
    }

    sat->simple_drag = (sat->a0dp * (1.0 - el->e0) < SGP4_SIMPLE_DRAG_PERIGEE_KM / SGP4_WGS72_RE + 1.0);

    double perigee_km = (sat->a0dp * (1.0 - el->e0) - 1.0) * SGP4_WGS72_RE;
    double s_km = SGP4_S_KM;

    if (perigee_km < SGP4_S_MOD_PERIGEE_KM) {
        s_km = perigee_km - SGP4_S_KM;
        if (perigee_km <= SGP4_S_LOW_PERIGEE_KM) {
            s_km = SGP4_S_FLOOR_KM;
        }
    }
    sat->qoms24 = pow((SGP4_Q0_KM - s_km) / SGP4_WGS72_RE, 4.0);
    sat->s_er   = s_km / SGP4_WGS72_RE + 1.0;

    /* Geometry of the drag profile. xi is the reciprocal of how far
     * the semimajor axis sits above the density anchor s. */
    double xi     = 1.0 / (sat->a0dp - sat->s_er);
    sat->eta      = sat->a0dp * el->e0 * xi;
    double etasq  = sat->eta * sat->eta;
    double eeta   = el->e0 * sat->eta;

    /* fabs() because eta exceeds 1 on very eccentric orbits, and the
     * next line raises this to a fractional power: without it, five of
     * the published test cases produce NaN. */
    double psisq  = fabs(1.0 - etasq);

    double coef   = sat->qoms24 * xi * xi * xi * xi;
    double coef1  = coef / pow(psisq, 3.5);
    double a3ovk2 = -SGP4_WGS72_J3 / SGP4_WGS72_K2;   /* aE^3 = 1 */
    double theta4 = theta2 * theta2;
    double x1m5th = 1.0 - 5.0 * theta2; 
    double pinvsq = 1.0 / (sat->a0dp * sat->a0dp * betao2 * betao2);

    sat->theta  = theta;
    sat->sinio  = sin(el->i0_deg * DEG2RAD);
    sat->betao2 = betao2;
    sat->x1mth2 = 1.0 - theta2;
    sat->x3thm1 = x3thm1;
    sat->x7thm1 = 7.0 * theta2 - 1.0;

    double c2 = coef1 * sat->n0dp *
                (sat->a0dp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                 0.75 * SGP4_WGS72_K2 * xi / psisq * sat->x3thm1 *
                 (8.0 + 3.0 * etasq * (8.0 + etasq)));

    sat->c1 = el->bstar * c2;

    double c3 = coef * xi * a3ovk2 * sat->n0dp * sat->sinio / el->e0;
    if (el->e0 <= SGP4_LOW_ECC) { c3 = 0.0; }

    sat->c4 = 2.0 * sat->n0dp * coef1 * sat->a0dp * betao2 *
              (sat->eta * (2.0 + 0.5 * etasq) + el->e0 * (0.5 + 2.0 * etasq) -
               2.0 * SGP4_WGS72_K2 * xi / (sat->a0dp * psisq) *
               (-3.0 * sat->x3thm1 *
                    (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                0.75 * sat->x1mth2 *
                    (2.0 * etasq - eeta * (1.0 + etasq)) *
                    cos(2.0 * el->argp0_deg * DEG2RAD)));

    sat->c5 = 2.0 * coef1 * sat->a0dp * betao2 *
              (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);

    double temp1 = 3.0 * SGP4_WGS72_K2 * pinvsq * sat->n0dp;
    double temp2 = temp1 * SGP4_WGS72_K2 * pinvsq;
    double temp3 = 1.25 * SGP4_WGS72_K4 * pinvsq * pinvsq * sat->n0dp;

    sat->xmdot = sat->n0dp + 0.5 * temp1 * betao * sat->x3thm1 +
                 0.0625 * temp2 * betao *
                 (13.0 - 78.0 * theta2 + 137.0 * theta4);

    sat->omgdot = -0.5 * temp1 * x1m5th +
                  0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                  temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);

    double xhdot1 = -temp1 * theta;

    sat->xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * theta;
                            
    sat->omgcof = el->bstar * c3 * cos(el->argp0_deg * DEG2RAD);
    sat->xmcof  = -(2.0 / 3.0) * coef * el->bstar / eeta;
    if (el->e0 <= SGP4_LOW_ECC) { sat->xmcof = 0.0; }

    sat->xnodcf = 3.5 * betao2 * xhdot1 * sat->c1;
    sat->t2cof  = 1.5 * sat->c1;
    sat->xlcof  = 0.125 * a3ovk2 * sat->sinio * (3.0 + 5.0 * theta) / (1.0 + theta);             
    sat->aycof  = 0.25 * a3ovk2 * sat->sinio;
    sat->delmo  = pow(1.0 + sat->eta * cos(el->m0_deg * DEG2RAD), 3.0);
    sat->sinmo  = sin(el->m0_deg * DEG2RAD);
    
    if (sat->simple_drag) {
        sat->d2 = sat->d3 = sat->d4 = 0.0;
        sat->t3cof = sat->t4cof = sat->t5cof = 0.0;
    } else {
        double c1sq = sat->c1 * sat->c1;
        double tmp;

        sat->d2 = 4.0 * sat->a0dp * xi * c1sq;
            tmp = sat->d2 * xi * sat->c1 / 3.0;
        sat->d3 = (17.0 * sat->a0dp + sat->s_er) * tmp;
        sat->d4 = 0.5 * tmp * sat->a0dp * xi * (221.0 * sat->a0dp + 31.0 * sat->s_er) * sat->c1;
        sat->t3cof = sat->d2 + 2.0 * c1sq;
        sat->t4cof = 0.25 * (3.0 * sat->d3 + sat->c1 * (12.0 * sat->d2 + 10.0 * c1sq));  
        sat->t5cof = 0.2 * (3.0 * sat->d4 + 12.0 * sat->c1 * sat->d3 +
                            6.0 * sat->d2 * sat->d2 +
                            15.0 * c1sq * (2.0 * sat->d2 + c1sq));
    }

    sat->initialised = 1;
    return SPODY_SGP4_OK;
}

int spody_sgp4_at(const spody_sgp4_sat *sat, double tsince_min, double r_teme_km[3], double v_teme_kms[3])
{
                  
    if (sat == NULL || !sat->initialised ||
        r_teme_km == NULL || v_teme_kms == NULL) {
        return SPODY_SGP4_ERR_ELEMENTS;
    }

    const spody_sgp4_elements *el = &sat->el;
    double xmdf, omgadf, xnoddf, omega, xmp, tsq, xnode;
    double tempa, tempe, templ, a, e, xl, beta, xn;
    double axn, ayn, xlt, xll, aynl, capu, epw, temp;
    double sinepw = 0.0, cosepw = 0.0;
    double temp3 = 0.0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0;
    int it;

    /* ---- 1. secular: gravity, then drag ---- */
    xmdf   = el->m0_deg    * DEG2RAD + sat->xmdot  * tsince_min;
    omgadf = el->argp0_deg * DEG2RAD + sat->omgdot * tsince_min;
    xnoddf = el->raan0_deg * DEG2RAD + sat->xnodot * tsince_min;
    omega  = omgadf;
    xmp    = xmdf;
    tsq    = tsince_min * tsince_min;
    xnode  = xnoddf + sat->xnodcf * tsq;
    tempa  = 1.0 - sat->c1 * tsince_min;
    tempe  = el->bstar * sat->c4 * tsince_min;
    templ  = sat->t2cof * tsq;

    if (!sat->simple_drag) {
        double delomg = sat->omgcof * tsince_min;
        double delm   = sat->xmcof *
                        (pow(1.0 + sat->eta * cos(xmdf), 3.0) - sat->delmo);
        double tcube  = tsq * tsince_min;
        double tfour  = tsince_min * tcube;

        temp  = delomg + delm;
        xmp   = xmdf + temp;
        omega = omgadf - temp;
        tempa = tempa - sat->d2 * tsq - sat->d3 * tcube - sat->d4 * tfour;
        tempe = tempe + el->bstar * sat->c5 * (sin(xmp) - sat->sinmo);
        templ = templ + sat->t3cof * tcube +
                tfour * (sat->t4cof + tsince_min * sat->t5cof);
    }

    a    = sat->a0dp * tempa * tempa;
    e    = el->e0 - tempe;
    if (e < SGP4_ECC_TRAP)  { return SPODY_SGP4_ERR_PERT_ECC; }
    if (e < SGP4_ECC_FLOOR) { e = SGP4_ECC_FLOOR; }
    xl   = xmp + omega + xnode + sat->n0dp * templ;
    beta = sqrt(1.0 - e * e);
    xn   = SGP4_WGS72_KE / pow(a, 1.5);

    /* ---- 2. long-period periodics ---- */
    axn  = e * cos(omega);
    temp = 1.0 / (a * beta * beta);
    xll  = temp * sat->xlcof * axn;
    aynl = temp * sat->aycof;
    xlt  = xl + xll;
    ayn  = e * sin(omega) + aynl;

    /* ---- 3. Kepler, by Newton. FMOD2P inlined: two lines beat a
     * helper used once. The loop leaves sinepw..temp6 holding the
     * SECOND-to-last iterate, which is what ecose/esine want -- see
     * the listing on p.22 before tidying this up. ---- */
    capu = fmod(xlt - xnode, 2.0 * PI);
    if (capu < 0.0) {
        capu += 2.0 * PI;
    }
    epw = capu;
    for (it = 0; it < SGP4_KEPLER_MAX_ITER; it++) {
        double prev = epw;

        sinepw = sin(prev);
        cosepw = cos(prev);
        temp3  = axn * sinepw;
        temp4  = ayn * cosepw;
        temp5  = axn * cosepw;
        temp6  = ayn * sinepw;

        double step = (capu - temp4 + temp3 - prev) / (1.0 - temp5 - temp6);
        if (step >  SGP4_KEPLER_MAX_STEP) { step =  SGP4_KEPLER_MAX_STEP; }
        if (step < -SGP4_KEPLER_MAX_STEP) { step = -SGP4_KEPLER_MAX_STEP; }
        epw = prev + step;
        if (fabs(step) <= SGP4_KEPLER_TOL) { break; }
    }

    /* ---- 4. short-period preliminaries ---- */
    {
        double ecose = temp5 + temp6;
        double esine = temp3 - temp4;
        double elsq  = axn * axn + ayn * ayn;
        double tem   = 1.0 - elsq;
        double pl    = a * tem;
        double r     = a * (1.0 - ecose);
        double betal = sqrt(tem);
        double t3    = 1.0 / (1.0 + betal);
        double aor   = a / r;
        double rdot  = SGP4_WGS72_KE * sqrt(a)  * esine / r;
        double rfdot = SGP4_WGS72_KE * sqrt(pl) / r;
        double cosu  = aor * (cosepw - axn + ayn * esine * t3);
        double sinu  = aor * (sinepw - ayn - axn * esine * t3);
        /* ACTAN returns [0,2pi); only sin and cos of u are used
         * afterwards, so atan2's branch makes no difference. */
        double u     = atan2(sinu, cosu);
        double sin2u = 2.0 * sinu * cosu;
        double cos2u = 2.0 * cosu * cosu - 1.0;

        /* ---- 5. short-period periodics ---- */
        double tp1 = SGP4_WGS72_K2 / pl;
        double tp2 = tp1 / pl;
        double rk  = r * (1.0 - 1.5 * tp2 * betal * sat->x3thm1) +
                     0.5 * tp1 * sat->x1mth2 * cos2u;
        if (rk < 1.0) { return SPODY_SGP4_ERR_DECAYED; }

        double uk     = u - 0.25 * tp2 * sat->x7thm1 * sin2u;
        double xnodek = xnode + 1.5 * tp2 * sat->theta * sin2u;
        double xinck  = el->i0_deg * DEG2RAD +
                        1.5 * tp2 * sat->theta * sat->sinio * cos2u;
        double rdotk  = rdot - xn * tp1 * sat->x1mth2 * sin2u;
        double rfdotk = rfdot +
                        xn * tp1 * (sat->x1mth2 * cos2u + 1.5 * sat->x3thm1);

        /* ---- 6. orientation vectors, then state ---- */
        double sinuk  = sin(uk),     cosuk  = cos(uk);
        double sinik  = sin(xinck),  cosik  = cos(xinck);
        double sinnok = sin(xnodek), cosnok = cos(xnodek);
        double xmx = -sinnok * cosik;
        double xmy =  cosnok * cosik;
        double ux  = xmx * sinuk + cosnok * cosuk;
        double uy  = xmy * sinuk + sinnok * cosuk;
        double uz  = sinik * sinuk;
        double vx  = xmx * cosuk - cosnok * sinuk;
        double vy  = xmy * cosuk - sinnok * sinuk;
        double vz  = sinik * cosuk;

        /* SGP4 works in Earth radii and Earth radii per minute; the
         * report's driver converts on the way out (p.74). */
        double kmpv = SGP4_WGS72_RE * MINUTESxDAY / SECONDSxDAY;

        r_teme_km[0] = rk * ux * SGP4_WGS72_RE;
        r_teme_km[1] = rk * uy * SGP4_WGS72_RE;
        r_teme_km[2] = rk * uz * SGP4_WGS72_RE;
        v_teme_kms[0] = (rdotk * ux + rfdotk * vx) * kmpv;
        v_teme_kms[1] = (rdotk * uy + rfdotk * vy) * kmpv;
        v_teme_kms[2] = (rdotk * uz + rfdotk * vz) * kmpv;
    }
    return SPODY_SGP4_OK;

}
