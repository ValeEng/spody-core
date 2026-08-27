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
 * SGP4/SDP4: the analytic theory that GP element sets (TLE, OMM) are
 * fitted inside of. Their elements are mean elements of this theory,
 * so this file is not an approximation of their motion -- it is the
 * definition of it, and anything else fed the same numbers is simply
 * answering a different question.
 *
 * Equations from Hoots & Roehrich, Spacetrack Report No. 3 (1980).
 * Where the published verification vectors demand a departure from the
 * report, the departure is made and argued at the line that makes it,
 * with the reasoning in Vallado et al., AIAA 2006-6753 -- except twice,
 * where the vectors and that paper disagree and the vectors win: the
 * Greenwich angle at epoch, and the recovery of the semimajor axis.
 *
 * Two branches, chosen once at initialisation from the orbital period.
 * Below 225 minutes the near-Earth model; at or above it SDP4, which
 * adds the lunisolar secular and periodic terms and, inside the two
 * resonance bands, a Taylor integration of the mean longitude. It is a
 * different model rather than the same one with extra terms, which is
 * why the choice is made once and stored rather than re-tested.
 *
 * Everything runs in the report's units -- distance in Earth radii,
 * time in minutes -- so that ke alone carries the gravitational
 * parameter and the drag coefficients keep the magnitudes the report
 * prints. Kilometres appear once, in the last lines of spody_sgp4_at.
 *
 * Output is TEME, the frame the theory is defined in. Rotating it
 * elsewhere is left to the caller on purpose: a GP state carries the
 * accuracy of a fitted analytic theory, and quietly delivering it as
 * ICRF would hide which of the two the number came from.
 */
#include <stddef.h>
#include <math.h>

#include "spody_sgp4.h"
#include "spody_const.h"
#include "spody_earth_orientation.h"

/* Angle reduced to [0, 2pi). STR#3 calls this FMOD2P. */
static double fmod2p(double x)
{
    double r = fmod(x, 2.0 * PI);

    return (r < 0.0) ? r + 2.0 * PI : r;
}

/* Deep-space initialisation: the lunisolar geometry at epoch, the
 * secular rates it produces, and the resonance coefficients when the
 * mean motion falls in one of the two resonant bands. STR#3 calls this
 * DPINIT.*/
static void dpinit(spody_sgp4_sat *sat)
{
    const spody_sgp4_elements *el = &sat->el;
    spody_sgp4_deep *dp = &sat->deep;

    /* The report reaches DS50 by decoding the TLE card's two-digit year
     * and counting leap days; we are handed an MJD, so it is a
     * subtraction and the whole decode disappears. */
    double ds50 = el->epoch_mjd - SGP4_MJD_1950;
    double day  = ds50 + SGP4_DAY_1950_TO_1900;

    /* The TLE epoch goes in unadjusted: no dUT1. That is what produced
     * the element set, and STR#3's own linear expression in days since
     * 1950 -- 3.9e-04 deg away from sidereal time -- is the one thing
     * here the verification vectors will not tolerate, because this
     * angle cancels between xlamo and temp everywhere except inside
     * sin(xli - fasx2), where it sets the phase of the resonance
     * forcing and the integrator multiplies it. */
    dp->thgr = spody_gmst1982(el->epoch_mjd + JD_MJD_EPOCH);


    dp->xnq    = sat->n0dp;
    dp->xqncl  = el->i0_deg * DEG2RAD;
    dp->omegaq = el->argp0_deg * DEG2RAD;

    /* The Moon. Its node regresses around the ecliptic in 18.6 years,
     * which tilts the lunar orbit between 18.3 and 28.6 degrees from
     * the equator; the rest is that geometry referred to the equator. */
    double xnodce = SGP4_MOON_NODE_0 + SGP4_MOON_NODE_DOT * day;
    double stem   = sin(xnodce);
    double ctem   = cos(xnodce);
    double zcosil = SGP4_MOON_COSI_0 - SGP4_MOON_COSI_AMP * ctem;
    double zsinil = sqrt(1.0 - zcosil * zcosil);
    double zsinhl = SGP4_SIN_I_MOON * stem / zsinil;
    double zcoshl = sqrt(1.0 - zsinhl * zsinhl);

    double c   = SGP4_MOON_LON_0  + SGP4_MOON_LON_DOT  * day;
    double gam = SGP4_MOON_PERI_0 + SGP4_MOON_PERI_DOT * day;

    dp->zmol = fmod2p(c - gam);   /* mean longitude minus perigee */

    /* Argument of lunar perigee measured from the equatorial node. Only
     * its sine and cosine are wanted, so atan2's branch is irrelevant --
     * the same argument as on u in spody_sgp4_at. */
    double zx = atan2(SGP4_SIN_OBLIQ * stem / zsinil,
                      zcoshl * ctem + SGP4_COS_OBLIQ * zsinhl * stem);
    zx += gam - xnodce;

    double zcosgl = cos(zx);
    double zsingl = sin(zx);

    /* The Sun. */
    dp->zmos = fmod2p(SGP4_SUN_ANOM_0 + SGP4_SUN_ANOM_DOT * day);

    double xnoi   = 1.0 / dp->xnq;
    double siniq  = sat->sinio, cosiq = sat->theta;
    double eq     = el->e0, eqsq = eq * eq;
    double rteqsq = sqrt(sat->betao2), bsq = sat->betao2;
    double sinomo = sin(el->argp0_deg * DEG2RAD);
    double cosomo = cos(el->argp0_deg * DEG2RAD);
    double sinq   = sin(el->raan0_deg * DEG2RAD);
    double cosq   = cos(el->raan0_deg * DEG2RAD);
    int pass;

    for (pass = 0; pass < 2; pass++) {
        double zcosg, zsing, zcosi, zsini, zcosh, zsinh, cc, zn, ze;

        if (pass == 0) {          /* the Sun: its node IS the equinox */
            zcosg = SGP4_SUN_COSG;  zsing = SGP4_SUN_SING;
            zcosi = SGP4_COS_OBLIQ; zsini = SGP4_SIN_OBLIQ;
            zcosh = cosq;           zsinh = sinq;
            cc = SGP4_C1SS; zn = SGP4_ZNS; ze = SGP4_ZES;
        } else {                  /* the Moon */
            zcosg = zcosgl;         zsing = zsingl;
            zcosi = zcosil;         zsini = zsinil;
            zcosh = zcoshl * cosq + zsinhl * sinq;
            zsinh = sinq * zcoshl - cosq * zsinhl;
            cc = SGP4_C1L;  zn = SGP4_ZNL;  ze = SGP4_ZEL;
        }

        double a1  =  zcosg * zcosh + zsing * zcosi * zsinh;
        double a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
        double a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
        double a8  =  zsing * zsini;
        double a9  =  zsing * zsinh + zcosg * zcosi * zcosh;
        double a10 =  zcosg * zsini;
        double a2  =  cosiq * a7 + siniq * a8;
        double a4  =  cosiq * a9 + siniq * a10;
        double a5  = -siniq * a7 + cosiq * a8;
        double a6  = -siniq * a9 + cosiq * a10;

        double x1 =  a1 * cosomo + a2 * sinomo;
        double x2 =  a3 * cosomo + a4 * sinomo;
        double x3 = -a1 * sinomo + a2 * cosomo;
        double x4 = -a3 * sinomo + a4 * cosomo;
        double x5 =  a5 * sinomo;
        double x6 =  a6 * sinomo;
        double x7 =  a5 * cosomo;
        double x8 =  a6 * cosomo;

        double z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        double z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        double z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;

        double z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * eqsq;
        double z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * eqsq;
        double z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * eqsq;

        double z11 = -6.0 * a1 * a5 + eqsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        double z12 = -6.0 * (a1 * a6 + a3 * a5) + eqsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        double z13 = -6.0 * a3 * a6 + eqsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        double z21 =  6.0 * a2 * a5 + eqsq * ( 24.0 * x1 * x5 - 6.0 * x3 * x7);
        double z22 =  6.0 * (a4 * a5 + a2 * a6) + eqsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));           
        double z23 =  6.0 * a4 * a6 + eqsq * ( 24.0 * x2 * x6 - 6.0 * x4 * x8);

        z1 = z1 + z1 + bsq * z31;
        z2 = z2 + z2 + bsq * z32;
        z3 = z3 + z3 + bsq * z33;

        double s3 = cc * xnoi;
        double s2 = -0.5 * s3 / rteqsq;
        double s4 = s3 * rteqsq;
        double s1 = -15.0 * eq * s4;
        double s5 = x1 * x3 + x2 * x4;
        double s6 = x2 * x3 + x1 * x4;
        double s7 = x2 * x4 - x1 * x3;

        double se  =  s1 * zn * s5;
        double si  =  s2 * zn * (z11 + z13);
        double sl  = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eqsq);
        double sgh =  s4 * zn * (z31 + z33 - 6.0);
        double sh  = -zn * s2 * (z21 + z23);

        /* The node rate divides by sin i0 two lines down, and a
         * near-equatorial node carries nothing to divide. */
        if (dp->xqncl < SGP4_LOW_INCL) {
            sh = 0.0;
        }

        double ee2  =   2.0 * s1 * s6;
        double e3   =   2.0 * s1 * s7;
        double xi2  =   2.0 * s2 * z12;
        double xi3  =   2.0 * s2 * (z13 - z11);
        double xl2  =  -2.0 * s3 * z2;
        double xl3  =  -2.0 * s3 * (z3 - z1);
        double xl4  =  -2.0 * s3 * (-21.0 - 9.0 * eqsq) * ze;
        double xgh2 =   2.0 * s4 * z32;
        double xgh3 =   2.0 * s4 * (z33 - z31);
        double xgh4 = -18.0 * s4 * ze;
        double xh2  =  -2.0 * s2 * z22;
        double xh3  =  -2.0 * s2 * (z23 - z21);

        if (pass == 0) {
            dp->sse = se;  dp->ssi = si;  dp->ssl = sl;
            dp->ssh = sh / siniq;
            dp->ssg = sgh - cosiq * dp->ssh;

            dp->se2  = ee2;  dp->se3  = e3;
            dp->si2  = xi2;  dp->si3  = xi3;
            dp->sl2  = xl2;  dp->sl3  = xl3;   dp->sl4  = xl4;
            dp->sgh2 = xgh2; dp->sgh3 = xgh3;  dp->sgh4 = xgh4;
            dp->sh2  = xh2;  dp->sh3  = xh3;
        } else {
            dp->sse += se;  dp->ssi += si;  dp->ssl += sl;
            /* The report writes the solar term as cosiq*(sh/siniq) and
             * the lunar one as (cosiq/siniq)*sh. Same value, different
             * bracketing -- and the vectors carry both. */
            dp->ssg += sgh - cosiq / siniq * sh;
            dp->ssh += sh / siniq;

            dp->ee2  = ee2;  dp->e3   = e3;
            dp->xi2  = xi2;  dp->xi3  = xi3;
            dp->xl2  = xl2;  dp->xl3  = xl3;   dp->xl4  = xl4;
            dp->xgh2 = xgh2; dp->xgh3 = xgh3;  dp->xgh4 = xgh4;
            dp->xh2  = xh2;  dp->xh3  = xh3;
        }
    }

    /* The two resonance sets are mutually exclusive and neither is
     * written for a non-resonant element set. In FORTRAN the untouched
     * ones held the driver's zeros; in C they would hold the caller's
     * stack. */
    dp->resonant = dp->synchronous = 0;
    dp->del1 = dp->del2 = dp->del3 = 0.0;
    dp->fasx2 = dp->fasx4 = dp->fasx6 = 0.0;
    dp->d2201 = dp->d2211 = dp->d3210 = dp->d3222 = 0.0;
    dp->d4410 = dp->d4422 = dp->d5220 = dp->d5232 = 0.0;
    dp->d5421 = dp->d5433 = 0.0;
    dp->xlamo = dp->xfact = 0.0;

    double aqnv   = 1.0 / sat->a0dp;
    double xmao   = el->m0_deg   * DEG2RAD;
    double xnodeo = el->raan0_deg * DEG2RAD;
    double cosq2  = cosiq * cosiq;
    double bfact  = 0.0;

        if (dp->xnq < SGP4_RES_24H_HI && dp->xnq > SGP4_RES_24H_LO) {
        double g200 = 1.0 + eqsq * (-2.5 + 0.8125 * eqsq);
        double g310 = 1.0 + 2.0 * eqsq;
        double g300 = 1.0 + eqsq * (-6.0 + 6.60937 * eqsq);
        double f220 = 0.75 * (1.0 + cosiq) * (1.0 + cosiq);
        double f311 = 0.9375 * siniq * siniq * (1.0 + 3.0 * cosiq) - 0.75 * (1.0 + cosiq);
        double f330 = 1.0 + cosiq;
        double del1;

        f330 = 1.875 * f330 * f330 * f330;

        dp->resonant = 1;
        dp->synchronous = 1;

        /* del1 is a scratch amplitude first and the first coefficient
         * second -- the report reuses the name, and del2 and del3 read
         * the scratch value. Order is not negotiable here. */
        del1     = 3.0 * dp->xnq * dp->xnq * aqnv * aqnv;
        dp->del2 = 2.0 * del1 * f220 * g200 * SGP4_Q22;
        dp->del3 = 3.0 * del1 * f330 * g300 * SGP4_Q33 * aqnv;
        dp->del1 = del1 * f311 * g310 * SGP4_Q31 * aqnv;

        dp->fasx2 = SGP4_FASX2;
        dp->fasx4 = SGP4_FASX4;
        dp->fasx6 = SGP4_FASX6;

        dp->xlamo = xmao + xnodeo + dp->omegaq - dp->thgr;
        bfact = sat->xmdot + (sat->omgdot + sat->xnodot) - SGP4_THDT;
        bfact = bfact + dp->ssl + dp->ssg + dp->ssh;
    }
    else if (dp->xnq >= SGP4_RES_12H_LO && dp->xnq <= SGP4_RES_12H_HI &&
             eq >= SGP4_RES_12H_ECC) {
        double eoc  = eq * eqsq;
        double g201 = -0.306 - (eq - 0.64) * 0.440;
        double g211, g310, g322, g410, g422, g520, g521, g532, g533;
        double sini2, f220, f221, f321, f322, f441, f442, f522, f523, f542, f543;
        double xno2, ainv2, temp, temp1;

        dp->resonant = 1;

        if (eq <= SGP4_G_ECC_1) {
            g211 = 3.616 - 13.247 * eq + 16.290 * eqsq;
            g310 = -19.302 + 117.390 * eq - 228.419 * eqsq + 156.591 * eoc;
            g322 = -18.9068 + 109.7927 * eq - 214.6334 * eqsq + 146.5816 * eoc;
            g410 = -41.122 + 242.694 * eq - 471.094 * eqsq + 313.953 * eoc;
            g422 = -146.407 + 841.880 * eq - 1629.014 * eqsq + 1083.435 * eoc;
            /* The report prints -5740 here: the only coefficient of the
             * twenty-eight in this block written without decimals, and
             * the value the reference vectors are built on is
             * -5740.032. This branch is reached below e = 0.65, which
             * one published element set exercises and no other, so the
             * three missing digits are worth a factor 466 on it and
             * nothing anywhere else. */
            g520 = -532.114 + 3017.977 * eq - 5740.032 * eqsq + 3708.276 * eoc;
        } else {
            g211 = -72.099 + 331.819 * eq - 508.738 * eqsq + 266.724 * eoc;
            g310 = -346.844 + 1582.851 * eq - 2415.925 * eqsq + 1246.113 * eoc;
            g322 = -342.585 + 1554.908 * eq - 2366.899 * eqsq + 1215.972 * eoc;
            g410 = -1052.797 + 4758.686 * eq - 7193.992 * eqsq + 3651.957 * eoc;
            g422 = -3581.69 + 16178.11 * eq - 24462.77 * eqsq + 12422.52 * eoc;
            if (eq <= SGP4_G_ECC_2) {
                g520 = 1464.74 - 4664.75 * eq + 3763.64 * eqsq;
            } else {
                g520 = -5149.66 + 29936.92 * eq - 54087.36 * eqsq + 31324.56 * eoc;
            }
        }

        if (eq < SGP4_G_ECC_3) {
            g533 = -919.2277 + 4988.61 * eq - 9064.77 * eqsq + 5542.21 * eoc;
            g521 = -822.71072 + 4568.6173 * eq - 8491.4146 * eqsq + 5337.524 * eoc;
            g532 = -853.666 + 4690.25 * eq - 8624.77 * eqsq + 5341.4 * eoc;
        } else {
            g533 = -37995.78 + 161616.52 * eq - 229838.2 * eqsq + 109377.94 * eoc;
            g521 = -51752.104 + 218913.95 * eq - 309468.16 * eqsq + 146349.42 * eoc;
            g532 = -40023.88 + 170470.89 * eq - 242699.48 * eqsq + 115605.82 * eoc;
        }

        sini2 = siniq * siniq;
        f220  = 0.75 * (1.0 + 2.0 * cosiq + cosq2);
        f221  = 1.5 * sini2;
        f321  =  1.875 * siniq * (1.0 - 2.0 * cosiq - 3.0 * cosq2);
        f322  = -1.875 * siniq * (1.0 + 2.0 * cosiq - 3.0 * cosq2);
        f441  = 35.0 * sini2 * f220;
        f442  = 39.3750 * sini2 * sini2;
        f522  = 9.84375 * siniq * (sini2 * (1.0 - 2.0 * cosiq - 5.0 * cosq2) +
                                   0.33333333 * (-2.0 + 4.0 * cosiq + 6.0 * cosq2));
        f523  = siniq * (4.92187512 * sini2 * (-2.0 - 4.0 * cosiq + 10.0 * cosq2) +
                         6.56250012 * (1.0 + 2.0 * cosiq - 3.0 * cosq2));
        f542  = 29.53125 * siniq * (2.0 - 8.0 * cosiq +
                                    cosq2 * (-12.0 + 8.0 * cosiq + 10.0 * cosq2));
        f543  = 29.53125 * siniq * (-2.0 - 8.0 * cosiq +
                                    cosq2 * (12.0 + 8.0 * cosiq - 10.0 * cosq2));

        xno2  = dp->xnq * dp->xnq;
        ainv2 = aqnv * aqnv;

        temp1 = 3.0 * xno2 * ainv2;
        temp  = temp1 * SGP4_ROOT22;
        dp->d2201 = temp * f220 * g201;
        dp->d2211 = temp * f221 * g211;
        temp1 = temp1 * aqnv;
        temp  = temp1 * SGP4_ROOT32;
        dp->d3210 = temp * f321 * g310;
        dp->d3222 = temp * f322 * g322;
        temp1 = temp1 * aqnv;
        temp  = 2.0 * temp1 * SGP4_ROOT44;
        dp->d4410 = temp * f441 * g410;
        dp->d4422 = temp * f442 * g422;
        temp1 = temp1 * aqnv;
        temp  = temp1 * SGP4_ROOT52;
        dp->d5220 = temp * f522 * g520;
        dp->d5232 = temp * f523 * g532;
        temp  = 2.0 * temp1 * SGP4_ROOT54;
        dp->d5421 = temp * f542 * g521;
        dp->d5433 = temp * f543 * g533;

        dp->xlamo = xmao + xnodeo + xnodeo - dp->thgr - dp->thgr;
        bfact = sat->xmdot + sat->xnodot + sat->xnodot - SGP4_THDT - SGP4_THDT;
        bfact = bfact + dp->ssl + dp->ssh + dp->ssh;
    }

    if (dp->resonant) {
        dp->xfact = bfact - dp->xnq;
    }


}

/* Lunisolar secular perturbations, plus the resonance integration when
 * the element set is resonant. Rewrites the six mean elements in place.
 * STR#3 calls this DPSEC. */
static void dpsec(const spody_sgp4_sat *sat, double *xll, double *omgasm,
                  double *xnodes, double *em, double *xinc, double *xn,
                  double tsince_min)
{
    const spody_sgp4_deep *dp = &sat->deep;

    *xll    += dp->ssl * tsince_min;
    *omgasm += dp->ssg * tsince_min;
    *xnodes += dp->ssh * tsince_min;
    *em      = sat->el.e0 + dp->sse * tsince_min;
    *xinc    = dp->xqncl  + dp->ssi * tsince_min;

    if (!dp->resonant) {
        return;
    }

    /* The resonance is integrated from epoch on every call. STR#3 keeps
     * atime, xni and xli between calls and restarts only when crossing
     * the epoch: correct while the caller asks for monotonically
     * increasing times, wrong the moment it goes out and comes back,
     * because the answer then depends on the order the times were
     * asked in. Starting over costs one Taylor step per twelve hours
     * from epoch and buys a function of t alone. Which is also what
     * keeps this routine callable on one satellite from several
     * threads. */
    {
        double atime = 0.0;
        double xni   = dp->xnq;
        double xli   = dp->xlamo;
        double delt  = (tsince_min >= 0.0) ? SGP4_RES_STEP : -SGP4_RES_STEP;
        double xndot = 0.0, xnddt = 0.0, xldot = 0.0;
        double ft, xl, temp;

        for (;;) {
            if (dp->synchronous) {
                xndot = dp->del1 * sin(xli - dp->fasx2) +
                        dp->del2 * sin(2.0 * (xli - dp->fasx4)) +
                        dp->del3 * sin(3.0 * (xli - dp->fasx6));
                xnddt = dp->del1 * cos(xli - dp->fasx2) +
                        2.0 * dp->del2 * cos(2.0 * (xli - dp->fasx4)) +
                        3.0 * dp->del3 * cos(3.0 * (xli - dp->fasx6));
            } else {
                double xomi  = dp->omegaq + sat->omgdot * atime;
                double x2omi = xomi + xomi;
                double x2li  = xli + xli;

                xndot = dp->d2201 * sin(x2omi + xli - SGP4_G22) +
                        dp->d2211 * sin(xli - SGP4_G22) +
                        dp->d3210 * sin(xomi + xli - SGP4_G32) +
                        dp->d3222 * sin(-xomi + xli - SGP4_G32) +
                        dp->d4410 * sin(x2omi + x2li - SGP4_G44) +
                        dp->d4422 * sin(x2li - SGP4_G44) +
                        dp->d5220 * sin(xomi + xli - SGP4_G52) +
                        dp->d5232 * sin(-xomi + xli - SGP4_G52) +
                        dp->d5421 * sin(xomi + x2li - SGP4_G54) +
                        dp->d5433 * sin(-xomi + x2li - SGP4_G54);

                /* The four doubled terms are grouped at the end exactly
                 * as the report groups them: same value, and the sum
                 * lands on the same bits. */
                xnddt = dp->d2201 * cos(x2omi + xli - SGP4_G22) +
                        dp->d2211 * cos(xli - SGP4_G22) +
                        dp->d3210 * cos(xomi + xli - SGP4_G32) +
                        dp->d3222 * cos(-xomi + xli - SGP4_G32) +
                        dp->d5220 * cos(xomi + xli - SGP4_G52) +
                        dp->d5232 * cos(-xomi + xli - SGP4_G52) +
                        2.0 * (dp->d4410 * cos(x2omi + x2li - SGP4_G44) +
                               dp->d4422 * cos(x2li - SGP4_G44) +
                               dp->d5421 * cos(xomi + x2li - SGP4_G54) +
                               dp->d5433 * cos(-xomi + x2li - SGP4_G54));
            }

            /* xnddt comes out as d(xndot)/d(xli); the chain rule turns
             * it into a time derivative. */
            xldot = xni + dp->xfact;
            xnddt = xnddt * xldot;

            if (fabs(tsince_min - atime) < SGP4_RES_STEP) {
                break;
            }

            xli   += xldot * delt + xndot * SGP4_RES_STEP2;
            xni   += xndot * delt + xnddt * SGP4_RES_STEP2;
            atime += delt;
        }

        ft   = tsince_min - atime;
        *xn  = xni + xndot * ft + xnddt * ft * ft * 0.5;
        xl   = xli + xldot * ft + xndot * ft * ft * 0.5;
        temp = -*xnodes + dp->thgr + tsince_min * SGP4_THDT;

        *xll = dp->synchronous ? (xl - *omgasm + temp) : (xl + temp + temp);
    }
}

/* Lunisolar periodic perturbations, applied to the mean elements just
 * before Kepler's equation. STR#3 calls this DPPER. */
static void dpper(const spody_sgp4_sat *sat, double tsince_min, double *em,
                  double *xinc, double *omgasm, double *xnodes, double *xll)
{
    const spody_sgp4_deep *dp = &sat->deep;
    double zm, zf, sinzf, f2, f3;

    /* STR#3 recomputes these only when the propagation time has moved by
     * more than 30 minutes, caching them across calls in SAVTSN. That
     * cache is dropped: it is state between calls, which this module
     * does not have by design, and the paper dropped it too -> "thus
     * resulting in smoother behavior for deep-space orbits with small
     * time steps". */
    zm    = dp->zmos + SGP4_ZNS * tsince_min;
    zf    = zm + 2.0 * SGP4_ZES * sin(zm);   /* equation of the centre */
    sinzf = sin(zf);
    f2    =  0.5 * sinzf * sinzf - 0.25;
    f3    = -0.5 * sinzf * cos(zf);

    double ses  = dp->se2  * f2 + dp->se3  * f3;
    double sis  = dp->si2  * f2 + dp->si3  * f3;
    double sls  = dp->sl2  * f2 + dp->sl3  * f3 + dp->sl4  * sinzf;
    double sghs = dp->sgh2 * f2 + dp->sgh3 * f3 + dp->sgh4 * sinzf;
    double shs  = dp->sh2  * f2 + dp->sh3  * f3;

    zm    = dp->zmol + SGP4_ZNL * tsince_min;
    zf    = zm + 2.0 * SGP4_ZEL * sin(zm);
    sinzf = sin(zf);
    f2    =  0.5 * sinzf * sinzf - 0.25;
    f3    = -0.5 * sinzf * cos(zf);

    double sel  = dp->ee2  * f2 + dp->e3   * f3;
    double sil  = dp->xi2  * f2 + dp->xi3  * f3;
    double sll  = dp->xl2  * f2 + dp->xl3  * f3 + dp->xl4  * sinzf;
    double sghl = dp->xgh2 * f2 + dp->xgh3 * f3 + dp->xgh4 * sinzf;
    double shl  = dp->xh2  * f2 + dp->xh3  * f3;

    double pe   = ses  + sel;
    double pinc = sis  + sil;
    double pl   = sls  + sll;
    double pgh  = sghs + sghl;
    double ph   = shs  + shl;

    *xinc += pinc;
    *em   += pe;

        /* Tested on the perturbed inclination, not the one at epoch. STR#3
     * used the epoch value; AIAA 2006-6753 p.10 says the code shipped
     * with the paper "included Option (b)", which is testing at each
     * propagation time. The vectors were made by the code, not by the
     * recommendation. */
    if (*xinc >= SGP4_LYDDANE_INCL) {
        /* The perturbed inclination here too, not the epoch SINIQ and
         * COSIQ of the report: the same rule as the branch below, and
         * the same source. Every term based on the Keplerian orbit is
         * rebuilt from the perturbed values. Worth four orders of
         * magnitude on the inclined cases. */
        double sinip = sin(*xinc), cosip = cos(*xinc);

        ph   = ph / sinip;
        pgh  = pgh - cosip * ph;

        *omgasm += pgh;
        *xnodes += ph;
        *xll    += pl;
    } else {
        /* The node is reduced here, not merely used. Two lines below it
         * enters xls and dls as a NUMBER and not as an angle, times
         * cos i, and times pinc sin i, so whether it reads -6e-06 or
         * 6.28318 moves the answer by pinc * 2pi * sin i. Reducing it
         * also cancels out of xlt - xnode further downstream, so
         * nothing else notices. Case 23599 is exact to nanometres up to
         * t=400 and a kilometre out from t=420, which is exactly where
         * its node crosses zero. */
        *xnodes = fmod2p(*xnodes);
        /* Recomputed here rather than at entry: the report takes sin i
         * and cos i before the periodic correction reaches the
         * inclination, and the corrected implementations rebuild every
         * term based on the Keplerian orbit from the perturbed values.
         * This one line is worth six orders of magnitude on the
         * low-inclination cases -> 9998 goes from 8.2 km to 4 mm. */
        double sinis = sin(*xinc), cosis = cos(*xinc);
        double sinok = sin(*xnodes), cosok = cos(*xnodes);

        double alfdp = sinis * sinok;
        double betdp = sinis * cosok;
        double dalf  =  ph * cosok + pinc * cosis * sinok;
        double dbet  = -ph * sinok + pinc * cosis * cosok;
        double xls, dls;

        alfdp += dalf;
        betdp += dbet;

        xls = *xll + *omgasm + cosis * *xnodes;
        dls = pl + pgh - pinc * *xnodes * sinis;
        xls = xls + dls;

        /* ACTAN, and the reduction to [0,2pi) is load-bearing here.
         * Elsewhere the report's ACTAN can be replaced by a bare atan2,
         * because only the sine and cosine of the result are ever used
         * and a turn either way cancels. Not here: xnodes goes straight
         * into a linear combination two lines down, so the turn it
         * lands on shifts omgasm by 2 pi cos i. Cases 23177 and 23599
         * exist to catch exactly this, one for each way of missing
         * it. */
        {
            double xnoh = *xnodes;

            *xnodes = fmod2p(atan2(alfdp, betdp));

            /* The reduction to [0,2pi) can land the node a whole turn
             * from where it was an instant earlier, and omgasm is built
             * by subtracting it two lines down. The fix is to keep the
             * new node on the same turn as the old one. Case 23599 is
             * 1122 km wrong without this and 0.96 km with it. */
            if (fabs(xnoh - *xnodes) > PI) {
                if (*xnodes < xnoh) {
                    *xnodes += 2.0 * PI;
                } else {
                    *xnodes -= 2.0 * PI;
                }
            }
        }

        *xll   += pl;
        *omgasm = xls - *xll - cos(*xinc) * *xnodes;
    }
}


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
    
    /* The recovered pair must satisfy Kepler's third law exactly:
     * n0''^2 a0''^3 = ke^2. STR#3 instead reaches a0'' by a series in
     * d truncated after the third order, so the residual grows as d^4:
     * 1.7e-12 at the d = -7e-4 of an ordinary element set, 7.1e-4 at
     * the d = -0.108 that e = 0.995 produces, because d carries
     * (1 - e^2)^1.5 in its denominator. It is eccentricity that breaks
     * the series, not altitude. */
    sat->a0dp = pow(SGP4_WGS72_KE / sat->n0dp, 2.0 / 3.0);  

    sat->deep_space = ((2.0 * PI / sat->n0dp / MINUTESxDAY) >= SGP4_DEEP_SPACE_PERIOD_DAY);

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

    sat->c4 = 2.0 * sat->n0dp * coef1 * sat->a0dp * betao2 *
              (sat->eta * (2.0 + 0.5 * etasq) + el->e0 * (0.5 + 2.0 * etasq) -
               2.0 * SGP4_WGS72_K2 * xi / (sat->a0dp * psisq) *
               (-3.0 * sat->x3thm1 *
                    (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                0.75 * sat->x1mth2 *
                    (2.0 * etasq - eeta * (1.0 + etasq)) *
                    cos(2.0 * el->argp0_deg * DEG2RAD)));

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

    sat->xnodcf = 3.5 * betao2 * xhdot1 * sat->c1;
    sat->t2cof  = 1.5 * sat->c1;
    sat->xlcof  = 0.125 * a3ovk2 * sat->sinio * (3.0 + 5.0 * theta) / (1.0 + theta);             
    sat->aycof  = 0.25 * a3ovk2 * sat->sinio;
    
    if (sat->deep_space) {
        /* None of the near-Earth drag machinery exists in this branch:
         * SDP4 has no simplified-drag case, no C3, no C5 and no powers
         * of time past the second. Zeroing what the propagator still
         * reads is not tidiness -- in C those fields would otherwise be
         * whatever the caller's stack held. Same argument as the GO TO
         * 90 of the report, where the zeros came from the driver. */
        sat->simple_drag = 0;
        sat->c5 = sat->omgcof = sat->xmcof = 0.0;
        sat->delmo = sat->sinmo = 0.0;
        sat->d2 = sat->d3 = sat->d4 = 0.0;
        sat->t3cof = sat->t4cof = sat->t5cof = 0.0;
        dpinit(sat);
    } else {
        double c3 = coef * xi * a3ovk2 * sat->n0dp * sat->sinio / el->e0;
        if (el->e0 <= SGP4_LOW_ECC) { c3 = 0.0; }

        sat->simple_drag = (sat->a0dp * (1.0 - el->e0) < SGP4_SIMPLE_DRAG_PERIGEE_KM / SGP4_WGS72_RE + 1.0);

        sat->c5 = 2.0 * coef1 * sat->a0dp * betao2 *
                  (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
        sat->omgcof = el->bstar * c3 * cos(el->argp0_deg * DEG2RAD);
        sat->xmcof  = -(2.0 / 3.0) * coef * el->bstar / eeta;
        if (el->e0 <= SGP4_LOW_ECC) { sat->xmcof = 0.0; }
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
    double tempa, tempe, templ, xinc, a, e, xl, beta, xn;
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

    xinc = el->i0_deg * DEG2RAD;

    if (sat->deep_space) {
        xn = sat->n0dp;
        e  = el->e0;
        dpsec(sat, &xmp, &omega, &xnode, &e, &xinc, &xn, tsince_min);

        /* The resonance moves the mean motion, so the semimajor axis
         * follows it instead of starting from its value at epoch. This
         * is the one line SDP4 replaces rather than adds to. */
        a = pow(SGP4_WGS72_KE / xn, 2.0 / 3.0) * tempa * tempa;
        e = e - tempe;
        if (e < SGP4_ECC_TRAP)  { return SPODY_SGP4_ERR_PERT_ECC; }
        if (e < SGP4_ECC_FLOOR) { e = SGP4_ECC_FLOOR; }

        xmp = xmp + sat->n0dp * templ;
        dpper(sat, tsince_min, &e, &xinc, &omega, &xnode, &xmp);
        
        /* The periodics move the eccentricity too, and on a badly
         * conditioned element set they move it a long way: case 33334,
         * whose mean motion is 1e-5 rev/day, comes out of dpper with
         * e = -122. Checking only before them leaves sqrt(1 - e*e) to
         * return a NaN and the propagator to report it as a position. */
        if (e < SGP4_ECC_TRAP)  { return SPODY_SGP4_ERR_PERT_ECC; }
        if (e < SGP4_ECC_FLOOR) { e = SGP4_ECC_FLOOR; }

        /* A mean inclination the Sun and Moon have pushed below zero is
         * the same orbit reflected: flip it and turn both angles through
         * pi. STR#3 does this inside the secular routine, which is too
         * early -- the periodics can carry the inclination back across
         * zero, and correcting before they run leaves a step in the z
         * component. AIAA 2006-6753 p.15 moves it here; 25954 past 274
         * min, 28626 past 1130 min and 26900 at 9313 min are the three
         * cases built to show the difference. */
        if (xinc < 0.0) {
            xinc   = -xinc;
            xnode += PI;
            omega -= PI;
        }
        xl  = xmp + omega + xnode;
    } else {
        if (!sat->simple_drag) {
            double delomg = sat->omgcof * tsince_min;
            double delm   = sat->xmcof * (pow(1.0 + sat->eta * cos(xmdf), 3.0) - sat->delmo);                         
            double tcube  = tsq * tsince_min;
            double tfour  = tsince_min * tcube;

            temp  = delomg + delm;
            xmp   = xmdf + temp;
            omega = omgadf - temp;
            tempa = tempa - sat->d2 * tsq - sat->d3 * tcube - sat->d4 * tfour;
            tempe = tempe + el->bstar * sat->c5 * (sin(xmp) - sat->sinmo);
            templ = templ + sat->t3cof * tcube + tfour * (sat->t4cof + tsince_min * sat->t5cof);   
        }

        a = sat->a0dp * tempa * tempa;
        e = el->e0 - tempe;
        if (e < SGP4_ECC_TRAP)  { return SPODY_SGP4_ERR_PERT_ECC; }
        if (e < SGP4_ECC_FLOOR) { e = SGP4_ECC_FLOOR; }

        xl = xmp + omega + xnode + sat->n0dp * templ;
    }

    /* Initialisation computed these from the inclination at epoch. The
     * deep-space terms have since moved it, and the corrected
     * implementations recompute any term based on the Keplerian orbit
     * from the perturbed values (AIAA 2006-6753 p.9). In the near-Earth
     * branch xinc is still i0, so these stay copies. */
    double theta  = sat->theta;
    double sinio  = sat->sinio;
    double x1mth2 = sat->x1mth2;
    double x3thm1 = sat->x3thm1;
    double x7thm1 = sat->x7thm1;
    double xlcof  = sat->xlcof;
    double aycof  = sat->aycof;

    if (sat->deep_space) {
        double a3ovk2 = -SGP4_WGS72_J3 / SGP4_WGS72_K2;   /* aE^3 = 1 */
        double t2;

        theta  = cos(xinc);
        sinio  = sin(xinc);
        t2     = theta * theta;
        x1mth2 = 1.0 - t2;
        x3thm1 = 3.0 * t2 - 1.0;
        x7thm1 = 7.0 * t2 - 1.0;
        xlcof  = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * theta) / (1.0 + theta);
        aycof  = 0.25 * a3ovk2 * sinio;
    }

    beta = sqrt(1.0 - e * e);
    xn   = SGP4_WGS72_KE / pow(a, 1.5);

    /* ---- 2. long-period periodics ---- */
    axn  = e * cos(omega);
    temp = 1.0 / (a * beta * beta);
    xll  = temp * xlcof * axn;
    aynl = temp * aycof;
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
        /* elsq is the squared magnitude of the eccentricity vector after
         * the periodics, and nothing has kept it below one: case 33333
         * crosses at t=25 and the orbit stops being an orbit. */
        if (pl < 0.0) {
            return SPODY_SGP4_ERR_SEMI_LATUS;
        }

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
        double rk  = r * (1.0 - 1.5 * tp2 * betal * x3thm1) +
                     0.5 * tp1 * x1mth2 * cos2u;
        if (rk < 1.0) { return SPODY_SGP4_ERR_DECAYED; }

        double uk     = u - 0.25 * tp2 * x7thm1 * sin2u;
        double xnodek = xnode + 1.5 * tp2 * theta * sin2u;
        double xinck  = xinc +  1.5 * tp2 * theta * sinio * cos2u;                      
        double rdotk  = rdot - xn * tp1 * x1mth2 * sin2u;
        double rfdotk = rfdot +  xn * tp1 * (x1mth2 * cos2u + 1.5 * x3thm1);
                       
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