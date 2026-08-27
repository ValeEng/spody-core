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
 * SGP4/SDP4 -- the analytic propagator that GP mean element sets are
 * fitted inside of.
 *
 * A GP element set (a TLE, or the OMM that supersedes it) is NOT an
 * osculating state. Its elements are *mean* elements of a specific
 * analytic theory, and they only mean anything when evaluated by that
 * theory. Handing them to a numerical integrator is a physical error,
 * not an approximation. This module is the theory: give it an element
 * set and a time offset, get position and velocity in TEME.
 *
 * Sources. The equations are those of Spacetrack Report No. 3 (Hoots &
 * Roehrich, USAF, 1980), whose symbols are defined in the companion
 * Spacetrack Report No. 2 (Lane & Hoots); the corrections applied to
 * them since are documented in Vallado, Crawford, Hujsak & Kelso,
 * "Revisiting Spacetrack Report #3", AIAA 2006-6753. The model itself
 * descends from Brouwer (1959) for the geopotential and Lane &
 * Cranford (1969) for drag, with the deep-space extension of Hujsak
 * (1979). Conformance is verified against the test cases published
 * with AIAA 2006-6753 -- see tvb/validations/val_sgp4_vectors.c.
 *
 * Constants live in spody_const.h under SGP4_*, and are WGS-72. They
 * are part of the model's definition, not measurements of the Earth:
 * see the comment on that block before being tempted to improve them.
 *
 * Frame and units. Output is TEME (True Equator, Mean Equinox), which
 * is the frame SGP4 works in and is *not* ICRF; rotating to an
 * inertial frame is a separate step, deliberately kept out of this
 * interface so that conformance tests exercise the model alone.
 * Positions are km, velocities km/s, times minutes from the element
 * set's epoch, signed -- propagating backwards is part of the
 * published contract, and so is propagating years away from epoch.
 */

#ifndef SPODY_SGP4_H
#define SPODY_SGP4_H

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes. Zero is success; the non-zero values follow the
 * numbering used by the reference implementations, so that a failure
 * can be compared against published behaviour. Three of the published
 * test cases exist specifically to provoke these. */
enum {
    SPODY_SGP4_OK              = 0,
    SPODY_SGP4_ERR_ELEMENTS    = 1, /* eccentricity or a out of range  */
    SPODY_SGP4_ERR_MEAN_MOTION = 2, /* mean motion went negative       */
    SPODY_SGP4_ERR_PERT_ECC    = 3, /* perturbed eccentricity out of   */
                                    /* [0,1)                           */
    SPODY_SGP4_ERR_SEMI_LATUS  = 4, /* semi-latus rectum went negative */
    SPODY_SGP4_ERR_SUBORBITAL  = 5, /* epoch elements are sub-orbital  */
    SPODY_SGP4_ERR_DECAYED     = 6, /* satellite has decayed           */
    /* Ours, not the report's: numbered after its codes so the two
     * cannot be confused. Seven was the deep-space placeholder and is
     * gone now that SDP4 is here; eight keeps its number, so a code
     * that has already shipped still means what it meant. */
    SPODY_SGP4_ERR_INCLINATION = 8  /* i0 near 180 deg: 1 + cos i0 will*/ 
                                    /* be 0                            */

};

/* Mean elements exactly as a GP source delivers them: the seven
 * orbital quantities plus the drag term. Deliberately the OMM's own
 * fields in the OMM's own units, so that feeding this from a catalogue
 * record is a copy and not a conversion.
 *
 * epoch_mjd is needed only by the deep-space branch, which has to
 * place the Sun and the Moon; the near-Earth branch ignores it. It is
 * still required, because whether a set is deep space is decided by
 * the elements, not by the caller.
 *
 * ndot and nddot are carried by the format but are not used by SGP4:
 * drag enters through bstar alone. They are absent here on purpose --
 * a field that no equation reads is a field someone will eventually
 * try to make matter. */
typedef struct {
    double epoch_mjd;    /* UTC MJD of the element set                */
    double n0_rev_day;   /* Brouwer mean motion [rev/day]             */
    double e0;           /* eccentricity [-]                          */
    double i0_deg;       /* inclination [deg]                         */
    double raan0_deg;    /* right ascension of ascending node [deg]   */
    double argp0_deg;    /* argument of perigee [deg]                 */
    double m0_deg;       /* mean anomaly [deg]                        */
    double bstar;        /* drag term [1/(Earth radii)]               */
} spody_sgp4_elements;

/* Everything the deep-space branch precomputes once per element set.
 * Nested rather than appended so a near-Earth reader can skip it in one
 * jump, and always present rather than allocated: it costs about 350
 * bytes on an element set that never uses it, against a pointer, a
 * lifetime and a failure mode it would otherwise cost on every one.
 * Filled by the deep-space initialisation, read by the two deep-space
 * perturbation steps; nothing outside this module touches it. */
typedef struct {
    double thgr;                 /* Greenwich sidereal angle at epoch  */
    double xnq, xqncl, omegaq;   /* n0'', i0, w0 frozen at epoch       */

    double ssl, ssg, ssh, sse, ssi;          /* lunisolar secular      */
    double se2, se3, si2, si3, sl2, sl3, sl4;      /* lunisolar        */
    double sgh2, sgh3, sgh4, sh2, sh3;             /* periodic         */
    /* The same twelve quantities for the Moon. The report leaves them
     * in the working names because the lunar pass runs second; naming
     * them apart is the difference between reading DPPER and guessing
     * at it. */
    double ee2, e3, xi2, xi3, xl2, xl3, xl4;
    double xgh2, xgh3, xgh4, xh2, xh3;
    double zmos, zmol;           /* solar and lunar mean anomalies     */

    int    resonant;             /* IRESFL                             */
    int    synchronous;          /* ISYNFL: the 24 h kind              */
    double del1, del2, del3;                 /* 24 h amplitudes        */
    double fasx2, fasx4, fasx6;              /* 24 h phases            */
    double d2201, d2211, d3210, d3222;       /* 12 h                   */
    double d4410, d4422, d5220, d5232;
    double d5421, d5433;
    double xlamo, xfact;         /* integration start and frequency    */
} spody_sgp4_deep;

/* Everything the model precomputes once per element set.
 *
 * Callers set nothing here and read nothing here: it is filled by
 * spody_sgp4_init and consumed by spody_sgp4_at. It is spelled out
 * rather than hidden behind an opaque handle because the whole library
 * is built from source in one tree -- there is no ABI to protect, and
 * a concrete type can be put on the stack or held in a flat array of
 * tens of thousands, which is what conjunction screening will want.
 *
 * It grows as the branches land: the drag coefficients, then the
 * deep-space workspace. Growing it costs a recompile and nothing
 * else. */
typedef struct {
    /* ---- as supplied ---- */
    spody_sgp4_elements el;
    spody_sgp4_deep deep;
    /* ---- implementation-owned; do not read from outside ---- */
    double n0dp;            /* recovered mean motion [rad/min]          | dp = double prime */
    double a0dp;            /* recovered semimajor axis [Earth radii]   | dp = double prime */
    double s_er;            /* s in Earth radii                         | er = Earth radii  */
    double qoms24;          /* (q0-s)^4 in Earth radii^4                                    */
    double theta;           /* cos i0                                       */
    double eta;             /* a0'' * e0 * xi                               */
    double sinio;           /* sin i0                                       */
    double betao2;          /* 1 - e0^2                                     */
    double x1mth2;          /* 1 - theta^2                                  */
    double x3thm1;          /* 3 theta^2 - 1                                */
    double x7thm1;          /* 7 theta^2 - 1                                */
    double c1, c4, c5;      /* drag coefficients read at propagation time   */
    double xmdot, omgdot, xnodot;      /* secular rates                */
    double omgcof, xmcof, xnodcf;      /* drag corrections to omega, M */
    double t2cof, t3cof, t4cof, t5cof; /* powers of time since epoch   */
    double xlcof, aycof;               /* long-period periodics        */
    double delmo, sinmo;               /* frozen at epoch              */
    double d2, d3, d4;                 /* higher-order drag            */


    int    deep_space;      /* branch tag: period >= 225 min                                */
    int    simple_drag;     /* branch tag: perigee < 220 km                                 */
    int    initialised;

} spody_sgp4_sat;

/* Initialise *sat* from *el*. Returns SPODY_SGP4_OK, or a non-zero
 * code when the element set cannot be propagated at all (a sub-orbital
 * perigee, for instance). On failure *sat* is left unusable. */
int spody_sgp4_init(const spody_sgp4_elements *el, spody_sgp4_sat *sat);

/* Evaluate at *tsince_min* minutes from epoch, signed. Writes TEME
 * position [km] and velocity [km/s] and returns SPODY_SGP4_OK, or
 * leaves them untouched and returns a non-zero code when the model has
 * left its domain of validity at that time -- which is a legitimate
 * answer, not a bug: a decayed satellite has no position.
 *
 * *sat* is const: the initialised satellite is immutable, and error
 * state is the return value rather than something cached inside. One
 * satellite may therefore be evaluated from several threads at once. */
int spody_sgp4_at(const spody_sgp4_sat *sat, double tsince_min,
                  double r_teme_km[3], double v_teme_kms[3]);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_SGP4_H */
