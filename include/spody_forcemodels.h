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
#ifndef SPODY_FORCEMODELS_H
#define SPODY_FORCEMODELS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spody_const.h"
#include "spody_eclipse.h"
#include "spody_ephemeris.h"
#include "spody_harmonics.h"
#include "spody_integrators.h"   /* spody_rhs_fn */

/* Maximum number of third bodies that fit in a ForceBreakdown. The
 * ForceModelContext itself accepts an unbounded list, but the diagnostic
 * struct stores per-body accelerations in a fixed-size array to keep
 * itself trivially copyable / serializable. Real missions rarely exceed 3-5. */
#define SPODY_FM_MAX_THIRD 8

/* ============================================================
 * Spacecraft parameters
 * ============================================================
 * Mass and surface properties used by drag and SRP. The two
 * area/mass ratios are pre-computed by spody_init_Spacecraft so
 * the per-step force functions can avoid a division in the hot
 * path. Fill mass, area_*, Cd, Cr; then call the init.
 */
typedef struct {
    double mass;        /* kg                                          */
    double area_drag;   /* m^2 - cross section for drag                */
    double area_srp;    /* m^2 - cross section for SRP                 */
    double Cd;          /* drag coefficient                            */
    double Cr;          /* SRP reflectivity (1=absorb, 2=mirror)       */

    /* derived (do not set by hand; populated by spody_init_Spacecraft) */
    double am_drag;     /* = area_drag / mass  [m^2/kg]                */
    double am_srp;      /* = area_srp  / mass  [m^2/kg]                */
} Spacecraft;

/* Populate am_drag and am_srp as area / mass. Call once after filling
 * mass and the area_* fields, before passing the Spacecraft to a
 * force-model context. Safe to call again if you change mass/area. */
void spody_init_Spacecraft(Spacecraft *sc);

/* ============================================================
 * Force-model context
 * ============================================================
 * Aggregates everything the default RHS needs to evaluate dy/dt.
 * Built once by the caller (typically on the main thread) and then
 * passed as the opaque `user` pointer to the integrator.
 *
 * Threading:
 *   - The Spacecraft, third-body arrays, and scalar fields are
 *     read-only after setup -- safe to share across threads.
 *   - HarmonicGravity and MappedEphemeris are per-thread handles.
 *     Each worker thread must have its own ForceModelContext (or
 *     a copy of the shared fields) carrying its own hg/eph handles.
 *
 * Force toggling is purely runtime: a NULL pointer (hg, eph) or a
 * zero counter (n_third) or zero flag (enable_*) disables the force.
 * Branch overhead measured on the bench was within noise (<1% at
 * harmonics N>=100), so this is preferred over compile-time #if.
 */
/* Body-fixed orientation provider.
 *
 * Returns the rotation from ICRF to the central body's body-fixed
 * frame (and its inverse) at the given Ephemeris Time. Used by
 * `spody_force_sphericalharmonics` to rotate the satellite state
 * into the frame the harmonic coefficients are expressed in,
 * evaluate the field, and rotate the resulting acceleration back
 * to ICRF.
 *
 * Each known central body has its own implementation that "knows"
 * the body it serves (no need for a separate context pointer):
 *   - Moon  -> spody_bf_rotation_moon  (DE440 libration in slot 12)
 *   - Earth -> [TODO] GMST/IAU 2006 based
 *   - Mars  -> [TODO] IAU 2009 based
 *
 * `eph` is the only runtime input. Bodies whose orientation lives
 * in the ephemeris file (Moon librations) read it here; bodies
 * whose orientation comes from an analytic model independent of
 * the ephemeris (Earth via GMST) simply ignore the argument. */
typedef void (*spody_bf_rotation_fn)(MappedEphemeris *eph, double et,
                                      double R_icrf_to_bf[3][3],
                                      double R_bf_to_icrf[3][3]);

/* Concrete provider: Moon Principal Axes from DE440 libration
 * angles. Equivalent to chaining
 *   spody_get_lunarlibrationangles(eph, et, angles);
 *   spody_getrotmatrix_icrf2moonpa(angles..., R_icrf_to_bf);
 *   spody_getrotmatrix_moonpa2icrf(angles..., R_bf_to_icrf);
 * but exposed as a function pointer so the force-model layer can
 * call it generically. */
void spody_bf_rotation_moon(MappedEphemeris *eph, double et,
                             double R_icrf_to_bf[3][3],
                             double R_bf_to_icrf[3][3]);

typedef struct {
    /* central body (the body the satellite orbits) */
    double  mu_central;          /* km^3/s^2                          */
    double  R_central;           /* km - mean radius (impact check)   */
    int     naif_central;        /* NAIF id, e.g. 301 (Moon), 399 (Earth) */

    /* Body-fixed orientation provider for `naif_central`. Used by
     * the spherical-harmonics force to rotate state ICRF <-> body
     * fixed frame at every RHS evaluation. Application fills this
     * based on the parsed central-body name; for Moon use
     * `spody_bf_rotation_moon`. MUST be non-NULL when `hg` is set. */
    spody_bf_rotation_fn get_bf_rotation;

    /* spacecraft */
    const Spacecraft *sat;

    /* spherical harmonics on the central body (NULL = disabled) */
    HarmonicGravity *hg;

    /* ephemeris-driven perturbations (NULL = disabled). Must be
     * non-NULL whenever hg, n_third > 0, or enable_srp are active. */
    MappedEphemeris *eph;

    /* third bodies: parallel arrays of NAIF ids and GMs.
     * n_third == 0 disables third-body perturbations. */
    const int    *third_naif;
    const double *third_mu;      /* km^3/s^2 per body */
    int           n_third;

    /* solar radiation pressure (cannonball + cylindrical eclipse) */
    int     enable_srp;
    int     srp_occulter_naif;       /* body that can shadow the sat */
    double  srp_occulter_radius;     /* km                           */
    double  sun_radius;              /* km                           */

    /* atmospheric drag - placeholder, not yet implemented            */
    int     enable_drag;

    /* Time anchor: Ephemeris Time (seconds past J2000 TDB) at integrator
     * t = 0. The ephemeris query argument is simply
     *   et = et0 + t
     * (no unit conversion: integrator time and ET share the same units).
     * Use ET_FROM_JD(jd) from spody_const.h to convert from a JD epoch.
     * et0 = 0 corresponds to the J2000 epoch itself. */
    double  et0;
} ForceModelContext;

/* ============================================================
 * Atomic force functions
 * ============================================================
 * Each writes its contribution into acc_out. They are stateless,
 * reentrant, and may be called individually if you want to build
 * a custom RHS instead of using the default one.
 */

/* Two-body gravitational acceleration around a point mass at the origin.
 *   a = -mu * r / |r|^3                                                  */
void spody_force_twobody(double mu, const double r_sat[3], double acc_out[3]);

/* Disturbing part of the central-body spherical-harmonic gravity.
 * Calls `get_R` to obtain the ICRF<->body-fixed rotations at `et`,
 * transforms `r_sat` into the body-fixed frame, evaluates the
 * harmonic series with `hg`, and transforms the resulting
 * acceleration back to ICRF. `get_R` MUST be non-NULL; for the
 * Moon use `spody_bf_rotation_moon`.
 *
 * Note: returns ONLY the J2-and-higher disturbing acceleration; the
 * 2-body term must be summed separately (this matches the contract
 * of spody_get_hgaccbodyfixed_hpc, which starts the recurrence at n=2). */
void spody_force_sphericalharmonics(HarmonicGravity *hg, MappedEphemeris *eph,
                                    double et, const double r_sat[3],
                                    spody_bf_rotation_fn get_R,
                                    double acc_out[3]);

/* Cowell formulation of the third-body perturbation (Battin):
 *   a = -mu_3 * ( (r_sat - r_3) / |r_sat - r_3|^3  +  r_3 / |r_3|^3 )
 * where r_3 is the third body's position relative to the central body. */
void spody_force_thirdbody_cowell(double mu_3, const double r_3[3],
                                  const double r_sat[3], double acc_out[3]);

/* Cannonball SRP model. The caller supplies the eclipse fraction (1 in
 * full sunlight, 0 in full umbra) and the satellite-to-sun displacement.
 * Acceleration points away from the Sun (i.e. from sun toward sat):
 *   a = - SOLAR_LUMINOSITY_4PIC * Cr * (A/m) * fraction / |r|^3 * r_sat_to_sun */
void spody_force_srp(const Spacecraft *sat, double fraction_sunlight,
                     const double r_sat_to_sun[3], double acc_out[3]);

/* Atmospheric drag - PLACEHOLDER.
 * Currently writes acc_out = {0,0,0}. The real implementation will need
 * an atmospheric density model and the central body's rotation rate to
 * compute the air-relative velocity. The signature may change when the
 * model is wired in. */
void spody_force_drag(const Spacecraft *sat,
                      const double r_sat[3], const double v_sat[3],
                      double acc_out[3]);

/* ============================================================
 * Composite RHS (plugs into the integrator)
 * ============================================================
 * Fits the spody_rhs_fn signature: pass its address to
 * spody_setup_integrator together with a ForceModelContext as the
 * opaque user pointer. State layout: y = [r(3), v(3)], dim = 6.
 *
 * Sums perturbations in ascending order of typical magnitude on
 * a low-orbit (SRP -> drag -> 3rd body -> harmonics -> 2-body) to
 * minimise round-off accumulation in long propagations.
 */
int spody_force_rhs_default(double t, const double *y, double *dy, void *user);

/* ============================================================
 * Force breakdown (post-step diagnostic)
 * ============================================================
 * Decomposed accelerations for a given (t, y), produced by
 * spody_force_breakdown(). Intended to be filled AFTER an accepted
 * step (typically right after spody_propagate_onestep), so that the
 * diagnostic is bit-coherent with the integrator state.
 *
 * Cost: roughly one extra RHS evaluation per call. Keep it OUTSIDE
 * the integrator hot loop (i.e. do not call from within rhs_default
 * or from any per-stage callback). Calling it once per accepted step
 * adds ~12% to the wall-time of the propagation -- negligible if
 * you only log every k steps.
 *
 * Bit-equivalence: acc_total reproduces the result of rhs_default at
 * the same (t, y), with the same summation order
 * (SRP + drag + third_total + harmonics + 2body).
 *
 * The whole struct is written as one record into the breakdown binary
 * log -- including n_third and the per-body array. Internal padding
 * (typically 4 bytes between n_third and acc_thirdbody) is also
 * written verbatim and ignored by readers; this avoids any need for
 * field-by-field copies between two parallel layouts.
 */
typedef struct {
    double t;                                          /* sim time (s from t=0)   */

    double acc_total[3];                               /* sum of all forces       */
    double acc_2body[3];                               /* central two-body        */
    double acc_sphericalharmonics[3];                  /* J2+ disturbing          */
    double acc_thirdbody_total[3];                     /* sum of all third bodies */
    int    n_third;                                    /* # per-body entries used */
    double acc_thirdbody[SPODY_FM_MAX_THIRD][3];       /* per-body breakdown      */
    double acc_srp[3];                                 /* SRP                     */
    double acc_drag[3];                                /* drag (placeholder)      */
    double eclipse_fraction;                           /* 1=full sun, 0=full umbra */
} ForceBreakdown;

/* Re-evaluate the force decomposition on the given (t, y) and write the
 * per-component accelerations + eclipse fraction into `bd`. Same gating
 * rules as spody_force_rhs_default (NULL pointers / zero counters disable
 * a force). Intended for use AFTER spody_propagate_onestep, never inside
 * the RHS hot path. */
void spody_force_breakdown(const ForceModelContext *ctx,
                           double t, const double *y,
                           ForceBreakdown *bd);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_FORCEMODELS_H */
