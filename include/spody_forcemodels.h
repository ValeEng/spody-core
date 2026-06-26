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

#include "spody_atmosphere.h"
#include "spody_const.h"
#include "spody_eclipse.h"
#include "spody_eop.h"
#include "spody_ephemeris.h"
#include "spody_harmonics.h"
#include "spody_integrators.h"   /* spody_rhs_fn */

/* Forward declaration for Earth orientation per-thread handle.
 * The full type lives in spody_earth_orientation.h (P2.2c); we
 * only need the pointer here to keep the rotation callback
 * generic without dragging the IAU 2006 tables into every TU. */
struct MappedIAU2006;
typedef struct MappedIAU2006 MappedIAU2006;

/* Forward declaration for ForceModelContext: the rotation callback
 * typedef below takes a pointer to it, and the struct definition
 * itself stores a function pointer of that type, so we need both
 * names available before either is complete. */
struct ForceModelContext;
typedef struct ForceModelContext ForceModelContext;

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
 * The callback receives the full ForceModelContext so each body
 * can pick the runtime inputs it needs:
 *   - Moon  -> reads ctx->eph (DE440 libration in slot 12)
 *   - Earth -> reads ctx->eop + ctx->iau2006
 *               (IERS EOP + IAU 2006/2000A series)
 *   - Mars  -> [TODO] IAU 2009 based, no runtime inputs
 *
 * Bodies whose orientation comes from purely analytic models with
 * no per-step inputs simply ignore ctx. The contract is "pick what
 * you need, leave the rest", which keeps the kernel kernel
 * body-agnostic without forcing every body to declare a wrapper
 * struct of its required inputs. */
typedef void (*spody_bf_rotation_fn)(const ForceModelContext *ctx, double et,
                                      double R_icrf_to_bf[3][3],
                                      double R_bf_to_icrf[3][3]);

/* Concrete provider: Moon Principal Axes from DE440 libration
 * angles. Equivalent to chaining
 *   spody_get_lunarlibrationangles(ctx->eph, et, angles);
 *   spody_getrotmatrix_icrf2moonpa(angles..., R_icrf_to_bf);
 *   spody_getrotmatrix_moonpa2icrf(angles..., R_bf_to_icrf);
 * but exposed as a function pointer so the force-model layer can
 * call it generically. Reads `ctx->eph` (MUST be non-NULL); the
 * other ctx fields are ignored. */
void spody_bf_rotation_moon(const ForceModelContext *ctx, double et,
                             double R_icrf_to_bf[3][3],
                             double R_bf_to_icrf[3][3]);

struct ForceModelContext {
    /* central body (the body the satellite orbits) */
    double  mu_central;          /* km^3/s^2                          */
    double  R_central;           /* km - mean radius (impact check)   */
    int     naif_central;        /* NAIF id, e.g. 301 (Moon), 399 (Earth) */

    /* Body-fixed orientation provider for `naif_central`. Used by
     * the spherical-harmonics force to rotate state ICRF <-> body
     * fixed frame at every RHS evaluation. Application fills this
     * based on the parsed central-body name; for Moon use
     * `spody_bf_rotation_moon`, for Earth `spody_bf_rotation_earth`.
     * MUST be non-NULL when `hg` is set. */
    spody_bf_rotation_fn get_bf_rotation;

    /* spacecraft */
    const Spacecraft *sat;

    /* spherical harmonics on the central body (NULL = disabled) */
    HarmonicGravity *hg;

    /* ephemeris-driven perturbations (NULL = disabled). Must be
     * non-NULL whenever hg, n_third > 0, or enable_srp are active.
     * Also consumed by `spody_bf_rotation_moon` to read the lunar
     * libration angles from DE440 slot 12. */
    MappedEphemeris *eph;

    /* Per-thread Earth-orientation handles, consumed by
     * `spody_bf_rotation_earth` (P2.2). NULL when the central body
     * is not Earth: the Moon callback ignores them, so leaving
     * them unset has no effect on lunar runs. The application sets
     * them up from the TOML's `force_model.eop_file` and
     * `force_model.iau2006_dir` fields, which are written only
     * when central_body = "Earth". */
    MappedEOP     *eop;
    MappedIAU2006 *iau2006;

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

    /* atmospheric drag.
     *
     * `enable_drag` is the runtime toggle (parsed from the TOML).
     * `atmosphere` holds the density callback registered by the
     * central body (NRLMSISE-00 for Earth, MCD for Mars when it
     * ships, ...); a NULL pointer means the body has no atmosphere
     * and `spody_force_drag` returns zero unconditionally. The
     * space-weather handle is the per-thread query cache over the
     * shared MappedSpaceWeatherData parsed by the app once at
     * startup; NULL when drag is off or when the chosen model does
     * not consume space weather.
     *
     * `body_spin_rad_s` is the central body's sidereal rotation rate
     * about its body-fixed +Z axis (Earth: EARTH_ROT_RATE_RADPS).
     * Used to derive omega x r for the air-relative velocity in the
     * drag formula. Bodies without a rotation rate produce zero
     * drag (the velocity relative to a stationary atmosphere is the
     * inertial velocity itself, which is fine, but in practice the
     * drag toggle should be off in that case). */
    int                 enable_drag;
    SpodyAtmosphere    *atmosphere;
    MappedSpaceWeather *space_weather;
    double              body_spin_rad_s;

    /* Time anchor: Ephemeris Time (seconds past J2000 TDB) at integrator
     * t = 0. The ephemeris query argument is simply
     *   et = et0 + t
     * (no unit conversion: integrator time and ET share the same units).
     * Use ET_FROM_JD(jd) from spody_const.h to convert from a JD epoch.
     * et0 = 0 corresponds to the J2000 epoch itself. */
    double  et0;

    /* ---- CR3BP fields ----
     * Used only when the RHS is `spody_force_rhs_cr3bp`. Zero (or
     * uninitialised) in high-fidelity runs -- the HF RHS never reads
     * them. State y in CR3BP runs is interpreted in the synodic
     * rotating frame anchored on the barycenter of the two primaries;
     * positions in km, velocities in km/s. The system is autonomous
     * so `et0` is not consulted by the CR3BP RHS.
     *
     *   cr3bp_mu1, cr3bp_mu2 : gravitational parameters of the two
     *                          primaries (km^3/s^2). Convention is
     *                          cr3bp_mu1 >= cr3bp_mu2 (primary 1 is
     *                          the bigger one), but the RHS is
     *                          symmetric so the order is not enforced.
     *   cr3bp_L              : primary-primary separation (km),
     *                          assumed fixed.
     *   cr3bp_omega, _x1, _x2: derived caches populated by
     *                          spody_init_CR3BPContext from the three
     *                          inputs above. Do not set them by hand.
     *                          omega = sqrt((mu1+mu2)/L^3),
     *                          x1 = -(mu2/(mu1+mu2)) * L,
     *                          x2 = +(mu1/(mu1+mu2)) * L.
     */
    double cr3bp_mu1;
    double cr3bp_mu2;
    double cr3bp_L;
    double cr3bp_omega;
    double cr3bp_x1;
    double cr3bp_x2;
};

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
 * Pulls the harmonic-gravity handle from `ctx->hg`, evaluates the
 * ICRF<->body-fixed rotations at `et` via `ctx->get_bf_rotation`
 * (passing `ctx` through so the callback can pick its body-specific
 * inputs from ctx->eph / ctx->eop / ctx->iau2006 as needed),
 * transforms `r_sat` into the body-fixed frame, evaluates the
 * harmonic series, and transforms the resulting acceleration back
 * to ICRF.
 *
 * `ctx->hg` and `ctx->get_bf_rotation` MUST be non-NULL.
 *
 * Note: returns ONLY the J2-and-higher disturbing acceleration; the
 * 2-body term must be summed separately (this matches the contract
 * of spody_get_hgaccbodyfixed_hpc, which starts the recurrence at n=2). */
void spody_force_sphericalharmonics(const ForceModelContext *ctx,
                                    double et, const double r_sat[3],
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

/* Atmospheric drag in ICRF (km/s^2).
 *
 *   a_drag = -0.5 * rho * |v_rel|^2 * Cd * (A/m) * v_rel_hat
 *   v_rel  = v_sat - omega_central_icrf x r_sat
 *
 * Queries `ctx->atmosphere->density(...)` for rho at the satellite
 * position (mapped from ICRF to body-fixed via `ctx->get_bf_rotation`
 * so the density model sees coordinates in its native frame), and
 * builds omega_central_icrf as `body_spin_rad_s * (R_bf_to_icrf @
 * +z_bf)`.
 *
 * Returns acc_out = {0,0,0} when ANY of (ctx->sat, ctx->atmosphere,
 * ctx->get_bf_rotation, ctx->body_spin_rad_s > 0) is missing, or
 * when the density callback fails -- so the call site can wire it
 * unconditionally under the drag toggle. */
void spody_force_drag(const ForceModelContext *ctx, double et,
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

/* CR3BP RHS in the synodic rotating frame, dimensional units.
 *
 * Reads ONLY the cr3bp_* fields of the context (all HF fields may be
 * NULL/zero). State layout: y = [r(3), v(3)] in km, km/s in the
 * rotating frame whose +x axis points from the bigger primary toward
 * the smaller one, +z aligned with the orbital angular momentum.
 * Time `t` is irrelevant (autonomous system); `et0` is not consulted.
 *
 * Uses the pre-cached derived quantities (cr3bp_omega, cr3bp_x1,
 * cr3bp_x2) populated by spody_init_CR3BPContext. Caller must have
 * called that init once after filling cr3bp_mu1 / cr3bp_mu2 / cr3bp_L. */
int spody_force_rhs_cr3bp(double t, const double *y, double *dy, void *user);

/* Pre-compute the derived caches (omega, x1, x2) from cr3bp_mu1,
 * cr3bp_mu2, cr3bp_L. Call once after filling the three input fields.
 * Sets the derived caches to zero if any input is non-positive, so an
 * uninitialised CR3BP slot in an HF run stays harmlessly zero. */
void spody_init_CR3BPContext(ForceModelContext *ctx);

/* Transform a state expressed in one primary's local inertial frame
 * (origin at the primary, axes non-rotating, primary treated as
 * stationary) into the CR3BP synodic rotating frame at t = 0.
 *
 * At t = 0 the synodic axes coincide with the underlying inertial
 * axes by convention, so the rotation reduces to identity and only
 * a translation by the primary's synodic position plus the omega-x-r
 * velocity correction need to be applied. This is the conversion the
 * TOML input layer uses to seed the CR3BP integrator from a Keplerian
 * initial condition referred to one of the primaries.
 *
 *   r_primary_inertial, v_primary_inertial : input state expressed in
 *       the chosen primary's inertial frame (km, km/s)
 *   mu1_km3_s2, mu2_km3_s2 : GM of primary_1 (bigger) and primary_2
 *   L_km                   : primary-primary separation
 *   primary_index          : 1 = primary_1, 2 = primary_2
 *   r_synodic, v_synodic   : output state in synodic frame (km, km/s) */
void spody_inertial_to_cr3bp_synodic(
        const double r_primary_inertial[3],
        const double v_primary_inertial[3],
        double mu1_km3_s2, double mu2_km3_s2, double L_km,
        int    primary_index,
        double r_synodic[3], double v_synodic[3]);

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
