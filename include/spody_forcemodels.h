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

/* The SRP occulter list is "central body + third bodies", so the
 * eclipse machinery must be able to take one more body than the
 * third-body cap. */
#if SPODY_FM_MAX_THIRD + 1 > SPODY_ECL_MAX_OCCULTERS
#error "SPODY_ECL_MAX_OCCULTERS too small for SPODY_FM_MAX_THIRD + 1"
#endif

/* ============================================================
 * Empirical acceleration in RIC
 * ============================================================
 * An acceleration given as a table of nodes in the satellite's own
 * radial / in-track / cross-track axes (r_hat = r/|r|, c_hat =
 * r x v/|r x v|, i_hat = c_hat x r_hat, rotation only:
 * spody_getrotmatrix_ric2icrf), linearly interpolated in ET between
 * the nodes and held at the end values outside them. The classical
 * "empirical acceleration" of orbit determination (Montenbruck & Gill,
 * "Satellite Orbits", Springer, 2000); the Monte Carlo fills it with
 * Gauss-Markov process noise. Not owned by the context. */
typedef struct {
    const double *et;      /* node epochs, ET s, strictly ascending   */
    const double *a_ric;   /* 3 per node (R, I, C), km/s^2            */
    size_t        n;       /* >= 1                                    */
} SpodyEmpiricalAccel;

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

/* Angular velocity of the central body's body-fixed frame, in ICRF
 * (rad/s), at `et` -- the omega of the transport theorem
 *   v_icrf = R_bf_to_icrf v_rot + omega x r_icrf
 * that turns a velocity measured in the rotating frame (ECEF-style)
 * into an inertial one.
 *
 * Every body, the Earth included: the rotation R(t+h) R(t-h)^T of
 * ctx->get_bf_rotation over +-SPODY_BF_OMEGA_FD_STEP_S, read off in
 * axis-angle form (exact for a fixed axis). For the Earth that is the
 * true rotation of the engine's ITRS chain -- about the CIP, tilted
 * from the ITRS z axis by polar motion, at the EOP rate -- not a
 * nominal spin about z; for the Moon it follows the DE440 libration.
 * The GNSS converters call this same function. Requires
 * ctx->get_bf_rotation. */
void spody_bf_angular_velocity_icrf(const ForceModelContext *ctx, double et,
                                    double omega_icrf[3]);

/* Solid-body tide of the central body: the potential its elastic
 * deformation adds under the pull of the tide-raising bodies (IERS
 * Conventions 2010 sec. 6.2.1, frequency-independent step, eq.
 * 6.6-6.7). Every evaluation turns the raisers' body-fixed positions
 * into corrections of the normalized coefficients,
 *   dC_nm - i dS_nm = k_nm/(2n+1) sum_j (GM_j/GM)(R/r_j)^(n+1)
 *                     Pbar_nm(sin phi_j) exp(-i m lambda_j),
 * for n = 2..max_degree, plus the degree-4 terms a degree-2 forcing
 * induces through k^(+)_2m (all zero: no coupling), and returns the
 * acceleration of that small field. The direct pull of the raisers
 * on the satellite is not here: it is the third-body force.
 *
 * Filled by the application from the central body's registry row;
 * the engine holds no per-body knowledge. `gm` and `r_ref` are the
 * harmonics file's, the ones its coefficients are normalized with.
 * `dc20_perm` is the permanent part a zero-tide field already holds
 * and the correction must leave out; 0 for a tide-free field. */
#define SPODY_TIDE_MAX_RAISERS 2
typedef struct {
    int    n_raisers;
    int    raiser_naif[SPODY_TIDE_MAX_RAISERS];
    double raiser_mu[SPODY_TIDE_MAX_RAISERS];   /* km^3/s^2           */
    int    max_degree;                          /* 2 or 3             */
    double k_re[4][4];                          /* [n][m], n = 2..3   */
    double k_im[4][4];
    double kplus[3];                            /* k^(+)_2m, m = 0..2 */
    double gm;                                  /* km^3/s^2           */
    double r_ref;                               /* km                 */
    double dc20_perm;
} SpodySolidTides;

/* Earth radiation pressure: albedo (reflected sunlight, dayside only)
 * and infrared (thermal emission) of the Earth on a cannonball
 * satellite, with the latitude- and season-dependent albedo and
 * emissivity of Knocke et al. (1988). Each surface element is a
 * Lambertian emitter of exitance
 *   M = a F cos(zeta_sun) [if lit] + e F / 4,  F = solar flux at Earth,
 * seen by the satellite with irradiance (M/pi) cos(theta) dA / d^2,
 * theta = angle between the element's normal and the direction to the
 * satellite. Integrated over the visible disk in solid angle.
 * Acceleration = Cr (A/m) E / c, along each element -> satellite ray.
 * Requires ctx->eph (Sun), ctx->get_bf_rotation (latitude) and a
 * spacecraft with am_srp / Cr. Earth only (the application refuses
 * any other central body). */
void spody_force_earthradiation(const ForceModelContext *ctx, double et,
                                const double r[3], double acc[3]);

/* General relativity, Schwarzschild term of the central body (IERS
 * Conventions 2010 eq. 10.12, first line):
 *   a = GM/(c^2 r^3) [ (2(beta+gamma) GM/r - gamma v.v) r
 *                      + 2(1+gamma) (r.v) v ]
 * r, v relative to the central body (ICRF, km, km/s), GM = mu_central.
 * The Lense-Thirring and de Sitter lines of the same equation are not
 * modelled (1e-11..1e-12 of gravity, IERS). */
void spody_force_relativity(double mu, const double r[3], const double v[3],
                            double acc[3]);

/* Acceleration (ICRF, km/s^2) of the solid tide at `et` for a
 * satellite at `r` (ICRF, central-body centred). Requires ctx->tides,
 * ctx->eph and ctx->get_bf_rotation. */
void spody_force_solidtides(const ForceModelContext *ctx, double et,
                            const double r[3], double acc[3]);

struct ForceModelContext {
    /* central body (the body the satellite orbits) */
    double  mu_central;          /* km^3/s^2                          */
    double  R_central;           /* km - equatorial radius            */
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

    /* solid-body tide of the central body (NULL = disabled) */
    const SpodySolidTides *tides;

    /* general relativity, Schwarzschild term of the central body */
    int     enable_relativity;

    /* Earth radiation pressure (albedo + infrared, Knocke) */
    int     enable_earthradiation;

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

    /* Solar radiation pressure (cannonball) with a multi-occulter
     * eclipse. The occulter list is application policy -- typically
     * the central body plus every third body, never the Sun itself --
     * and is built once, before the run: it is epoch-independent, so
     * the RHS never rebuilds it. The force combines the whole list
     * into a single lit fraction (see spody_get_satlitfraction), so
     * the Earth shadowing a lunar orbiter and the Moon transiting the
     * Sun for an Earth orbiter both land in the acceleration.
     *
     * srp_n_occulters == 0 disables shadow modelling: the satellite
     * is then lit at all times. Entries with r_eq <= 0 are ignored
     * individually, which is how a run can keep a body in the list
     * without letting it cast a shadow. Each shape is a sphere or a
     * spheroid (force_model.body_shape), its pole in the integration
     * frame, fixed for the run. */
    int     enable_srp;
    int     srp_n_occulters;
    int     srp_occulter_naif [SPODY_ECL_MAX_OCCULTERS];
    SpodyBodyShape srp_occulter_shape[SPODY_ECL_MAX_OCCULTERS];
    double  sun_radius;                                    /* km */

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

    /* Optional density calibration k(t), multiplied onto the model
     * density inside `spody_force_drag` (see MappedDensityScale in
     * spody_atmosphere.h). NULL means k = 1 (uncalibrated model),
     * the default for a zero-initialised context. */
    const MappedDensityScale *density_scale;

    /* Optional empirical acceleration in RIC (SpodyEmpiricalAccel),
     * added after relativity in spody_force_rhs_default and in the
     * total of spody_force_breakdown. NULL (the zero-initialised
     * default) = none, and the sums are bit-identical to before. */
    const SpodyEmpiricalAccel *empirical_accel;

    /* Shadow contact found by spody_next_force_discontinuity inside a
     * step that has to be redone (integrator t, per thread), kept until
     * the integration passes it. A value <= t means none; the stepping
     * loop sets it to INFINITY before each run. */
    double  disc_contact;

    /* Time anchor: Ephemeris Time (seconds past J2000 TDB) at integrator
     * t = 0. The ephemeris query argument is simply
     *   et = et0 + t
     * (no unit conversion: integrator time and ET share the same units).
     * Use ET_FROM_JD(jd) from spody_const.h to convert from a JD epoch.
     * et0 = 0 corresponds to the J2000 epoch itself. */
    double  et0;

    /* Integration time scale. 0 (default): the integrator's t is TDB
     * seconds past et0, et = et0 + t, as above. 1: t is TT seconds past
     * the same instant, tt = tt0 + t with tt0 = et0 - (TDB-TT)(et0), and
     * et comes through spody_ctx_et -- the time coordinate IERS 2010
     * (TN36 sec. 10.3) prescribes for geocentric satellite equations of
     * motion. Every consumer of absolute time goes through spody_ctx_et;
     * every time written to a file goes through spody_ctx_label, so the
     * files keep "ET - et0" in both modes. */
    int     time_scale_tt;
    double  tt0;

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

/* Empirical acceleration in ICRF (km/s^2): the table `ea` at `et`
 * (linear between nodes, end values outside), rotated from the RIC
 * axes of (r, v). Zero when ea is NULL or empty, or when the RIC axes
 * are undefined (r = 0 or r parallel to v). */
void spody_force_empirical(const SpodyEmpiricalAccel *ea, double et,
                           const double r[3], const double v[3],
                           double acc_out[3]);

/* Atmospheric drag in ICRF (km/s^2).
 *
 *   a_drag = -0.5 * rho * |v_rel|^2 * Cd * (A/m) * v_rel_hat
 *   v_rel  = v_sat - omega_central_icrf x r_sat
 *
 * Queries `ctx->atmosphere->density(...)` for rho at the satellite
 * position (mapped from ICRF to body-fixed via `ctx->get_bf_rotation`
 * so the density model sees coordinates in its native frame), scales
 * it by the optional calibration table `ctx->density_scale` (k = 1
 * when NULL), and builds omega_central_icrf as `body_spin_rad_s *
 * (R_bf_to_icrf @ +z_bf)`.
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

/* Time of the integrator (see ctx->time_scale_tt) as ET, the argument
 * of every ephemeris, rotation and space-weather query. TDB mode:
 * et0 + t, unchanged. TT mode: tt = tt0 + t, et = tt + (TDB-TT)(tt):
 * the periodic term is evaluated at tt instead of et, 1.7 ms away, and
 * moves by at most 3.3e-10 s per second, so et is off by < 6e-13 s
 * (5 nm on a LEO) -- one deltet per RHS call instead of two. */
double spody_ctx_et(const ForceModelContext *ctx, double t);

/* The integrator's t as written to files: ET - et0, so a file means
 * the same thing whatever the integration time scale. TDB mode returns
 * t itself (bit for bit). */
double spody_ctx_label(const ForceModelContext *ctx, double t);

/* Inverse of spody_ctx_label: the integrator's t at file time `label`
 * (ET - et0). TDB mode returns `label` itself. Used to put the output
 * grid, and the end of the run, on round ET labels. */
double spody_ctx_t_of_label(const ForceModelContext *ctx, double label);

/* Retunes the harmonics truncation degree for the coming step.
 *
 * Call it from the stepping loop, immediately BEFORE each
 * spody_propagate_onestep, passing the integrator's current t / y / h
 * and the same ForceModelContext the RHS uses. Never call it from
 * inside the RHS: the stages of one step would then sample different
 * vector fields, the embedded error estimate would read the model
 * jump as truncation error, and the controller would shrink h
 * fighting a discontinuity that is not in the dynamics. Once per
 * step, before any stage, keeps the whole step on ONE smooth field --
 * the property the Runge-Kutta order derivation assumes.
 *
 * A loop that never calls it keeps the fixed degree and reproduces
 * earlier results bit for bit; there is no state to initialise.
 *
 * Writes `ctx->hg->N_eval` and returns 1 when the degree differs from
 * the one the previous call left, 0 otherwise -- including when ctx->hg
 * is NULL (harmonics off), so it is safe to call unconditionally. A
 * changed degree is a changed vector field, and an integrator that
 * carries a derivative across steps (RK45's FSAL stage) must be told:
 * a loop that calls this feeds the return value to
 * spody_integrator_invalidate_fsal. That costs one RHS evaluation on
 * the steps where the degree moves and keeps every stage of every step
 * on the field the step was chosen for.
 *
 * The degree-n term of the potential decays as (R_ref/r)^n, so
 * requiring it below a relative threshold eps gives the degree
 * directly:
 *
 *     N(r) = ln(1/eps) / ln(r / R_ref)
 *
 * Only the ratio r/R_ref appears, so there is no per-body or
 * per-model constant to calibrate: the same expression serves any
 * central body and any coefficient set. It deliberately DROPS the
 * decay of the coefficients themselves (Kaula-type sigma_n ~ K/n^2).
 * Including that would lower the predicted degree and make the result
 * depend on how well a given field follows the power law; keeping
 * only the geometric factor makes the answer an upper bound instead,
 * which is the safe direction for a truncation rule.
 *
 * Being an upper bound is an argument, not a proof: a field whose
 * coefficients do not decay smoothly with degree (resonances, an
 * unusually rough body) can in principle need more. Check a new
 * coefficient set against a per-degree convergence walk on its real
 * Cnm/Snm before relying on this.
 *
 * The degree is picked from a LOWER BOUND on the radius the step can
 * reach, not from the radius it starts at:
 *
 *     r_bound = |r| - SPODY_HG_ADAPTIVE_STEP_MARGIN * |v| * |h|
 *
 * Without it, a step starting where the field is weak would enter its
 * stages with a degree too low for the part of the step that descends
 * toward the body. The bound is what makes the choice safe rather
 * than merely typical, and it is cheap: where steps are long the
 * velocity is low, so |v|*h stays a small fraction of the radius.
 *
 * The result is capped at the loaded degree (ctx->hg->hgd->N), so
 * this can only ever ask for LESS work than the fixed setting, never
 * for coefficients that were not loaded. */
int spody_adapt_hgdegree(double t, const double *y, double h, void *user);

/* First known discontinuity of the force model after the start of the
 * last accepted step (integrator time), or INFINITY. Called once per
 * step, right after it:
 *
 *   - NRLMSISE-00 inputs (drag on, space weather loaded): the 3-hour
 *     Ap bins, the daily F10.7 and the day of year change on the
 *     3-hour UTC grid; the next boundary after integ->t is returned.
 *   - Shadow contacts (SRP with occulters): the first penumbra or
 *     umbra contact of each occulter inside the last step
 *     [t_old, t] (spody_get_eclipse_residual, Brent on
 *     spody_dense_state_rv6), remembered in ctx->disc_contact until
 *     the integration passes it. Contacts within 2 SPODY_DISC_STOP_EPS_S
 *     of t_old belong to the crossing step and are not reported.
 *
 * A result <= integ->t means the last step crossed a discontinuity:
 * the caller redoes it from (t_old, y_old) and stops there. A
 * Runge-Kutta step across a jump in f loses its order (Hairer,
 * Norsett, Wanner, "Solving Ordinary Differential Equations I", 2nd
 * ed., Springer, 1993, Sect. II.6). Before the first step (h_old == 0
 * or t == t_old) only the known discontinuities are returned. */
double spody_next_force_discontinuity(ForceModelContext *ctx,
                                      const IntegratorAllData *integ);

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
 * (SRP, Earth radiation, drag, each third body in turn, harmonics,
 * solid tide, relativity, then 2body).
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
    double acc_solidtides[3];                          /* solid-body tide (appended:
                                                        * SPDYACC_ v2) */
    double acc_relativity[3];                          /* general relativity,
                                                        * Schwarzschild (appended:
                                                        * SPDYACC_ v3) */
    double acc_earthradiation[3];                      /* Earth albedo + infrared
                                                        * (appended: SPDYACC_ v4) */
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
