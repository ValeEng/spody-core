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
#ifndef SPODY_MISSION_H
#define SPODY_MISSION_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spody_forcemodels.h"   /* Spacecraft, ForceModelContext, ForceBreakdown */
#include "spody_integrators.h"   /* IntegratorAllData                            */
#include "spody_events.h"        /* SpodyEvent, EventRecord                      */
#include "spody_io.h"            /* SpodyLogBuffer                                */

/* ============================================================
 * Mission orchestration
 *
 * The Mission ties together the four moving parts of a propagation:
 *   - the Spacecraft descriptor (read-only physical parameters)
 *   - the ForceModelContext (which forces are active, and how)
 *   - the IntegratorAllData (numerical method + state)
 *   - two optional binary output streams (trajectory + breakdown)
 *
 * Lifecycle:
 *
 *   Spacecraft sat = ...; spody_init_Spacecraft(&sat);
 *   ForceModelContext ctx = ...;
 *   IntegratorAllData integ;
 *   spody_setup_integrator(&integ, SPODY_INTEG_RK45, NULL, 6,
 *                          spody_force_rhs_default, &ctx);
 *   spody_set_integrator_state(&integ, t_start, y0);
 *
 *   Mission m;
 *   spody_setup_Mission(&m, &sat, &ctx, &integ,
 *                       t_start, t_end,
 *                       "traj.spdytraj",   // or NULL = no trajectory log
 *                       "brkd.spdybrkd");  // or NULL = no breakdown log
 *
 *   spody_run_propagate(&m);     // loops onestep until t_end, saves each step
 *
 *   spody_free_Mission(&m);      // flush + close output files
 *
 * Threading: not safe to share a Mission across threads. Each worker
 * thread must own its own Mission with its own ctx (containing its
 * own MappedEphemeris + HarmonicGravity handles).
 * ============================================================ */

/* Mission uses the standard log magics and the default buffer capacity
 * declared in spody_io.h:
 *   SPODY_MAGIC_TRAJECTORY ("SPDYTRAJ"), SPODY_MAGIC_BREAKDOWN ("SPDYBRKD"),
 *   SPODY_LOG_BUFFER_CAPACITY (1024 records). */

/* ============================================================
 * On-disk record types
 * ============================================================ */

/* Trajectory record: 7 doubles = 56 bytes. */
typedef struct {
    double t;            /* simulation time, seconds from ctx->t0_jd */
    double y[6];         /* state: r(3), v(3) -- km, km/s            */
} TrajectoryRecord;

/* Breakdown record: spody_run_propagate writes one ForceBreakdown
 * struct per record into the breakdown binary log -- the full struct,
 * including n_third and the per-body array (and any padding). The
 * record_size of the breakdown buffer is therefore sizeof(ForceBreakdown).
 * No separate struct or field-by-field copy is needed. */

/* ============================================================
 * Mission struct
 * ============================================================ */
typedef struct {
    /* configuration (caller-owned, kept by reference) */
    const Spacecraft  *sat;
    ForceModelContext *ctx;
    IntegratorAllData *integ;

    double t_start;
    double t_end;

    /* events: caller-owned array, checked after each accepted step */
    SpodyEvent *events;                /* may be NULL if n_events == 0 */
    int         n_events;

    /* output streams */
    SpodyLogBuffer trajectory_log;     /* zeroed if no trajectory file */
    SpodyLogBuffer breakdown_log;      /* zeroed if no breakdown file  */
    SpodyLogBuffer events_log;         /* zeroed if no events file     */
    int            log_trajectory;     /* 1 if trajectory_log is open  */
    int            log_breakdown;      /* 1 if breakdown_log is open   */
    int            log_events;         /* 1 if events_log is open      */

    /* runtime bookkeeping */
    long n_steps_done;
    long n_records_saved;
    long n_events_triggered;
} Mission;

/* ============================================================
 * Lifecycle
 * ============================================================ */

/* Initialize a Mission. Each path may be NULL to disable that log.
 *
 *   sat          : caller-owned Spacecraft (must outlive the Mission).
 *   ctx          : caller-owned ForceModelContext (per-thread).
 *   integ        : caller-owned IntegratorAllData (per-thread, already set up
 *                  with rhs = spody_force_rhs_default and user = ctx).
 *   t_start, t_end : propagation horizon, seconds from ctx->et0.
 *   trajectory_path : output binary for TrajectoryRecord. NULL = disabled.
 *   breakdown_path  : output binary for ForceBreakdown.   NULL = disabled.
 *   events_path     : output binary for EventRecord.      NULL = disabled.
 *
 * The integrator is NOT (re)set internally; the caller is expected to call
 * spody_set_integrator_state(integ, t_start, y0) before spody_run_propagate
 * if the state needs initialisation. Events are configured separately by
 * setting m->events / m->n_events directly after this call (the array is
 * caller-owned). */
int spody_setup_Mission(Mission *m,
                        const Spacecraft *sat,
                        ForceModelContext *ctx,
                        IntegratorAllData *integ,
                        double t_start, double t_end,
                        const char *trajectory_path,
                        const char *breakdown_path,
                        const char *events_path);

/* Return codes from spody_run_propagate. */
#define SPODY_RUN_OK                 0   /* completed up to t_end           */
#define SPODY_RUN_STOPPED_BY_EVENT   1   /* an event with action=STOP fired */
/* (any negative value is forwarded from the integrator/log layer)         */

/* Drive the integrator from integ->t up to t_end, calling
 * spody_propagate_onestep() in a loop. After each accepted step:
 *   - if trajectory log is enabled, append a TrajectoryRecord;
 *   - if breakdown log is enabled, evaluate spody_force_breakdown
 *     and append a ForceBreakdown record;
 *   - for every configured event, evaluate the predicate. On a trigger:
 *       * if action includes LOG and the events file is enabled, append
 *         an EventRecord;
 *       * if action includes STOP, return SPODY_RUN_STOPPED_BY_EVENT
 *         after writing pending logs.
 *
 * Returns SPODY_RUN_OK on full completion. */
int spody_run_propagate(Mission *m);

/* Flush + close any open log files. Always safe to call (no-op on
 * a Mission that didn't enable any log). */
int spody_free_Mission(Mission *m);

/* ============================================================
 * Convenience: dump the binary outputs to CSV
 *
 * Wrappers over spody_log_dump_csv that pre-fill the column names
 * for the two known record types. Useful for quick analysis.
 * ============================================================ */
int spody_mission_trajectory_to_csv(const char *bin_path, const char *csv_path,
                                    int decimate);
int spody_mission_breakdown_to_csv (const char *bin_path, const char *csv_path,
                                    int decimate);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_MISSION_H */
