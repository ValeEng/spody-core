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
#include <string.h>
#include "spody_mission.h"

/* ============================================================
 * Setup / teardown
 * ============================================================ */
int spody_setup_Mission(Mission *m,
                        const Spacecraft *sat,
                        ForceModelContext *ctx,
                        IntegratorAllData *integ,
                        double t_start, double t_end,
                        const char *trajectory_path,
                        const char *breakdown_path,
                        const char *events_path) {
    if (!m || !ctx || !integ) return -1;

    memset(m, 0, sizeof(*m));
    m->sat   = sat;
    m->ctx   = ctx;
    m->integ = integ;
    m->t_start = t_start;
    m->t_end   = t_end;

    if (trajectory_path) {
        int rc = spody_setup_LogBuffer(&m->trajectory_log,
                                       trajectory_path,
                                       SPODY_LOG_BUFFER_CAPACITY,
                                       sizeof(TrajectoryRecord),
                                       SPODY_MAGIC_TRAJECTORY);
        if (rc != 0) {
            fprintf(stderr, "spody_setup_Mission: trajectory log setup failed (rc=%d)\n", rc);
            return rc;
        }
        m->log_trajectory = 1;
    }

    if (breakdown_path) {
        int rc = spody_setup_LogBuffer(&m->breakdown_log,
                                       breakdown_path,
                                       SPODY_LOG_BUFFER_CAPACITY,
                                       sizeof(ForceBreakdown),
                                       SPODY_MAGIC_BREAKDOWN);
        if (rc != 0) {
            fprintf(stderr, "spody_setup_Mission: breakdown log setup failed (rc=%d)\n", rc);
            spody_free_LogBuffer(&m->trajectory_log);
            m->log_trajectory = 0;
            return rc;
        }
        m->log_breakdown = 1;
    }

    if (events_path) {
        int rc = spody_setup_LogBuffer(&m->events_log,
                                       events_path,
                                       SPODY_LOG_BUFFER_CAPACITY,
                                       sizeof(EventRecord),
                                       SPODY_MAGIC_EVENTS);
        if (rc != 0) {
            fprintf(stderr, "spody_setup_Mission: events log setup failed (rc=%d)\n", rc);
            spody_free_LogBuffer(&m->trajectory_log);
            spody_free_LogBuffer(&m->breakdown_log);
            m->log_trajectory = 0;
            m->log_breakdown  = 0;
            return rc;
        }
        m->log_events = 1;
    }

    return 0;
}

int spody_free_Mission(Mission *m) {
    if (!m) return -1;
    int rc1 = 0, rc2 = 0, rc3 = 0;
    if (m->log_trajectory) rc1 = spody_free_LogBuffer(&m->trajectory_log);
    if (m->log_breakdown)  rc2 = spody_free_LogBuffer(&m->breakdown_log);
    if (m->log_events)     rc3 = spody_free_LogBuffer(&m->events_log);
    m->log_trajectory = 0;
    m->log_breakdown  = 0;
    m->log_events     = 0;
    if (rc1) return rc1;
    if (rc2) return rc2;
    return rc3;
}

/* ============================================================
 * Driver loop
 *
 * After each accepted step:
 *   - log trajectory (one record copy + buffer append, ~ns)
 *   - if breakdown log is enabled, re-evaluate the force decomposition
 *     on the current accepted (t, y) and log it. Cost ~1 RHS extra,
 *     ~12% on top of a single propagation step at N=100.
 * ============================================================ */
int spody_run_propagate(Mission *m) {
    if (!m || !m->integ) return -1;

    /* If the caller hasn't initialised t, set it now. */
    if (m->integ->t < m->t_start) {
        m->integ->t = m->t_start;
    }

    while (m->integ->t < m->t_end) {

        /* clip last step to land exactly on t_end */
        double remaining = m->t_end - m->integ->t;
        if (m->integ->h > remaining) m->integ->h = remaining;

        int rc = spody_propagate_onestep(m->integ);
        if (rc != 0) return rc;

        m->n_steps_done++;

        if (m->log_trajectory) {
            TrajectoryRecord r;
            r.t = m->integ->t;
            r.y[0] = m->integ->y[0]; r.y[1] = m->integ->y[1]; r.y[2] = m->integ->y[2];
            r.y[3] = m->integ->y[3]; r.y[4] = m->integ->y[4]; r.y[5] = m->integ->y[5];
            int wrc = spody_log_append(&m->trajectory_log, &r);
            if (wrc != 0) return wrc;
        }

        if (m->log_breakdown) {
            ForceBreakdown bd;
            spody_force_breakdown(m->ctx, m->integ->t, m->integ->y, &bd);
            int wrc = spody_log_append(&m->breakdown_log, &bd);
            if (wrc != 0) return wrc;
        }

        if (m->log_trajectory || m->log_breakdown) {
            m->n_records_saved++;
        }

        /* events: check after the step is logged so any trigger sees the
         * accepted state. A STOP action exits the loop after writing any
         * pending event record. */
        int stop_requested = 0;
        if (m->n_events > 0 && m->events) {
            for (int i = 0; i < m->n_events; i++) {
                SpodyEvent *ev = &m->events[i];
                int already = ev->triggered;
                /* Refined check: dense-output + Brent root finding inside the
                 * just-completed step. Falls back to coarse check if the
                 * integrator method does not provide dense_eval (e.g. RK4). */
                if (spody_event_check_refined(ev, m->ctx, m->integ) && !already) {
                    /* fired right now: bookkeeping + optional log + optional stop */
                    m->n_events_triggered++;

                    if (m->log_events &&
                        (ev->action == SPODY_EVENT_ACTION_LOG ||
                         ev->action == SPODY_EVENT_ACTION_LOG_AND_STOP)) {
                        EventRecord r;
                        r.t           = ev->t_trigger;
                        r.kind        = (int)ev->kind;
                        r.naif_id     = ev->naif_id;
                        r.radius_km   = ev->radius_km;
                        r.distance_km = ev->distance_at_trigger;
                        r.y[0] = ev->y_trigger[0]; r.y[1] = ev->y_trigger[1]; r.y[2] = ev->y_trigger[2];
                        r.y[3] = ev->y_trigger[3]; r.y[4] = ev->y_trigger[4]; r.y[5] = ev->y_trigger[5];
                        int wrc = spody_log_append(&m->events_log, &r);
                        if (wrc != 0) return wrc;
                    }

                    if (ev->action == SPODY_EVENT_ACTION_STOP ||
                        ev->action == SPODY_EVENT_ACTION_LOG_AND_STOP) {
                        stop_requested = 1;
                    }
                }
            }
        }
        if (stop_requested) return SPODY_RUN_STOPPED_BY_EVENT;
    }

    return SPODY_RUN_OK;
}

/* ============================================================
 * CSV converters (column names hard-coded for the two record types)
 * ============================================================ */
int spody_mission_trajectory_to_csv(const char *bin_path, const char *csv_path,
                                    int decimate) {
    static const char *cols[7] = {
        "t", "x", "y", "z", "vx", "vy", "vz"
    };
    FILE *out = fopen(csv_path, "w");
    if (!out) {
        perror("spody_mission_trajectory_to_csv fopen");
        return -1;
    }
    int rc = spody_log_dump_csv(bin_path, out, cols, 7, decimate);
    fclose(out);
    return rc;
}

/* Custom converter for the breakdown binary log: each record is one full
 * ForceBreakdown struct (mixed doubles + int n_third + padding + per-body
 * array), so the generic spody_log_dump_csv is not sufficient. We dump a
 * structured CSV with explicit column names for each field, including the
 * per-body slots up to SPODY_FM_MAX_THIRD. Unused per-body slots (i >=
 * n_third) are written as zeros (which is what spody_force_breakdown
 * leaves them, thanks to its memset). */
int spody_mission_breakdown_to_csv(const char *bin_path, const char *csv_path,
                                   int decimate) {
    if (decimate < 1) decimate = 1;

    FILE *bin = fopen(bin_path, "rb");
    if (!bin) {
        perror("spody_mission_breakdown_to_csv fopen bin");
        return -1;
    }

    SpodyLogHeader h;
    int rc = spody_log_read_header(bin, &h);
    if (rc != 0) {
        fprintf(stderr, "spody_mission_breakdown_to_csv: bad header (rc=%d)\n", rc);
        fclose(bin);
        return -2;
    }
    if (h.record_size != sizeof(ForceBreakdown)) {
        fprintf(stderr, "spody_mission_breakdown_to_csv: record_size=%u, expected %zu\n",
                h.record_size, sizeof(ForceBreakdown));
        fclose(bin);
        return -3;
    }

    FILE *out = fopen(csv_path, "w");
    if (!out) {
        perror("spody_mission_breakdown_to_csv fopen csv");
        fclose(bin);
        return -4;
    }

    /* CSV header */
    fputs("t,"
          "acc_total_x,acc_total_y,acc_total_z,"
          "acc_2body_x,acc_2body_y,acc_2body_z,"
          "acc_hg_x,acc_hg_y,acc_hg_z,"
          "acc_3rd_total_x,acc_3rd_total_y,acc_3rd_total_z,"
          "n_third", out);
    for (int i = 0; i < SPODY_FM_MAX_THIRD; i++) {
        fprintf(out, ",acc_3rd_%d_x,acc_3rd_%d_y,acc_3rd_%d_z", i, i, i);
    }
    fputs(",acc_srp_x,acc_srp_y,acc_srp_z,"
          "acc_drag_x,acc_drag_y,acc_drag_z,"
          "eclipse_fraction\n", out);

    /* records */
    long n_in = 0, n_out = 0;
    ForceBreakdown bd;
    while (fread(&bd, sizeof(bd), 1, bin) == 1) {
        if ((n_in % decimate) == 0) {
            fprintf(out, "%.18e,", bd.t);
            fprintf(out, "%.18e,%.18e,%.18e,",
                    bd.acc_total[0], bd.acc_total[1], bd.acc_total[2]);
            fprintf(out, "%.18e,%.18e,%.18e,",
                    bd.acc_2body[0], bd.acc_2body[1], bd.acc_2body[2]);
            fprintf(out, "%.18e,%.18e,%.18e,",
                    bd.acc_sphericalharmonics[0],
                    bd.acc_sphericalharmonics[1],
                    bd.acc_sphericalharmonics[2]);
            fprintf(out, "%.18e,%.18e,%.18e,",
                    bd.acc_thirdbody_total[0],
                    bd.acc_thirdbody_total[1],
                    bd.acc_thirdbody_total[2]);
            fprintf(out, "%d", bd.n_third);
            for (int i = 0; i < SPODY_FM_MAX_THIRD; i++) {
                fprintf(out, ",%.18e,%.18e,%.18e",
                        bd.acc_thirdbody[i][0],
                        bd.acc_thirdbody[i][1],
                        bd.acc_thirdbody[i][2]);
            }
            fprintf(out, ",%.18e,%.18e,%.18e",
                    bd.acc_srp[0], bd.acc_srp[1], bd.acc_srp[2]);
            fprintf(out, ",%.18e,%.18e,%.18e",
                    bd.acc_drag[0], bd.acc_drag[1], bd.acc_drag[2]);
            fprintf(out, ",%.18e\n", bd.eclipse_fraction);
            n_out++;
        }
        n_in++;
    }

    fclose(out);
    fclose(bin);
    fprintf(stderr, "spody_mission_breakdown_to_csv: %ld records read, %ld rows written (decimate=%d)\n",
            n_in, n_out, decimate);
    return 0;
}
