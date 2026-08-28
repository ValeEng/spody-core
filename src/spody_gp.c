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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */
#include <math.h>
#include <stdio.h>

#include "spody_gp.h"
#include "spody_const.h"
#include "spody_earth_orientation.h"
#include "spody_eop.h"
#include "spody_forcemodels.h"
#include "spody_time.h"

int spody_gp_to_state_icrf(const spody_sgp4_elements *el,
                           const char *eop_file,
                           const char *iau2006_dir,
                           double *et_s,
                           double r_icrf_km[3],
                           double v_icrf_kms[3])
{
    if (!el || !eop_file || !iau2006_dir || !et_s || !r_icrf_km || !v_icrf_kms) {
        return 1;
    }

    spody_sgp4_sat sat;
    int rc = spody_sgp4_init(el, &sat);
    if (rc != SPODY_SGP4_OK) {
        fprintf(stderr, "gp: element set rejected by sgp4, rc=%d\n", rc);
        return 1;
    }

    double r_teme[3], v_teme[3];
    rc = spody_sgp4_at(&sat, 0.0, r_teme, v_teme);
    if (rc != SPODY_SGP4_OK) {
        fprintf(stderr, "gp: sgp4 failed at epoch, rc=%d\n", rc);
        return 1;
    }

    /* The epoch is anchored on midnight and the day fraction is added
     * afterwards, rather than scaling the whole MJD at once. The ulp of
     * a modern Julian Date is some 40 us, which is 30 cm along a LEO
     * track, and the half-integer below keeps the day exact so only the
     * fraction carries rounding. Then the timescale chain, one step at
     * a time: leap seconds standing at this date, the fixed 32.184 to
     * TT, the periodic offset to TDB. */
    double mjd_day  = floor(el->epoch_mjd);
    double mjd_frac = el->epoch_mjd - mjd_day;
    double utc_sec  = (mjd_day - (JD_J2000 - JD_MJD_EPOCH)) * SECONDSxDAY
                    + mjd_frac * SECONDSxDAY;
    double tt_sec   = utc_sec + spody_tai_minus_utc(el->epoch_mjd) - TT2TAI_SEC;
    double et       = tt_sec + spody_tdb_minus_tt(tt_sec);

    MappedEOPData eop_data = {0};
    if (spody_setup_MappedEOPData(&eop_data, eop_file) != 0) {
        fprintf(stderr, "gp: cannot load EOP from '%s'\n", eop_file);
        return 1;
    }
    MappedIAU2006Data iau_data = {0};
    if (spody_setup_MappedIAU2006Data(&iau_data, iau2006_dir) != 0) {
        fprintf(stderr, "gp: cannot load IAU 2006 from '%s'\n", iau2006_dir);
        spody_free_MappedEOPData(&eop_data);
        return 1;
    }
    MappedEOP     eop_map = {0};
    MappedIAU2006 iau_map = {0};
    spody_setup_MappedEOP(&eop_map, &eop_data);
    spody_setup_MappedIAU2006(&iau_map, &iau_data);
    ForceModelContext ctx = { .eop = &eop_map, .iau2006 = &iau_map };

    double R_teme2i[3][3], R_i2teme[3][3];
    spody_teme2icrf_rotation(&ctx, et, R_teme2i, R_i2teme);

    for (int a = 0; a < 3; ++a) {
        r_icrf_km[a]  = R_teme2i[a][0] * r_teme[0]
                      + R_teme2i[a][1] * r_teme[1]
                      + R_teme2i[a][2] * r_teme[2];
        /* TEME does not turn with the Earth, so the same matrix carries
         * the velocity and there is no omega x r term. */
        v_icrf_kms[a] = R_teme2i[a][0] * v_teme[0]
                      + R_teme2i[a][1] * v_teme[1]
                      + R_teme2i[a][2] * v_teme[2];
    }
    *et_s = et;

    spody_free_MappedIAU2006(&iau_map);
    spody_free_MappedEOP(&eop_map);
    spody_free_MappedIAU2006Data(&iau_data);
    spody_free_MappedEOPData(&eop_data);
    return 0;
}
