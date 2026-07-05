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
#include "spody_ephemeris.h"

static int read_ephemeris_file_header(FILE *file, EphemerisFile_Header *ep ) {
    //parsing for de440

    int found_group = 0; 
    int data_row = 0;    
    char line[BUFFER_SIZE_EPH];

    while (fgets(line, sizeof(line), file)) {
        
        // looking for GROUP 1050 
        if (!found_group) {
            if (strstr(line, "GROUP") || strstr(line, "KSIZE")){
                if(strstr(line, "KSIZE")) {
                    #if DEBUG_EPHEMERIS == 1
                    printf("RECORD SIZE found (normally the first row)\n");
                    #endif
                    char *token = strtok(line, " \t\r\n");
                    while (token) {
                        if (strstr(token, "KSIZE")) {
                            token = strtok(NULL, " \t\r\n");
                            ep->bytes_per_record = strtod(token, NULL);
                            #if DEBUG_EPHEMERIS == 1
                            printf("load KSIZE : %d\n",ep->bytes_per_record);
                            #endif
                        } else if (strstr(token, "NCOEFF")) {
                            token = strtok(NULL, " \t\r\n");
                            ep->number_coefficients_per_record = strtod(token, NULL); ;
                            #if DEBUG_EPHEMERIS == 1
                            printf("load NCOEFF : %d\n",ep->number_coefficients_per_record);
                            #endif
                        }
                        token = strtok(NULL, " \t\r\n");
                    }
                    #if DEBUG_EPHEMERIS == 1
                    printf("RECORD SIZE %d | N COEFF %d\n", ep->bytes_per_record, ep->number_coefficients_per_record);
                    #endif
                    continue;
                }else if (strstr(line, "1030")) {
                    found_group = 1030;
                    #if DEBUG_EPHEMERIS == 1
                    printf("GROUP 1030 found\n");
                    #endif
                }else if (strstr(line, "1050")) {
                    found_group = 1050;
                    #if DEBUG_EPHEMERIS == 1
                    printf("GROUP 1050 found\n");
                    #endif
                }
            }
            continue; //next line 
        }

        char *token = strtok(line, " \t\r\n"); 
        
        if (token == NULL) {
            continue; // skip empty line post GROUP 1050 
        }
        if (found_group == 1030){
            int ele = 0;
            while (token) {
                double val = strtod(token, NULL);

                if (ele == 0) {
                    ep->start_epoch = val;
                } else if (ele == 1) {
                    ep->end_epoch = val;
                } else if (ele == 2) {
                    /* ASCII gives days; we keep the value here temporarily
                     * and convert to seconds (and ET for start/end) just
                     * before writing the binary header. */
                    ep->seconds_per_record = (int)val;
                }

                #if DEBUG_EPHEMERIS == 1
                printf("Ele %d: %f\n", ele, val);
                #endif
                ele++;
                token = strtok(NULL, " \t\r\n");
            }
            found_group = 0;

        }else if (found_group == 1050){
            int ele = 0;
            while (token) {
                double val = strtod(token, NULL);

                if (data_row == 0) {
                    ep->location[ele] = val;
                } else if (data_row == 1) {
                    ep->number_coefficients_per_component[ele] = val;
                } else if (data_row == 2) {
                    ep->number_complete_sets_coefficients_per_record[ele] = val;
                }

                #if DEBUG_EPHEMERIS == 1
                printf("row %d, Ele %d: %f\n", data_row, ele, val);
                #endif

                ele++;
                token = strtok(NULL, " \t\r\n");
            }
            // count for the 3 rows of data that we need
            data_row++;
            if (data_row > 2) {
                found_group = 0;
                data_row = 0;
                //break;
            }
        }
    }
    
    return 0;
}

static double fds2cd(char *str) {
    //parse fortran "double" string to c double 

    char temp[64];
    strncpy(temp, str, 63);
    temp[63] = '\0';
    
    //change from D (or d) to E
    char *p = strchr(temp, 'D');
    if (!p) p = strchr(temp, 'd');
    if (p) *p = 'E';
    
    return strtod(temp, NULL);
}

static int read_record_block(FILE *fp, EphemerisFile_Header *ep, EphemerisFile_Record *eprec) {
   
    char line[BUFFER_SIZE_EPH];
    int n_coeff_expected = 0;
    #if DEBUG_EPHEMERIS == 1
    printf("in read record\n");
    #endif

    if (!fgets(line, sizeof(line), fp)) return 0; // it's end of file
    //printf("Reading line: %s\n", line);
    if (sscanf(line, "%d %d", &eprec->record_number, &eprec->number_coefficients_per_record) != 2) return 0; // Parsing of number of coefficients
    //printf("Reading record %d with %d coefficients\n", eprec->record_number, eprec->number_coefficients_per_record);
    if (eprec->number_coefficients_per_record != ep->number_coefficients_per_record) {
        /* Real parser error -- left unconditional so the user is told
         * the on-disk ASCII is inconsistent with its own header. */
        printf("Warning: Expected %d coefficients, but found %d in block %d\n", n_coeff_expected, eprec->number_coefficients_per_record, eprec->record_number);
        printf("HEADER info diverge from data blocks read.\n");
        printf("Check the ephemeris file integrity.\n");
        return 0;
    }else{
        n_coeff_expected = eprec->number_coefficients_per_record;
        #if DEBUG_EPHEMERIS == 1
        printf("n_coeff_expected : %d\n",n_coeff_expected);
        #endif
    }

    #if DEBUG_EPHEMERIS == 1
    printf("cleaning ok\n");
    #endif

    int coeff_idx = 0;
    char *token;
    while (coeff_idx < n_coeff_expected) {
        if (fgets(line, sizeof(line), fp)){
            #if DEBUG_EPHEMERIS == 1
            printf("Reading coefficients line: %s\n", line);
            #endif
            token = strtok(line, " \t\r\n");
            while (token != NULL && coeff_idx < n_coeff_expected) {
                //printf("Token: %s\n", token);
                eprec->record[coeff_idx] = fds2cd(token);
                coeff_idx++;
                token = strtok(NULL, " \t\r\n");
            }
        }else{
            return -100; //end file
        }
    }
    #if DEBUG_EPHEMERIS == 1
    printf("token loaded : %d\n",coeff_idx);
    #endif
    /* Convert the JD time bracket from the ASCII source to ET (seconds past
     * J2000 TDB). The on-disk binary format (SPDEET) keeps record epochs in
     * ET so the runtime never has to convert back. */
    eprec->start_epoch = ET_FROM_JD(eprec->record[0]);
    eprec->end_epoch   = ET_FROM_JD(eprec->record[1]);

    return 1; //good
}

static int create_binary_ephemeris_file(EphemerisFile_Header *ep, int64_t *old_epoch, const char *ascp_filename, const char *bin_filename) {
    #if DEBUG_EPHEMERIS == 1
    printf("in create binary\n");
    #endif
    FILE *fp_ascp = fopen(ascp_filename, "r");
    if (!fp_ascp) { perror("cannot open ascp file"); return -1; }

    FILE *fp_bin = fopen(bin_filename, "ab");  // "ab" for append mode
    if (!fp_bin) { perror("cannot open bin file"); fclose(fp_ascp); return -1; }
    #if DEBUG_EPHEMERIS == 1
    printf("file loaded\n");
    #endif

    double t_start, t_end;
    size_t record_size = ep->bytes_per_record;
    #if DEBUG_EPHEMERIS == 1
    printf("sizeof(EphemerisFile_Record) : %zu\n", sizeof(EphemerisFile_Record));
    printf("ep->number_coefficients_per_record * sizeof(double) : %zu\n",ep->number_coefficients_per_record * sizeof(double));
    printf("read -> bytes_per_record : %zu\n",record_size);
    #endif

    EphemerisFile_Record *eprec = malloc(record_size);
    if (!eprec) {
        perror("Error memory allocation for EphemerisFile_Record");
        fclose(fp_ascp);
        fclose(fp_bin);
        return -1;
    }
    #if DEBUG_EPHEMERIS == 1
    printf("malloc ok\n");
    #endif

    int block_count = 0;
    // Read each record block and write to binary file
    while (read_record_block(fp_ascp, ep, eprec)) {
        block_count++;
        #if DEBUG_EPHEMERIS == 1
        printf("block %d readed\n",block_count);
        printf("old epoch : %d\n",*old_epoch);
        #endif

        /* Detect duplicate records across consecutive ASCII chunks. ET values
         * for adjacent records differ by at least seconds_per_record (~32 days
         * = 2.7e6 s), so an int64 truncation is more than safe to compare. */
        if (*old_epoch == (int64_t)eprec->start_epoch) {
            /* Real diagnostic about overlapping chunks -- left unconditional
             * so the user knows we silently dropped a record. */
            printf("\n\n! a clone record foud ! record %d start date %.2f \n\n",eprec->record_number, eprec->start_epoch);
            continue;
        }
        *old_epoch = (int64_t)eprec->start_epoch;
        #if DEBUG_EPHEMERIS == 1
        printf("new old epoch : %lld\n", (long long)*old_epoch);
        printf("Writing record %d to binary file\n", eprec->record_number);
        printf("size of : %zu\n",record_size);
        #endif
        fwrite(eprec, record_size, 1, fp_bin);
        #if DEBUG_EPHEMERIS == 1
        printf("Written record %04d with %04d coefficients from %.6f to %.6f\n", eprec->record_number, eprec->number_coefficients_per_record, eprec->start_epoch, eprec->end_epoch);
        #endif

    }
    free(eprec);
    fclose(fp_ascp);
    fclose(fp_bin);
    return 0;
}

/**
 * Evaluate a Chebyshev polynomial at a given point.
 *
 * Uses the Clenshaw-Smith method (a modified Horner scheme) for efficient evaluation.
 *
 * @param time_scaled Normalized and scaled time in the range [-1.0, 1.0].
 * @param coefficients Array of Chebyshev coefficients.
 * @param n_coeffs Number of coefficients (polynomial order).
 * @return The polynomial value at the given time_scaled.
 */

static double chebyshev_evaluate(double time_scaled, const double *coefficients, int n_coeffs) {
    if (n_coeffs <= 0) return 0.0;
    if (n_coeffs == 1) return coefficients[0];

    // t_scaled = x --> T_n(x)
    double x = time_scaled;
    double two_x = 2.0 * x;

    // v_k initialization
    double v_kp1 = 0.0; // v_{k+1}
    double v_k = 0.0;   // v_k
    double v_km1;       // v_{k-1}

    // backward sum
    for (int k = n_coeffs - 1; k >= 1; k--) {
        // v_{k-1} = c_k + 2x * v_k - v_{k+1}
        v_km1 = coefficients[k] + (two_x * v_k) - v_kp1;
        #if DEBUG_CHEBYSHEV == 1
        printf("k=%02d\t c_k=%+.21e\t v_km1=%+.21e\n", k, coefficients[k], v_km1);
        #endif
        // update for next iter 
        v_kp1 = v_k;
        v_k = v_km1;
    }

    // final result --> c_0 + x * v_0 - v_1.
    // in the loop, v_k ---> v_0 | v_kp1 ---> v_1.
    // return is for k=0 ---> c_0 + x*v_0 - v_1  
    // that is the standard formula for c_0 + x*v_k - v_{k+1}, k=0.
    return coefficients[0] + (x * v_k) - v_kp1;
}

/**
 * Compute the position (X, Y, Z) of a celestial body at a given time.
 *
 * @param target_idx Index of the celestial body (0=Mercury, 1=Venus, ...) respect to @param ef.
 * @param et Target Ephemeris Time (seconds past J2000 TDB).
 * @param ef Ephemeris JPL file structure (from ASCII Header).
 * @param coeffs Array of coefficients from the ASCP file.
 * @param t_start Julian Date at the start of the current block.
 * @param t_end Julian Date at the end of the current block.
 * @param result Output array for the position [X, Y, Z].
 * @return 1 if successful, 0 if the date is out of range or data is invalid.
 */

static int calculate_body_position(MappedEphemeris *map, int target_idx, double et, double result[3]) {

    // cache hit: same body, same epoch -> skip record/set lookup and Chebyshev
    if (target_idx >= 0 && target_idx < EPH_CACHE_SLOTS &&
        map->cache_valid[target_idx] && map->cache_jd[target_idx] == et) {
        result[0] = map->cache_pos[target_idx][0];
        result[1] = map->cache_pos[target_idx][1];
        result[2] = map->cache_pos[target_idx][2];
        return 1;
    }

    const MappedEphemerisData *med = map->med;

    // 1. Retrieve body parameters
    int numeber_coefficients_per_component = med->header->number_coefficients_per_component[target_idx];
    int start_index = med->header->location[target_idx] - 1; // 1-based -> 0-based
    int n_components = 3; // X, Y, Z

    // 2. Subdivision identification (all values already in ET seconds since
    //    the SPDEET file format is ET-native: header epochs and per-record
    //    epochs are stored as `seconds past J2000 TDB`). Subtractions stay
    //    between magnitudes ~1e9 with ULP ~200 ns, ~250x better than JD.
    int record_id = (int)floor((et - med->header->start_epoch) / med->header->seconds_per_record);

    EphemerisFile_Record *eprec = med->records[record_id];
    double rec_start = eprec->start_epoch;
    double rec_end   = eprec->end_epoch;
    double block_duration = rec_end - rec_start;
    int n_sets = med->header->number_complete_sets_coefficients_per_record[target_idx];
    double set_duration = block_duration / (double)n_sets;
    int set_id = (int)floor((et - rec_start) / set_duration);

    // floating point error check
    if (set_id >= n_sets) {
        set_id = n_sets - 1;
        printf("Adjusted Set ID due to floating point error: %d\n", set_id);
    }

    // 3. Tau evaluation
    double t_gran_start = rec_start + set_id * set_duration;
    double t_gran_end   = t_gran_start + set_duration;

    // normalizing of tau --> [-1.0, 1.0]
    double tau = (2.0 * et - t_gran_start - t_gran_end) / (t_gran_end - t_gran_start);

    // 4. Buffer param
    // a complete set (X, Y, Z) per set_id is 3 * n_coeffs.
    int set_length = n_components * numeber_coefficients_per_component;
    int offset = start_index + (set_id * set_length);

    // 5. Polynomial evaluation for X, Y, Z
    for (int i = 0; i < n_components; i++) {
        const double *coeff_ptr = eprec->record + offset + (i * numeber_coefficients_per_component);
        result[i] = chebyshev_evaluate(tau, coeff_ptr, numeber_coefficients_per_component);
    }

    // store in cache (cache_jd field is now keyed by ET, name kept for layout stability)
    if (target_idx >= 0 && target_idx < EPH_CACHE_SLOTS) {
        map->cache_jd[target_idx] = et;
        map->cache_pos[target_idx][0] = result[0];
        map->cache_pos[target_idx][1] = result[1];
        map->cache_pos[target_idx][2] = result[2];
        map->cache_valid[target_idx] = 1;
    }

    return 1;
}

//MAPPING FUNCTIONS--------------------------------------------------------------------------------------------------

static int ephemeris_map_file(MappedEphemerisData *med, const char *filename) {

    if (!med) return -1; //if NULL exit

    if (mf_map_file(&med->mf, filename) != 0) return -2; //mapping error

    /* Private heap copy of the header: the file mapping is read-only
     * (PAGE_READONLY / FILE_MAP_READ), and subset files need their
     * coverage epochs reconciled below -- patching the mapped bytes
     * would fault. Records stay pointers into the read-only map. */
    med->header = malloc(sizeof *med->header);
    if (!med->header) { mf_unmap_file(&med->mf); return -3; }
    memcpy(med->header, med->mf.ptr, sizeof *med->header);

    #if DEBUG_EPHEMERIS == 1
    printf("mf_size : %zu\n",med->mf.size);
    printf("sizeof(EphemerisFile_Header) : %zu\n",sizeof(EphemerisFile_Header));
    #endif

    /* Validate the on-disk format: magic SPDEET + supported version. */
    if (memcmp(med->header->magic, SPODY_EPH_MAGIC_ET, SPODY_EPH_MAGIC_LEN) != 0) {
        printf("ephemeris_map_file: bad magic (got '%.8s', expected '%.8s'). "
               "Regenerate the .spody binary with the current spody_createfile_*.\n",
               med->header->magic, SPODY_EPH_MAGIC_ET);
        free(med->header); med->header = NULL;
        mf_unmap_file(&med->mf);
        return -10;
    }
    if (med->header->format_version != SPODY_EPH_FORMAT_VERSION) {
        printf("ephemeris_map_file: unsupported format_version %u (expected %u)\n",
               (unsigned)med->header->format_version,
               (unsigned)SPODY_EPH_FORMAT_VERSION);
        free(med->header); med->header = NULL;
        mf_unmap_file(&med->mf);
        return -11;
    }

    size_t remaining_bytes = med->mf.size - sizeof(EphemerisFile_Header);
    #if DEBUG_EPHEMERIS == 1
    printf("remaning bytes : %zu \n",remaining_bytes);
    #endif

    med->num_records = remaining_bytes / med->header->bytes_per_record; //we hope it is exact division

    #if DEBUG_EPHEMERIS == 1
    printf("bytes per record: %d\n", med->header->bytes_per_record);
    printf("Number of records mapped: %zu\n", med->num_records);
    #endif

    med->records = (EphemerisFile_Record**)malloc(sizeof(EphemerisFile_Record*) * med->num_records);
    if (!med->records) {
        free(med->header); med->header = NULL;
        mf_unmap_file(&med->mf);
        return -3;
    }

    uint8_t *ptr = (uint8_t*)med->mf.ptr + sizeof(EphemerisFile_Header);
    for (size_t i = 0; i < med->num_records; i++) {
        med->records[i] = (EphemerisFile_Record*)ptr;
        #if DEBUG_EPHEMERIS == 1
        printf("Mapped record %zu: record number %d, n_coeff %d, start_epoch %.6f, end_epoch %.6f\n", i+1, med->records[i]->record_number, med->records[i]->number_coefficients_per_record, med->records[i]->start_epoch, med->records[i]->end_epoch);
        #endif
        ptr += med->header->bytes_per_record;
    }

    /* Subset files (e.g. a partial DE440 conversion covering only the
     * chunks the user downloaded) may carry the full-span epochs the
     * converter read from header.440 before knowing which chunks it
     * would write. The records are the truth: reconcile the private
     * header copy so the record-index arithmetic in
     * get_body_position() stays exact. */
    if (med->num_records > 0) {
        double rec_start = med->records[0]->start_epoch;
        double rec_end   = med->records[med->num_records - 1]->end_epoch;
        if (med->header->start_epoch != rec_start ||
            med->header->end_epoch   != rec_end) {
            printf("ephemeris: subset file -- header claims %.3f .. %.3f ET "
                   "but records cover %.3f .. %.3f ET; using the records' "
                   "range.\n",
                   med->header->start_epoch, med->header->end_epoch,
                   rec_start, rec_end);
            med->header->start_epoch = rec_start;
            med->header->end_epoch   = rec_end;
        }
    }

    /* No runtime conversion needed: epochs in the file are already ET. */
    return 0;
}

static int ephemeris_unmap_file(MappedEphemerisData *med) {
    if (!med) return -1;

    free(med->records); med->records = NULL;
    /* The header is always a private heap copy (see ephemeris_map_file
     * and spody_setup_partialMappedEphemerisData). */
    free(med->header); med->header = NULL;
    med->num_records = 0;

    return mf_unmap_file(&med->mf);
}

int spody_createfile_MappedEphemerisData(const char *path, const char **file_names, const int n_files, const char *de){
    
    char header_path[1000]; //TBD 
    char bin_filename[1000];
    char ascp_filename[1000];

    int returnNumber;

    sprintf(header_path, "./%s/header.%s",path,de); 
    sprintf(bin_filename, "./%s/de%s.spody",path,de); 

    FILE *file = fopen(header_path,"r");
    if (!file) { perror("cannot open header file"); return -1; }

    /* Truncate the destination on first open: avoids accidental append
     * to a previous run. Subsequent record writes use "ab". */
    FILE *fp_bin = fopen(bin_filename, "wb");
    if (!fp_bin) { perror("cannot open bin file"); fclose(file); return -1; }

    EphemerisFile_Header ep = {0};
    returnNumber = read_ephemeris_file_header(file, &ep);
    fclose(file);

    size_t n = ep.number_coefficients_per_record;
    size_t record_size = sizeof(EphemerisFile_Record) + n * sizeof(double);
    ep.bytes_per_record = (int)record_size;
    #if DEBUG_EPHEMERIS == 1
    printf("ep.bytes_per_record : %d\n",ep.bytes_per_record);
    #endif

    /* Convert the header epochs from the ASCII (JD, days) source to the
     * SPDEET on-disk format (ET, seconds). */
    memcpy(ep.magic, SPODY_EPH_MAGIC_ET, SPODY_EPH_MAGIC_LEN);
    ep.format_version    = SPODY_EPH_FORMAT_VERSION;
    ep.reserved          = 0;
    ep.start_epoch       = ET_FROM_JD(ep.start_epoch);
    ep.end_epoch         = ET_FROM_JD(ep.end_epoch);
    /* days -> s. Explicit cast: the field is int and SECONDSxDAY is a
     * double; the product is an exactly-representable integer. */
    ep.seconds_per_record = (int)(ep.seconds_per_record * SECONDSxDAY);

    /* first write header info */
    fwrite(&ep, sizeof(EphemerisFile_Header), 1, fp_bin);
    fclose(fp_bin);

    #if DEBUG_EPHEMERIS == 1
    printf("header writed, size : %zu (magic=%.8s, version=%u)\n",
           sizeof(EphemerisFile_Header), ep.magic, (unsigned)ep.format_version);
    #endif

    int64_t old_epoch = 1; /* necessary to avoid duplicate; ET values fit easily in int64 */

    for (int i = 0; i < n_files; i++){

        sprintf(ascp_filename, "./%s/ascp%s.%s",path,file_names[i],de);

        returnNumber = create_binary_ephemeris_file(&ep,&old_epoch,ascp_filename,bin_filename);
        #if DEBUG_EPHEMERIS == 1
        printf("old epoch : %lld\n",(long long)old_epoch);
        #endif

        /* One-line progress per file -- kept unconditional so the user
         * (and the GUI's wizard convert window) sees the loop tick. */
        printf("file %d writed\n",i);

    }

    /* The header was written before any chunk was converted, carrying
     * the full-DE440 epoch span read from header.440. When only a
     * subset of the ASCII chunks is converted (e.g. the GUI wizard's
     * modern-era profile), refresh the on-disk epochs from the records
     * actually written so the file is self-consistent. */
    if (returnNumber == 0) {
        FILE *fp_fix = fopen(bin_filename, "rb+");
        if (fp_fix) {
            EphemerisFile_Header hdr;
            long file_size = 0;
            if (fread(&hdr, sizeof hdr, 1, fp_fix) == 1 &&
                fseek(fp_fix, 0, SEEK_END) == 0 &&
                (file_size = ftell(fp_fix)) > (long)sizeof hdr &&
                hdr.bytes_per_record > 0) {
                long n_rec = (file_size - (long)sizeof hdr)
                             / hdr.bytes_per_record;
                long first_off = (long)sizeof hdr
                    + (long)offsetof(EphemerisFile_Record, start_epoch);
                long last_off  = (long)sizeof hdr
                    + (n_rec - 1) * (long)hdr.bytes_per_record
                    + (long)offsetof(EphemerisFile_Record, end_epoch);
                double first_start = 0.0, last_end = 0.0;
                if (n_rec > 0 &&
                    fseek(fp_fix, first_off, SEEK_SET) == 0 &&
                    fread(&first_start, sizeof first_start, 1, fp_fix) == 1 &&
                    fseek(fp_fix, last_off, SEEK_SET) == 0 &&
                    fread(&last_end, sizeof last_end, 1, fp_fix) == 1 &&
                    (hdr.start_epoch != first_start ||
                     hdr.end_epoch   != last_end)) {
                    hdr.start_epoch = first_start;
                    hdr.end_epoch   = last_end;
                    if (fseek(fp_fix, 0, SEEK_SET) == 0 &&
                        fwrite(&hdr, sizeof hdr, 1, fp_fix) == 1) {
                        printf("header epochs refreshed to the converted "
                               "range (%.3f .. %.3f ET)\n",
                               first_start, last_end);
                    }
                }
            }
            fclose(fp_fix);
        }
    }

    return returnNumber;

}

int spody_setup_MappedEphemerisData(MappedEphemerisData *med, const char *filename){
    /* Subset-coverage reconciliation happens inside ephemeris_map_file
     * (on the private header copy) so the partial loader gets healed
     * epochs too. */
    return ephemeris_map_file(med, filename);
}

int spody_setup_MappedEphemeris(MappedEphemeris *map, const MappedEphemerisData *med){
    if (!map || !med) return -1;
    map->med = med;
    for (int i = 0; i < EPH_CACHE_SLOTS; i++) map->cache_valid[i] = 0;
    return 0;
}

int spody_free_MappedEphemeris(MappedEphemeris *map){
    if (!map) return -1;
    map->med = NULL;
    for (int i = 0; i < EPH_CACHE_SLOTS; i++) map->cache_valid[i] = 0;
    return 0;
}

int spody_free_MappedEphemerisData(MappedEphemerisData *med){
    return ephemeris_unmap_file(med);
}

/*****************NAIF ID********************
    SSB     = 0
    MERCURY = 1 or 199
    VENUS   = 2 or 299
    EMB     = 3
    MARS    = 4 or 499
    JUPITER = 5 or 599
    SATURN  = 6 or 699
    URANUS  = 7 or 799
    NEPTUNE = 8 or 899
    PLUTO   = 9 or 999
    SUN     = 10
    MOON    = 301
    EARTH   = 399
    (planet barycenter and planet center collapse to the same record in DE440
     except for Earth/Moon, handled explicitly via EMRAT)
********************************************/

static int get_body_position_ssb(MappedEphemeris *map, int naif_id, double et, double result[3]){
    double temp[3];
    switch (naif_id) {
    case 0:
        result[0] = 0.0; result[1] = 0.0; result[2] = 0.0;
        return 0;
    case 1: case 199: return calculate_body_position(map, 0,  et, result);
    case 2: case 299: return calculate_body_position(map, 1,  et, result);
    case 3:           return calculate_body_position(map, 2,  et, result);
    case 4: case 499: return calculate_body_position(map, 3,  et, result);
    case 5: case 599: return calculate_body_position(map, 4,  et, result);
    case 6: case 699: return calculate_body_position(map, 5,  et, result);
    case 7: case 799: return calculate_body_position(map, 6,  et, result);
    case 8: case 899: return calculate_body_position(map, 7,  et, result);
    case 9: case 999: return calculate_body_position(map, 8,  et, result);
    case 10:          return calculate_body_position(map, 10, et, result);
    case 399: {
        // Earth_ssb = EMB_ssb - 1/(1+EMRAT) * r_moon_earth
        calculate_body_position(map, 2, et, result);
        calculate_body_position(map, 9, et, temp);
        double f = -1.0 / (1.0 + EMRAT);
        result[0] += f * temp[0];
        result[1] += f * temp[1];
        result[2] += f * temp[2];
        return 1;
    }
    case 301: {
        // Moon_ssb = EMB_ssb + EMRAT/(1+EMRAT) * r_moon_earth
        calculate_body_position(map, 2, et, result);
        calculate_body_position(map, 9, et, temp);
        double f = EMRAT / (1.0 + EMRAT);
        result[0] += f * temp[0];
        result[1] += f * temp[1];
        result[2] += f * temp[2];
        return 1;
    }
    default:
        result[0] = 0.0; result[1] = 0.0; result[2] = 0.0;
        return -1;
    }
}

int spody_get_ephposition(MappedEphemeris *map, int central_idx, int target_idx, double et, double result[3]){

    // fast path: Earth <-> Moon uses a single Chebyshev evaluation
    if (central_idx == 399 && target_idx == 301) {
        return calculate_body_position(map, 9, et, result);
    }
    if (central_idx == 301 && target_idx == 399) {
        calculate_body_position(map, 9, et, result);
        result[0] = -result[0];
        result[1] = -result[1];
        result[2] = -result[2];
        return 0;
    }

    double central[3];
    if (get_body_position_ssb(map, target_idx, et, result) < 0) {
        printf("Target body not supported\n");
        return -1;
    }
    if (get_body_position_ssb(map, central_idx, et, central) < 0) {
        printf("Central body not supported\n");
        result[0] = 0.0; result[1] = 0.0; result[2] = 0.0;
        return -1;
    }
    result[0] -= central[0];
    result[1] -= central[1];
    result[2] -= central[2];
    return 0;
}

int spody_get_ephposition_batch(MappedEphemeris *map, int central_idx, const int *target_idx_array, int n_targets, double et, double *result){
    // Flat buffer layout: result[3*i + 0..2] is (x,y,z) for target i.
    // Central body is SSB-reduced once; all targets reuse it.
    // Further deduplication (e.g. EMB+Moon_geo shared across Earth/Moon requests)
    // is handled automatically by the per-body cache in calculate_body_position.
    double central[3];
    if (get_body_position_ssb(map, central_idx, et, central) < 0) {
        printf("Central body not supported\n");
        for (int i = 0; i < 3 * n_targets; i++) result[i] = 0.0;
        return -1;
    }

    for (int i = 0; i < n_targets; i++) {
        double target_ssb[3];
        double *out = result + 3 * i;
        if (get_body_position_ssb(map, target_idx_array[i], et, target_ssb) < 0) {
            printf("Target body %d not supported\n", target_idx_array[i]);
            out[0] = 0.0; out[1] = 0.0; out[2] = 0.0;
            continue;
        }
        out[0] = target_ssb[0] - central[0];
        out[1] = target_ssb[1] - central[1];
        out[2] = target_ssb[2] - central[2];
    }
    return 0;
}

int spody_get_lunarlibrationangles(MappedEphemeris *map, double et, double result[3]){
    
    calculate_body_position(map, 12, et, result);

    return 0; 

}

void spody_getrotmatrix_icrf2moonpa(double phi, double theta, double psi, double R[3][3]){
    const double cphi = cos(phi);
    const double sphi = sin(phi);
    const double cth  = cos(theta);
    const double sth  = sin(theta);
    const double cpsi = cos(psi);
    const double spsi = sin(psi);

    R[0][0] =  cpsi*cphi - spsi*cth*sphi;
    R[0][1] =  cpsi*sphi + spsi*cth*cphi;
    R[0][2] =  spsi*sth;

    R[1][0] = -spsi*cphi - cpsi*cth*sphi;
    R[1][1] = -spsi*sphi + cpsi*cth*cphi;
    R[1][2] =  cpsi*sth;

    R[2][0] =  sth*sphi;
    R[2][1] = -sth*cphi;
    R[2][2] =  cth;
}

void spody_getrotmatrix_moonpa2icrf(double phi, double theta, double psi, double R[3][3]){
    
    double Cfwd[3][3];
    spody_getrotmatrix_icrf2moonpa(phi, theta, psi, Cfwd);

    /* transpose */
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R[i][j] = Cfwd[j][i];
}


//-----Numeric error safer function ---------------------------------------------------------------

/*------only JD without the subtraction of in_start

int spody_setup_partialMappedEphemeris(MappedEphemeris *map, const char *filename, double in_start, double in_end){

    //TBD necessary a free function but now we end the program for memory free

    MappedEphemeris full = {0};
    if (ephemeris_map_file(&full, filename) != 0) return -1;

    int record_id_start = (int)floor((in_start - full.header->start_epoch)/full.header->days_per_record);
    int record_id_end = (int)floor((in_end - full.header->start_epoch)/full.header->days_per_record);
    int number_of_records = (record_id_end - record_id_start) + 1; // is alway + 1 wrt the difference 
    if (number_of_records <= 0) return -4;

    printf("record_id_start : %d | record_id_end : %d \nnumber_of_reccords : %d\n", record_id_start, record_id_end, number_of_records);

    map->header = malloc(sizeof(EphemerisFile_Header));
    *map->header = *full.header;

    map->records = (EphemerisFile_Record**)malloc(sizeof(EphemerisFile_Record*) * number_of_records); // * map->header->bytes_per_record
    if (!map->records) return -3; //malloc error

    for(int i = 0; i < number_of_records; i++ ) {
        //map->records[i] = full.records[ i + record_id_start ];
        size_t sz = map->header->bytes_per_record;
        map->records[i] = malloc(sz);
        memcpy(map->records[i], full.records[record_id_start + i], sz);
        printf("%zu bytes copied",sz);
    }
    
    map->num_records = number_of_records;
    map->header->start_epoch = map->records[0]->start_epoch;
    map->header->end_epoch = map->records[ map->num_records - 1 ]->end_epoch;

    printf("map->num_records : %zu \nmap->header->start_epoch : %f | map->records[0]->start_epoch : %f \nmap->header->end_epoch : %f | map->records[map->num_records - 1]->end_epoch : %f\n", map->num_records, map->header->start_epoch, map->records[0]->start_epoch, map->header->end_epoch, map->records[map->num_records - 1]->end_epoch);
    printf("size of entire allocatedd memory for the map : %zu \n", sizeof(MappedEphemeris) + sizeof(EphemerisFile_Header) + sizeof(EphemerisFile_Record*) * number_of_records + number_of_records * map->header->bytes_per_record);
    printf("size of map only : %zu\n", sizeof(MappedEphemeris));
    printf("size of header : %zu\n", sizeof(EphemerisFile_Header));
    printf("size of records pointers array : %zu\n", sizeof(EphemerisFile_Record*) * number_of_records);
    printf("size of all records stored : %d\n", number_of_records * map->header->bytes_per_record);
    ephemeris_unmap_file(&full);


    return 0;
}
*/

int spody_setup_partialMappedEphemerisData(MappedEphemerisData *med, const char *filename, double in_start_et, double in_end_et){

    /* Loads the same .spody file as spody_setup_MappedEphemerisData but
     * keeps in memory only the records covering [in_start_et, in_end_et]
     * (in ET seconds past J2000). Useful for missions of bounded duration
     * on memory-constrained targets. */

    MappedEphemerisData full = {0};
    if (ephemeris_map_file(&full, filename) != 0) return -1;

    int record_id_start = (int)floor((in_start_et - full.header->start_epoch) / (double)full.header->seconds_per_record);
    int record_id_end   = (int)floor((in_end_et   - full.header->start_epoch) / (double)full.header->seconds_per_record);
    int n = (record_id_end - record_id_start) + 1;
    if (n <= 0) { ephemeris_unmap_file(&full); return -4; }

    /* private copy of header + records (the full mmap will be unmapped) */
    med->header = malloc(sizeof(EphemerisFile_Header));
    if (!med->header) { ephemeris_unmap_file(&full); return -3; }
    *med->header = *full.header;

    med->records = (EphemerisFile_Record**)malloc(sizeof(EphemerisFile_Record*) * n);
    if (!med->records) { ephemeris_unmap_file(&full); return -3; }

    for (int i = 0; i < n; i++) {
        size_t sz = med->header->bytes_per_record;
        med->records[i] = malloc(sz);
        if (!med->records[i]) { ephemeris_unmap_file(&full); return -3; }
        memcpy(med->records[i], full.records[record_id_start + i], sz);
    }
    med->num_records = (size_t)n;

    /* refresh header epochs to the actual subset range (still in ET) */
    med->header->start_epoch = med->records[0]->start_epoch;
    med->header->end_epoch   = med->records[n - 1]->end_epoch;

    ephemeris_unmap_file(&full);
    return 0;
}