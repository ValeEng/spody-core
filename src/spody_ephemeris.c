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

int read_ephemeris_file_header(FILE *file, EphemerisFile_Header *ep ) {
    //parsing for de440

    int found_group = 0; 
    int data_row = 0;    
    char line[BUFFER_SIZE_EPH];

    while (fgets(line, sizeof(line), file)) {
        
        // looking for GROUP 1050 
        if (!found_group) {
            if (strstr(line, "GROUP") || strstr(line, "KSIZE")){
                if(strstr(line, "KSIZE")) {
                    printf("RECORD SIZE found (normally the first row)\n");  
                    char *token = strtok(line, " \t\r\n");
                    while (token) {
                        if (strstr(token, "KSIZE")) {
                            token = strtok(NULL, " \t\r\n");
                            ep->bytes_per_record = strtod(token, NULL); 
                            printf("load KSIZE : %d\n",ep->bytes_per_record);
                        } else if (strstr(token, "NCOEFF")) {
                            token = strtok(NULL, " \t\r\n");
                            ep->number_coefficients_per_record = strtod(token, NULL); ; 
                            printf("load NCOEFF : %d\n",ep->number_coefficients_per_record);
                        }
                        token = strtok(NULL, " \t\r\n");
                    }
                    printf("RECORD SIZE %d | N COEFF %d\n", ep->bytes_per_record, ep->number_coefficients_per_record);
                    continue;              
                }else if (strstr(line, "1030")) {
                    found_group = 1030;
                    printf("GROUP 1030 found\n"); 
                }else if (strstr(line, "1050")) {
                    found_group = 1050;
                    printf("GROUP 1050 found\n"); 
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
                    ep->days_per_record = val; 
                }
                
                printf("Ele %d: %f\n", ele, val);
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
                
                printf("row %d, Ele %d: %f\n", data_row, ele, val);

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

double fds2cd(char *str) {
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

int read_record_block(FILE *fp, EphemerisFile_Header *ep, EphemerisFile_Record *eprec) {
   
    char line[BUFFER_SIZE_EPH];
    int n_coeff_expected = 0;
    printf("in read record\n");

    if (!fgets(line, sizeof(line), fp)) return 0; // it's end of file 
    //printf("Reading line: %s\n", line);
    if (sscanf(line, "%d %d", &eprec->record_number, &eprec->number_coefficients_per_record) != 2) return 0; // Parsing of number of coefficients
    //printf("Reading record %d with %d coefficients\n", eprec->record_number, eprec->number_coefficients_per_record);
    if (eprec->number_coefficients_per_record != ep->number_coefficients_per_record) {
        printf("Warning: Expected %d coefficients, but found %d in block %d\n", n_coeff_expected, eprec->number_coefficients_per_record, eprec->record_number);
        printf("HEADER info diverge from data blocks read.\n");
        printf("Check the ephemeris file integrity.\n");
        return 0;
    }else{
        n_coeff_expected = eprec->number_coefficients_per_record;
        printf("n_coeff_expected : %d\n",n_coeff_expected);
    }

    printf("cleaning ok\n");

    int coeff_idx = 0;
    char *token;
    while (coeff_idx < n_coeff_expected) { 
        if (fgets(line, sizeof(line), fp)){
            printf("Reading coefficients line: %s\n", line);
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
    printf("token loaded : %d\n",coeff_idx);
    eprec->start_epoch = eprec->record[0];
    eprec->end_epoch = eprec->record[1];

    return 1; //good 
}

int create_binary_ephemeris_file(EphemerisFile_Header *ep, int *old_epoch, const char *ascp_filename, const char *bin_filename) {
    printf("in create binary\n");
    FILE *fp_ascp = fopen(ascp_filename, "r");
    if (!fp_ascp) { perror("Errore file ascp"); return -1; }

    FILE *fp_bin = fopen(bin_filename, "ab");  // "ab" for append mode
    if (!fp_bin) { perror("Errore file bin"); fclose(fp_ascp); return -1; }
    printf("file loaded\n");

    double t_start, t_end;
    size_t record_size = ep->bytes_per_record;
    printf("sizeof(EphemerisFile_Record) : %zu\n", sizeof(EphemerisFile_Record));
    printf("ep->number_coefficients_per_record * sizeof(double) : %zu\n",ep->number_coefficients_per_record * sizeof(double));
    printf("read -> bytes_per_record : %zu\n",record_size);

    EphemerisFile_Record *eprec = malloc(record_size);
    if (!eprec) {
        perror("Error memory allocation for EphemerisFile_Record");
        fclose(fp_ascp);
        fclose(fp_bin);
        return -1;
    }
    printf("malloc ok\n");
    
    int block_count = 0;
    // Read each record block and write to binary file
    while (read_record_block(fp_ascp, ep, eprec)) {
        block_count++;
        printf("block %d readed\n",block_count);
        printf("old epoch : %d\n",*old_epoch);



        if (*old_epoch == (int) eprec->start_epoch ){ //eprec->record_number == 1 && ep->start_epoch != eprec->start_epoch && strcmp(bin_filename, "de440.spody") == 0
            printf("\n\n! a clone record foud ! record %d start date %.2f \n\n",eprec->record_number, eprec->start_epoch);
            continue;
        } 
        *old_epoch = eprec->start_epoch;
        printf("new old epoch : %d\n",*old_epoch);


        printf("Writing record %d to binary file\n", eprec->record_number);
        printf("size of : %zu\n",record_size);
        fwrite(eprec, record_size, 1, fp_bin);
        printf("Written record %04d with %04d coefficients from %.6f to %.6f\n", eprec->record_number, eprec->number_coefficients_per_record, eprec->start_epoch, eprec->end_epoch);
    
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

double chebyshev_evaluate(double time_scaled, const double *coefficients, int n_coeffs) {
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
 * @param jd_epoch Target Julian Date.
 * @param ef Ephemeris JPL file structure (from ASCII Header).
 * @param coeffs Array of coefficients from the ASCP file.
 * @param t_start Julian Date at the start of the current block.
 * @param t_end Julian Date at the end of the current block.
 * @param result Output array for the position [X, Y, Z].
 * @return 1 if successful, 0 if the date is out of range or data is invalid.
 */

int calculate_body_position(MappedEphemeris *map, int target_idx, double jd_epoch, double result[3]) {
    
    #if DEBUG_EPHEMERIS == 1
    printf("\n");
    for(int i=0;i<15;i++){
        printf("bodies mapped idx %02d | location %04d | n coeff per comp %02d | n sets per record %d\n", i+1, map->header->location[i], map->header->number_coefficients_per_component[i], map->header->number_complete_sets_coefficients_per_record[i]);
    }
    #endif

    // 1. Retrieve body parameters
    int numeber_coefficients_per_component = map->header->number_coefficients_per_component[target_idx];
    int start_index = map->header->location[target_idx] - 1; // 1-based -> 0-based
    int n_components = 3; // X, Y, Z <--- not to flexible!!!!

    // 2. Subdivision identification 
    int record_id = (int)floor((jd_epoch - map->header->start_epoch)/map->header->days_per_record);
    //printf("Record ID calc: %d\n", record_id);
    EphemerisFile_Record *eprec = map->records[record_id];
    double block_duration = eprec->end_epoch - eprec->start_epoch;
    double set_duration = block_duration / map->header->number_complete_sets_coefficients_per_record[target_idx];
    int set_id = (int)floor((jd_epoch - eprec->start_epoch) / set_duration); //wich subdivision we need 
    
    #if DEBUG_EPHEMERIS == 1
    printf("Record ID: %d | Start epoch header: %f | start epoch record [%d]: %f | JD epoch now : %f \n", record_id, map->header->start_epoch, record_id, map->records[record_id]->start_epoch, jd_epoch);
    printf("Set ID: %d\n", set_id);
    #endif

    // floating point error check 
    if (set_id >= map->header->number_complete_sets_coefficients_per_record[target_idx]) {
        set_id = map->header->number_complete_sets_coefficients_per_record[target_idx] - 1;
        printf("Adjusted Set ID due to floating point error: %d\n", set_id);
    }  // we can do this because jd_epoch <= t_end

    // 3. Tau evaluation 
    double t_gran_start = eprec->start_epoch + set_id * set_duration;
    double t_gran_end = t_gran_start + set_duration;
    
    // normalizing of tau --> [-1.0, 1.0] 
    double tau = (2.0 * jd_epoch - t_gran_start - t_gran_end) / (t_gran_end - t_gran_start); // prob not the best for floting point error 
    double midpoint = 0.5 * (t_gran_start + t_gran_end);
    double half_range = 0.5 * (t_gran_end - t_gran_start);
    //double tau2 = (jd_epoch - midpoint) / half_range; //TBD with some tests
    
    #if DEBUG_EPHEMERIS == 1
    printf("Subdivision number: %d|%d, Tau: LP %.21f HP %.21f\n", set_id,map->header->number_complete_sets_coefficients_per_record[target_idx], tau, tau2);
    printf("Gran Start Time: %f | Gran End Time: %f\n",t_gran_start , t_gran_end);
    #endif

    // 4. Buffer param    
    // a complete set (X, Y, Z) per set_id is 3 * n_coeffs.
    int set_length = n_components * numeber_coefficients_per_component; 
    int offset = start_index + (set_id * set_length);

    // 5. Polinomial evaluation for X, Y, Z
    for (int i = 0; i < n_components; i++) {
        // X ---> offset
        // Y ---> offset + n_coeffs
        // Z ---> offset + 2 * n_coeffs

        #if DEBUG_EPHEMERIS == 1
        printf("Calculating component %d with coeffs starting at index %d\n", i, offset + (i * numeber_coefficients_per_component));
        #endif

        const double *coeff_ptr = eprec->record + offset + (i * numeber_coefficients_per_component);
        result[i] = chebyshev_evaluate(tau, coeff_ptr, numeber_coefficients_per_component);
    }

    return 1;
}

//MAPPING FUNCTIONS--------------------------------------------------------------------------------------------------   

int ephemeris_map_file(MappedEphemeris *map, const char *filename) {
    
    if (!map) return -1; //if NULL exit

    if (mf_map_file(&map->mf, filename) != 0) return -2; //mapping error

    map->header = (EphemerisFile_Header*)map->mf.ptr;
    //printf("Mapped header start epoch: %.6f\n", map->header->start_epoch);
    //printf("Mapped header end epoch: %.6f\n", map->header->end_epoch);
    //printf("bytes per record: %d\n", map->header->bytes_per_record);

    printf("mf_size : %zu\n",map->mf.size);
    printf("sizeof(EphemerisFile_Header) : %zu\n",sizeof(EphemerisFile_Header));

    size_t remaining_bytes = map->mf.size - sizeof(EphemerisFile_Header);
    printf("remaning bytes : %zu \n",remaining_bytes);
    


    map->num_records = remaining_bytes / map->header->bytes_per_record; //we hope it is exact division
    
    #if DEBUG_EPHEMERIS == 1
    printf("bytes per record: %d\n", map->header->bytes_per_record);
    printf("Number of records mapped: %zu\n", map->num_records);
    #endif

    map->records = (EphemerisFile_Record**)malloc(sizeof(EphemerisFile_Record*) * map->num_records); // * map->header->bytes_per_record
    if (!map->records) return -3; //malloc error

    uint8_t *ptr = (uint8_t*)map->mf.ptr + sizeof(EphemerisFile_Header); //work on bytes 
    for (size_t i = 0; i < map->num_records; i++) {
        //printf("Mapping record %zu at address %p\n", i, ptr);
        map->records[i] = (EphemerisFile_Record*)ptr;
        #if DEBUG_EPHEMERIS == 1
        printf("Mapped record %zu: record number %d, n_coeff %d, start_epoch %.6f, end_epoch %.6f\n", i+1, map->records[i]->record_number, map->records[i]->number_coefficients_per_record, map->records[i]->start_epoch, map->records[i]->end_epoch);
        #endif
        size_t record_bytes = map->header->bytes_per_record; //(sizeof(EphemerisFile_Record) + map->records[i]->number_coefficients_per_record * sizeof(double));
        //printf("record_bytes : %zu\n",record_bytes);
        ptr += record_bytes;
    }

    return 0;
}

int ephemeris_unmap_file(MappedEphemeris *map) {
    if (!map) return -1;

    free(map->records);
    map->records = NULL;
    map->header = NULL;
    map->num_records = 0;

    return mf_unmap_file(&map->mf);
}

int spody_createfile_MappedEphemeris(const char *path, const char **file_names, const int n_files, const char *de){
    
    char header_path[1000]; //TBD 
    char bin_filename[1000];
    char ascp_filename[1000];

    int returnNumber;

    sprintf(header_path, "./%s/header.%s",path,de); 
    sprintf(bin_filename, "./%s/de%s.spody",path,de); 

    FILE *file = fopen(header_path,"r");
    if (!file) { perror("Errore file header"); return -1; }
    
    FILE *fp_bin = fopen(bin_filename, "ab");  // "ab" for append mode
    if (!fp_bin) { perror("Errore file bin"); fclose(file); return -1; }

    EphemerisFile_Header ep = {0};
    returnNumber = read_ephemeris_file_header(file, &ep);
    fclose(file);

    size_t n = ep.number_coefficients_per_record;
    size_t record_size = sizeof(EphemerisFile_Record) + n * sizeof(double);
    ep.bytes_per_record = (int)record_size; //TBD update record size in header necessary understand KSIZE param
    printf("ep.bytes_per_record : %d\n",ep.bytes_per_record);

        //first write header info
    fwrite(&ep, sizeof(EphemerisFile_Header), 1, fp_bin);
    fclose(fp_bin);
    
    printf("header writed, size : %zu\n",sizeof(EphemerisFile_Header));

    int old_epoch = 1 ; //necessary to avoid duplicate 

    for (int i = 0; i < n_files; i++){
        
        sprintf(ascp_filename, "./%s/ascp%s.%s",path,file_names[i],de); 
        //FILE *file = fopen(ascp_filename,"r");
        //if (!file) { perror("Errore file ASCP"); return -1; }
        
        returnNumber = create_binary_ephemeris_file(&ep,&old_epoch,ascp_filename,bin_filename); //TBD the append mode in create_binary_ephemeris_file
        printf("old epoch : %d\n",old_epoch);

        fclose(file);
        printf("\n\nfile %d writed\n",i);

    }

    return returnNumber;

}

int spody_setup_MappedEphemeris(MappedEphemeris *map, const char *filename){

    int returnNumber =  ephemeris_map_file(map, filename);

    if (map->header->start_epoch != map->records[0]->start_epoch){
        printf("!It is a subset!\n");
        printf("Header and ephemeris file have different start epoch.\n");
        printf("Change header start epoch to the subset JD start.\n");
        map->header->start_epoch = map->records[0]->start_epoch;
    }
    if (map->header->end_epoch != map->records[map->num_records-1]->end_epoch){
        printf("!It is a subset!\n");
        printf("Header and ephemeris file have different end epoch.\n");
        printf("Change header end epoch to the subset JD end.\n");
        map->header->end_epoch = map->records[map->num_records-1]->end_epoch;
    }

    return returnNumber;
}

int spody_get_position(MappedEphemeris *map, int central_idx ,int target_idx, double jd_epoch, double result[3]){

/*****************NAIF ID********************
    SUN     = 10
    EARTH   = 399
    MOON    = 301
********************************************/

    //central_idx and target_idx are naif code 
    int returnNumber;
    double temp[3]={0.0};
    
    /*
    int central_body, target_body;

    if (strstr(central_idx, "Earth")||strstr(central_idx, "earth")||strstr(central_idx, "EARTH") ){ 
        central_body = 399; //Earth Barycenter
    }else if (strstr(central_idx, "Moon")||strstr(central_idx, "moon")||strstr(central_idx, "MOON") ){ 
        central_body = 301; //Moon
    }
    
    if (strstr(target_body, "Earth")||strstr(target_body, "earth")||strstr(target_body, "EARTH") ){ 
        target_body = 399; //Earth Barycenter
    }else if (strstr(central_idx, "Moon")||strstr(target_body, "moon")||strstr(target_body, "MOON") ){ 
        target_body = 301; //Moon
    }else if (strstr(central_idx, "Sun")||strstr(target_body, "sun")||strstr(target_body, "SUN") ){ 
        target_body = 10; //Moon
    }
    */

    switch (central_idx)
    {
    case 399:
        
        switch (target_idx)
        {
        case 301:
            returnNumber = calculate_body_position(map, 9, jd_epoch, result); //TBD parse the arrey of bodies to find the correct id of the map
            return 0;
        case 10:
            returnNumber = calculate_body_position(map, 9, jd_epoch, result); //earth-moon distance
            //earth distance from earth-moon barycenter
            result[0] = -(1/(1+EMRAT)) * result[0];
            result[1] = -(1/(1+EMRAT)) * result[1];
            result[2] = -(1/(1+EMRAT)) * result[2];

            returnNumber = calculate_body_position(map, 2, jd_epoch, temp); //erth-moon barycenter to solar system barycenter
            //earth distance from solar system barycenter 
            result[0] = result[0] + temp[0];
            result[1] = result[1] + temp[1];
            result[2] = result[2] + temp[2];

            returnNumber = calculate_body_position(map, 10, jd_epoch, temp); //sun from solar system barycenter
            result[0] = temp[0] - result[0];
            result[1] = temp[1] - result[1];
            result[2] = temp[2] - result[2];
            return 0;        
        default:
            printf("Target body not supported");
            result[0] = 0.0;
            result[1] = 0.0;
            result[2] = 0.0;
            return -1; // out from the switch 
        } 

        break;
    case 301:

        switch (target_idx)
        {
        case 399:
            returnNumber = calculate_body_position(map, 9, jd_epoch, result); //TBD parse the arrey of bodies to find the correct id of the map
            result[0] = -result[0];
            result[1] = -result[1];
            result[2] = -result[2];
            return 0;
        case 10: //TBD
            returnNumber = calculate_body_position(map, 9, jd_epoch, result); //earth-moon distance
            //earth distance from earth-moon barycenter
            result[0] = -MOON_MU / (EARTH_MU + MOON_MU) * result[0];
            result[1] = -MOON_MU / (EARTH_MU + MOON_MU) * result[1];
            result[2] = -MOON_MU / (EARTH_MU + MOON_MU) * result[2];

            returnNumber = calculate_body_position(map, 2, jd_epoch, temp); //erth-moon barycenter to solar system barycenter
            result[0] = result[0] + temp[0];
            result[1] = result[1] + temp[1];
            result[2] = result[2] + temp[2];

            returnNumber = calculate_body_position(map, 10, jd_epoch, temp); //sun from solar system barycenter
            result[0] = temp[0] - result[0];
            result[1] = temp[1] - result[1];
            result[2] = temp[2] - result[2];
            return 0;        
        default:
            result[0] = 0.0;
            result[1] = 0.0;
            result[2] = 0.0;
            return -1; // out from the switch 
        }

        break;
    default:
        printf("Target body not supported");
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
        return -1; // out from the switch 
    }

    printf("Central body not supported");
    return -2; // function starting error
}

int spody_get_lunarlibrationangles(MappedEphemeris *map, double jd_epoch, double result[3]){
    
    calculate_body_position(map, 12, jd_epoch, result);

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

int spody_setup_partialMappedEphemeris(MappedEphemeris *map, const char *filename, double in_start, double in_end){

    //TBD necessary a free function but now we end the program for memory free

    MappedEphemeris full = {0};
    if (ephemeris_map_file(&full, filename) != 0) return -1;

    int record_id_start = (int)floor((in_start - full.header->start_epoch)/full.header->days_per_record);
    int record_id_end = (int)floor((in_end - full.header->start_epoch)/full.header->days_per_record);
    int number_of_records = (record_id_end - record_id_start) + 1; // is alway + 1 wrt the difference 
    if (number_of_records <= 0) return -4;

    printf("start : %.6f | end : %.6f\n", in_start, in_end);
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
        printf("record start epoch : %.6f",map->records[i]->start_epoch);
        map->records[i]->start_epoch -= in_start; 
        map->records[i]->start_epoch *= SECONDSxDAY; 
        printf("record start epoch : %.6f ",map->records[i]->end_epoch);
        map->records[i]->end_epoch -= in_start; 
        map->records[i]->end_epoch *= SECONDSxDAY; 

        printf("%zu bytes copied\n",sz);
    }
    
    map->num_records = number_of_records;
    map->header->days_per_record *= SECONDSxDAY ; //now we are working with seconds 
    map->header->start_epoch = map->records[0]->start_epoch;
    map->header->end_epoch = map->records[ map->num_records - 1 ]->end_epoch;

    printf("map->header->days_per_record (seconds_per_record) : %d\n",map->header->days_per_record);
    printf("map->num_records : %zu \nmap->header->start_epoch : %f | map->records[0]->start_epoch : %f \nmap->header->end_epoch : %f | map->records[map->num_records - 1]->end_epoch : %f\n", map->num_records, map->header->start_epoch, map->records[0]->start_epoch, map->header->end_epoch, map->records[map->num_records - 1]->end_epoch);
    printf("size of entire allocatedd memory for the map : %zu \n", sizeof(MappedEphemeris) + sizeof(EphemerisFile_Header) + sizeof(EphemerisFile_Record*) * number_of_records + number_of_records * map->header->bytes_per_record);
    printf("size of map only : %zu\n", sizeof(MappedEphemeris));
    printf("size of header : %zu\n", sizeof(EphemerisFile_Header));
    printf("size of records pointers array : %zu\n", sizeof(EphemerisFile_Record*) * number_of_records);
    printf("size of all records stored : %d\n", number_of_records * map->header->bytes_per_record);
    ephemeris_unmap_file(&full);


    return 0;
}