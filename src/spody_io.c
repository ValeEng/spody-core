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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "spody_io.h"

/* ============================================================
 * Lifecycle
 * ============================================================ */
/* Copy a C string of arbitrary length into an 8-byte magic field:
 * zero-padded if shorter, truncated if longer. Result is NOT
 * zero-terminated (the field can use all 8 bytes for content). */
static void copy_magic(char dst[SPODY_IO_MAGIC_LEN], const char *src) {
    memset(dst, 0, SPODY_IO_MAGIC_LEN);
    if (!src) return;
    size_t n = strlen(src);
    if (n > SPODY_IO_MAGIC_LEN) n = SPODY_IO_MAGIC_LEN;
    memcpy(dst, src, n);
}

int spody_setup_LogBuffer(SpodyLogBuffer *buf,
                          const char *filename,
                          size_t capacity_records,
                          size_t record_size,
                          const char *magic) {
    if (!buf || !filename || capacity_records == 0 || record_size == 0) return -1;

    buf->data             = NULL;
    buf->capacity_records = 0;
    buf->record_size      = 0;
    buf->count            = 0;
    buf->fp               = NULL;
    memset(buf->magic, 0, SPODY_IO_MAGIC_LEN);

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("spody_setup_LogBuffer fopen");
        return -2;
    }

    /* write the on-disk header */
    SpodyLogHeader h;
    memset(&h, 0, sizeof(h));
    copy_magic(h.magic, magic);
    h.version       = SPODY_IO_HEADER_VERSION;
    h.header_size   = (uint32_t)sizeof(SpodyLogHeader);
    h.record_size   = (uint32_t)record_size;
    h.endian_marker = SPODY_IO_ENDIAN_MARKER;
    if (fwrite(&h, sizeof(h), 1, fp) != 1) {
        fclose(fp);
        return -3;
    }

    /* allocate ring */
    void *data = malloc(capacity_records * record_size);
    if (!data) {
        fclose(fp);
        return -4;
    }

    buf->data             = data;
    buf->capacity_records = capacity_records;
    buf->record_size      = record_size;
    buf->count            = 0;
    buf->fp               = fp;
    copy_magic(buf->magic, magic);

    return 0;
}

int spody_log_flush(SpodyLogBuffer *buf) {
    if (!buf) return -1;
    if (!buf->fp || !buf->data || buf->count == 0) return 0;
    size_t written = fwrite(buf->data, buf->record_size, buf->count, buf->fp);
    if (written != buf->count) return -2;
    buf->count = 0;
    return 0;
}

int spody_log_append(SpodyLogBuffer *buf, const void *record) {
    if (!buf || !record) return -1;
    if (!buf->fp || !buf->data) return -2;

    /* copy into the next slot */
    char *slot = (char*)buf->data + buf->count * buf->record_size;
    memcpy(slot, record, buf->record_size);
    buf->count++;

    /* flush if full */
    if (buf->count >= buf->capacity_records) {
        return spody_log_flush(buf);
    }
    return 0;
}

int spody_free_LogBuffer(SpodyLogBuffer *buf) {
    if (!buf) return -1;
    int rc = 0;
    if (buf->fp) {
        rc = spody_log_flush(buf);
        fclose(buf->fp);
        buf->fp = NULL;
    }
    if (buf->data) {
        free(buf->data);
        buf->data = NULL;
    }
    buf->capacity_records = 0;
    buf->record_size      = 0;
    buf->count            = 0;
    return rc;
}

/* ============================================================
 * Reading
 * ============================================================ */
int spody_log_read_header(FILE *fp, SpodyLogHeader *h_out) {
    if (!fp || !h_out) return -1;
    if (fseek(fp, 0, SEEK_SET) != 0) return -2;
    if (fread(h_out, sizeof(*h_out), 1, fp) != 1) return -3;
    if (h_out->endian_marker != SPODY_IO_ENDIAN_MARKER) return -4;
    if (h_out->version != SPODY_IO_HEADER_VERSION) return -5;
    if (h_out->record_size == 0) return -6;
    if (h_out->header_size < sizeof(*h_out)) return -7;

    /* leave fp positioned at the first record. If the file was written by a
     * future version with a larger header, header_size > sizeof(*h_out) and
     * we skip the trailing extra bytes; if it was written with this exact
     * version they are equal and the seek is a no-op. */
    if (fseek(fp, (long)h_out->header_size, SEEK_SET) != 0) return -8;
    return 0;
}

/* ============================================================
 * Binary -> CSV converter
 *
 * Walks the file record by record, casts each as double*, prints one
 * CSV row. The library doesn't know the field names: the caller passes
 * them as `column_names`.
 * ============================================================ */
int spody_log_dump_csv(const char *bin_path,
                       FILE *csv_out,
                       const char *const *column_names,
                       int n_columns,
                       int decimate) {
    if (!bin_path || !csv_out || !column_names || n_columns <= 0) return -1;
    if (decimate < 1) decimate = 1;

    FILE *bin = fopen(bin_path, "rb");
    if (!bin) {
        perror("spody_log_dump_csv fopen");
        return -2;
    }

    SpodyLogHeader h;
    int rc = spody_log_read_header(bin, &h);
    if (rc != 0) {
        fprintf(stderr, "spody_log_dump_csv: bad header in %s (rc=%d)\n", bin_path, rc);
        fclose(bin);
        return -3;
    }
    if (h.record_size != (uint32_t)(n_columns * (int)sizeof(double))) {
        fprintf(stderr, "spody_log_dump_csv: record_size=%u but n_columns*sizeof(double)=%zu\n",
                h.record_size, (size_t)(n_columns * sizeof(double)));
        fclose(bin);
        return -4;
    }

    /* CSV header */
    for (int i = 0; i < n_columns; i++) {
        fputs(column_names[i], csv_out);
        if (i + 1 < n_columns) fputc(',', csv_out);
    }
    fputc('\n', csv_out);

    /* read records into a stack buffer (record_size <= 4 KB in any sensible setup) */
    double row[1024];
    if ((size_t)n_columns > sizeof(row)/sizeof(row[0])) {
        fprintf(stderr, "spody_log_dump_csv: too many columns (%d)\n", n_columns);
        fclose(bin);
        return -5;
    }

    long n_in = 0, n_out = 0;
    while (fread(row, h.record_size, 1, bin) == 1) {
        if ((n_in % decimate) == 0) {
            for (int i = 0; i < n_columns; i++) {
                fprintf(csv_out, "%.18e", row[i]);
                if (i + 1 < n_columns) fputc(',', csv_out);
            }
            fputc('\n', csv_out);
            n_out++;
        }
        n_in++;
    }

    fclose(bin);
    fprintf(stderr, "spody_log_dump_csv: %ld records read, %ld rows written (decimate=%d)\n",
            n_in, n_out, decimate);
    return 0;
}
