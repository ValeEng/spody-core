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
#ifndef SPODY_IO_H
#define SPODY_IO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

/* ============================================================
 * Generic record-based binary log
 *
 * spody_io is intentionally domain-agnostic. It buffers fixed-size
 * binary records to disk with a small header that identifies the
 * file kind (magic + version + record_size + header_size). Higher-
 * level modules (Mission, custom user code) define their own record
 * struct types and call append() per record. The library doesn't
 * care what's in the records as long as their byte layout is fixed.
 *
 * Why a 1024-record buffer (default): a typical trajectory record
 * is ~56 bytes; 1024 records = ~56 KB per flush, well below any
 * sensible OS write granularity. Ten thousand simulation steps then
 * cost ~10 syscall round-trips instead of 10000. Configurable via
 * the capacity_records argument to spody_setup_LogBuffer.
 *
 * Endianness: the on-disk format is little-endian (the only platform
 * that matters in 2026). A marker is written in the header to detect
 * the mistake on read. No byte-swapping is performed.
 * ============================================================ */

#define SPODY_IO_HEADER_VERSION  1u
#define SPODY_IO_ENDIAN_MARKER   0x01020304u   /* read as 0x04030201 -> wrong endian */
#define SPODY_IO_MAGIC_LEN       8

/* Default capacity of the in-memory ring (records). Used by spody_io
 * itself and re-exported to higher-level modules (e.g. spody_mission). */
#define SPODY_LOG_BUFFER_CAPACITY 1024

/* ============================================================
 * Standard magic strings (8 bytes each, no terminator).
 *
 * Higher-level modules and user code should pick a magic from this
 * list (or define their own using the same 8-byte convention) and
 * pass it to spody_setup_LogBuffer. The constants below guarantee
 * a correctly-sized null-padded literal.
 * ============================================================ */
#define SPODY_MAGIC_TRAJECTORY  "SPDYTRAJ"     /* state vector log         */
#define SPODY_MAGIC_BREAKDOWN   "SPDYBRKD"     /* per-force breakdown log  */
#define SPODY_MAGIC_EVENTS      "SPDYEVNT"     /* triggered events log     */
#define SPODY_MAGIC_USER1       "SPDYUSR1"     /* user-defined log slot 1  */
#define SPODY_MAGIC_USER2       "SPDYUSR2"     /* user-defined log slot 2  */

/* ============================================================
 * On-disk header (24 bytes, written once at file open).
 *
 * Layout is fixed for v1. Future versions may extend the trailing
 * part: readers can skip exactly `header_size` bytes from the start
 * of the file to find the first record, even if they don't recognise
 * extra fields.
 *
 *   offset  size  field
 *     0      8    magic
 *     8      4    version           SPODY_IO_HEADER_VERSION
 *    12      4    header_size       sizeof(SpodyLogHeader) at write time
 *    16      4    record_size       bytes per record
 *    20      4    endian_marker     SPODY_IO_ENDIAN_MARKER
 * ============================================================ */
typedef struct {
    char     magic[SPODY_IO_MAGIC_LEN];   /* e.g. "SPDYTRAJ", "SPDYBRKD" */
    uint32_t version;                     /* SPODY_IO_HEADER_VERSION     */
    uint32_t header_size;                 /* bytes to skip before first record */
    uint32_t record_size;                 /* bytes per record            */
    uint32_t endian_marker;               /* SPODY_IO_ENDIAN_MARKER      */
} SpodyLogHeader;

/* ============================================================
 * Buffer + sink (record-agnostic)
 *
 * Threading: not safe to share a SpodyLogBuffer across threads.
 * Each thread that wants its own log must own its own buffer
 * (and its own FILE*).
 * ============================================================ */
typedef struct {
    void   *data;                /* raw buffer = capacity_records * record_size */
    size_t  capacity_records;
    size_t  record_size;
    size_t  count;               /* records currently buffered (not yet flushed) */
    FILE   *fp;                  /* output file (owned: closed on free)          */
    char    magic[SPODY_IO_MAGIC_LEN];
} SpodyLogBuffer;

/* ============================================================
 * Lifecycle
 * ============================================================ */

/* Open `filename` for writing in binary mode, write the SpodyLogHeader,
 * allocate an internal buffer of `capacity_records` records.
 * Returns 0 on success.
 *
 * `magic` is a C string of any length: it is copied into an 8-byte
 * field, zero-padded if shorter, truncated if longer. To avoid mistakes,
 * pass one of the SPODY_MAGIC_* constants from this header. */
int spody_setup_LogBuffer(SpodyLogBuffer *buf,
                          const char *filename,
                          size_t capacity_records,
                          size_t record_size,
                          const char *magic);

/* Append one record. Triggers an automatic flush if the buffer fills up.
 * Returns 0 on success, non-zero on I/O error. */
int spody_log_append(SpodyLogBuffer *buf, const void *record);

/* Force a flush of any buffered records to disk. Called automatically
 * by spody_free_LogBuffer; call it manually if you want disk durability
 * at a known point (e.g. periodic checkpoints). */
int spody_log_flush(SpodyLogBuffer *buf);

/* Flush, close the file, free the buffer. Safe to call on an unopened
 * (zeroed) buffer. */
int spody_free_LogBuffer(SpodyLogBuffer *buf);

/* ============================================================
 * Reading
 * ============================================================ */

/* Read and validate the SpodyLogHeader from `fp` (must be at offset 0).
 * Returns 0 on success, non-zero if the header is missing/corrupted/
 * has the wrong endian marker. After success, fp is positioned at the
 * start of the first record. */
int spody_log_read_header(FILE *fp, SpodyLogHeader *h_out);

/* ============================================================
 * Binary -> CSV converter
 *
 * Reads `bin_path`, validates the header, then for each record dumps
 * one CSV row containing all the doubles in the record (cast as double*).
 * `column_names` are written as the CSV header line. `n_columns` must
 * match record_size / sizeof(double); the function fails otherwise.
 *
 * `decimate` >= 1: write only one record every `decimate` records
 * (1 = all, 10 = every 10th, etc.).
 *
 * Returns 0 on success.
 * ============================================================ */
int spody_log_dump_csv(const char *bin_path,
                       FILE *csv_out,
                       const char *const *column_names,
                       int n_columns,
                       int decimate);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_IO_H */
