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
#ifndef  SPODY_MAPPING_H
#define  SPODY_MAPPING_H

#include <stddef.h>

#ifdef _WIN32
#include <windows.h>
#define MAPPING_HANDLE HANDLE
#else
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>
#define MAPPING_HANDLE int
#endif

typedef struct {
    MAPPING_HANDLE hFile;
    void *ptr;
    size_t size;
} MappedFile;

#ifdef __cplusplus
extern "C" {
#endif

int mf_map_file(MappedFile *mf, const char *filename);
int mf_unmap_file(MappedFile *mf);

#ifdef __cplusplus
}
#endif

#endif // SPODY_MAPPING_H
