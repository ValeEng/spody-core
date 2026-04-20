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
#include "spody_mapping.h"

#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#endif

int mf_map_file(MappedFile *mf, const char *filename) {
    if (!mf || !filename) return -1; //if NULL exit

#ifdef _WIN32
    HANDLE hFile = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) return -2;

    LARGE_INTEGER size;
    if (!GetFileSizeEx(hFile, &size)) {
        CloseHandle(hFile);
        return -3;
    }

    HANDLE hMap = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!hMap) {
        CloseHandle(hFile);
        return -4;
    }

    void *ptr = MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, 0);
    if (!ptr) {
        CloseHandle(hMap);
        CloseHandle(hFile);
        return -5;
    }

    mf->hFile = hFile;
    mf->ptr = ptr;
    mf->size = (size_t)size.QuadPart;

#else
    int fd = open(filename, O_RDONLY);
    if (fd < 0) return -2;

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        return -3;
    }

    void *ptr = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (ptr == MAP_FAILED) {
        close(fd);
        return -4;
    }

    mf->hFile = fd;
    mf->ptr = ptr;
    mf->size = (size_t)st.st_size;
#endif

    return 0;
}

int mf_unmap_file(MappedFile *mf) {
    if (!mf || !mf->ptr) return -1;

#ifdef _WIN32
    UnmapViewOfFile(mf->ptr);
    CloseHandle(mf->hFile);
#else
    munmap(mf->ptr, mf->size);
    close(mf->hFile);
#endif

    mf->ptr = NULL;
    mf->size = 0;

    return 0;
}
