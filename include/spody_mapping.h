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
