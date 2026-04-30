/*
 * One-shot helper: generate raw_data/DE440/de440.spody from the JPL ASCII files.
 *
 * Run from the repository root (so relative paths resolve correctly):
 *   build/tvb/Release/gen_de440_spody.exe
 *
 * The generated binary will be written to ./raw_data/DE440/de440.spody.
 */
#include <stdio.h>
#include "spody_core.h"

int main(void) {

    const char *path = "raw_data/DE440";
    const char *de   = "440";

    // ASCII chunks present under raw_data/DE440/ascpXXXXX.440
    const char *dates[] = {
        "01550", "01650", "01750", "01850", "01950",
        "02050", "02150", "02250", "02350", "02450", "02550"
    };
    const int n_files = (int)(sizeof(dates) / sizeof(dates[0]));

    printf("Generating %s/de%s.spody from %d ASCII chunk(s)...\n", path, de, n_files);

    int rc = spody_createfile_MappedEphemerisData(path, dates, n_files, de);
    if (rc != 0) {
        printf("FAIL: spody_createfile_MappedEphemerisData returned %d\n", rc);
        return 1;
    }

    printf("DONE: ./%s/de%s.spody generated.\n", path, de);
    return 0;
}
