/*
 * Runtime accessors for the version constants generated into
 * spody_version.h by CMake at configure time.
 *
 * The compile-time macros (e.g. SPODY_VERSION_STRING) are available to any
 * translation unit that includes <spody_version.h>. These three functions
 * return the SAME strings, but coming from the compiled library -- which
 * is the only reliable way to detect a header/library version skew at
 * load time.
 */
#include "spody_version.h"

const char *spody_version(void)         { return SPODY_VERSION_STRING; }
const char *spody_git_hash(void)        { return SPODY_GIT_HASH; }
const char *spody_build_timestamp(void) { return SPODY_BUILD_TIMESTAMP; }
