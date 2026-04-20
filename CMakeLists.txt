# =============================================================================
# Spody Core - CMake build configuration
# =============================================================================
# This file configures how Spody Core is built.
#
# QUICK START
# -----------
#   mkdir build && cd build
#   cmake ..
#   cmake --build . --config Release        # Visual Studio / multi-config
#   cmake --build .                         # Makefile / Ninja (single-config)
#
# CONFIGURABLE OPTIONS (pass with -D<NAME>=<VALUE> on the cmake command line)
# --------------------------------------------------------------------------
#   SPODY_FAST_MATH        ON/OFF (default OFF)
#       Enables aggressive floating-point optimizations (/fp:fast on MSVC,
#       -ffast-math on GCC/Clang). Faster but may introduce small numerical
#       drift over long integrations. Leave OFF for bit-reproducible results.
#
#   SPODY_WHOLE_PROGRAM_OPT ON/OFF (default ON)
#       Enables whole-program / link-time optimization in Release builds
#       (/GL + /LTCG on MSVC, -flto on GCC/Clang). Slower link, faster runtime.
#
#   SPODY_SILENCE_MSVC_CRT_WARNINGS ON/OFF (default ON)
#       On MSVC, silences warnings about "unsafe" standard C functions
#       (strtok, sprintf, strncpy, sscanf). These are portable ISO C functions
#       and the MSVC "_s" alternatives are non-portable Microsoft extensions.
# =============================================================================

cmake_minimum_required(VERSION 3.10)
project(spody_core VERSION 1.0.0 LANGUAGES C)

# -----------------------------------------------------------------------------
# User-configurable options
# -----------------------------------------------------------------------------
option(SPODY_FAST_MATH
    "Enable aggressive floating-point optimizations (may affect precision)" OFF)

option(SPODY_WHOLE_PROGRAM_OPT
    "Enable whole-program / link-time optimization in Release builds" ON)

option(SPODY_SILENCE_MSVC_CRT_WARNINGS
    "Silence MSVC warnings about 'unsafe' standard C functions" ON)

# -----------------------------------------------------------------------------
# C standard
# -----------------------------------------------------------------------------
# C99 is required. Scientific code relies on math.h features (fmax, fmin,
# isnan, ...) and compound literals that are not guaranteed in C89.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)  # Prefer strict standard over GNU extensions

# -----------------------------------------------------------------------------
# Default build type
# -----------------------------------------------------------------------------
# Single-config generators (Makefile, Ninja) pick one build type at configure
# time. If the user did not specify one, default to Release so that a plain
# `cmake ..` produces optimized output, not a slow debug build.
#
# Multi-config generators (Visual Studio, Xcode) ignore CMAKE_BUILD_TYPE and
# let the user pick per-build with `cmake --build . --config <Type>`.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Build type: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

# -----------------------------------------------------------------------------
# Source files
# -----------------------------------------------------------------------------
set(SOURCES
    src/spody_eclipse.c
    src/spody_ephemeris.c
    src/spody_harmonics.c
    src/spody_mapping.c
    src/spody_math.c
)

# -----------------------------------------------------------------------------
# Library target
# -----------------------------------------------------------------------------
# Built as STATIC by default. Produces:
#   - libspody_core.a      on Linux / macOS
#   - spody_core.lib       on Windows (MSVC)
add_library(spody_core STATIC ${SOURCES})

# -----------------------------------------------------------------------------
# Compiler flags
# -----------------------------------------------------------------------------
if(MSVC)
    # ----- MSVC (Visual Studio) -----

    # Silence warnings about "unsafe" ISO C functions.
    # These are portable and correct; the "_s" alternatives are Microsoft-only.
    if(SPODY_SILENCE_MSVC_CRT_WARNINGS)
        target_compile_definitions(spody_core PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif()

    # Release-only optimization flags.
    # MSVC already applies /O2 by default in Release; we list it explicitly
    # for clarity and add a few extras useful for numerical code.
    target_compile_options(spody_core PRIVATE
        $<$<CONFIG:Release>:/O2>   # Maximum optimization for speed
        $<$<CONFIG:Release>:/Oi>   # Enable intrinsic functions (memcpy, sqrt...)
    )

    if(SPODY_WHOLE_PROGRAM_OPT)
        # /GL enables whole-program optimization (cross-TU inlining).
        # It requires /LTCG at link time to actually produce output.
        target_compile_options(spody_core PRIVATE
            $<$<CONFIG:Release>:/GL>
        )
        set_target_properties(spody_core PROPERTIES
            STATIC_LIBRARY_OPTIONS_RELEASE "/LTCG"
        )
    endif()

    if(SPODY_FAST_MATH)
        target_compile_options(spody_core PRIVATE
            $<$<CONFIG:Release>:/fp:fast>
        )
    endif()

else()
    # ----- GCC / Clang -----

    # Release-only optimization flags.
    target_compile_options(spody_core PRIVATE
        $<$<CONFIG:Release>:-O2>
        $<$<CONFIG:Release>:-DNDEBUG>
    )

    if(SPODY_WHOLE_PROGRAM_OPT)
        target_compile_options(spody_core PRIVATE
            $<$<CONFIG:Release>:-flto>
        )
        target_link_options(spody_core PRIVATE
            $<$<CONFIG:Release>:-flto>
        )
    endif()

    if(SPODY_FAST_MATH)
        target_compile_options(spody_core PRIVATE
            $<$<CONFIG:Release>:-ffast-math>
        )
    endif()

endif()

# -----------------------------------------------------------------------------
# Include directories
# -----------------------------------------------------------------------------
# PUBLIC: consumers of the library need these paths on their include path.
#
# BUILD_INTERFACE  -> used when the library is built in-tree (add_subdirectory
#                     or direct build of this project).
# INSTALL_INTERFACE -> used when the library has been installed via
#                     `cmake --install .` and consumed as a package.
target_include_directories(spody_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# -----------------------------------------------------------------------------
# System libraries
# -----------------------------------------------------------------------------
# On Unix/macOS, math functions (sqrt, sin, cos, ...) live in libm and must
# be linked explicitly. On Windows they are part of the CRT automatically.
if(UNIX OR APPLE)
    target_link_libraries(spody_core PRIVATE m)
endif()

# -----------------------------------------------------------------------------
# Installation rules
# -----------------------------------------------------------------------------
# `cmake --install . --prefix <path>` will produce:
#   <path>/lib/libspody_core.a     (or spody_core.lib on Windows)
#   <path>/include/spody_*.h       (all public headers)
install(TARGETS spody_core
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Install all public headers. Using a directory install means new headers
# added under include/ are picked up automatically, no edit required here.
install(DIRECTORY include/ DESTINATION include)

# -----------------------------------------------------------------------------
# Configuration summary
# -----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Spody Core configuration:")
message(STATUS "  Build type              : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fast math               : ${SPODY_FAST_MATH}")
message(STATUS "  Whole-program opt       : ${SPODY_WHOLE_PROGRAM_OPT}")
if(MSVC)
    message(STATUS "  Silence MSVC CRT warns  : ${SPODY_SILENCE_MSVC_CRT_WARNINGS}")
endif()
message(STATUS "")
