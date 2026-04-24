# =============================================================================
# Tests, Validations, and Benchmarks (TVB) for spody_core
# =============================================================================
# Built only when the root configure sets -DSPODY_BUILD_TVB=ON.
#
# Each executable is an isolated main() that links against spody_core.
# No external framework is required.
#
# Some tests need a path to a compiled DE440 binary (*.spody). Pass it via
# -DSPODY_DE440_PATH=<path>; without it the ctest entries that need data
# are skipped.
# =============================================================================

enable_testing()

# ---- Benchmarks ----
add_executable(bench_ephemeris_batch benchmarks/bench_ephemeris_batch.c)
target_link_libraries(bench_ephemeris_batch PRIVATE spody_core)

# ---- Tests ----
add_executable(test_ephemeris_moon_sun tests/test_ephemeris_moon_sun.c)
target_link_libraries(test_ephemeris_moon_sun PRIVATE spody_core)

# Register ctest entries only if a data path was provided.
if(SPODY_DE440_PATH)
    add_test(NAME ephemeris_moon_sun
             COMMAND test_ephemeris_moon_sun "${SPODY_DE440_PATH}")
else()
    message(STATUS "SPODY_DE440_PATH not set - ctest entries skipped "
                   "(executables still built, run them manually)")
endif()
