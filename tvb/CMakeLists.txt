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

add_executable(bench_integrators benchmarks/bench_integrators.c)
target_link_libraries(bench_integrators PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(bench_integrators PRIVATE m)
endif()

add_executable(bench_integrator_physics benchmarks/bench_integrator_physics.c)
target_link_libraries(bench_integrator_physics PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(bench_integrator_physics PRIVATE m)
endif()

add_executable(bench_harmonics_hpc benchmarks/bench_harmonics_hpc.c)
target_link_libraries(bench_harmonics_hpc PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(bench_harmonics_hpc PRIVATE m)
endif()

# ---- One-shot tools ----
add_executable(gen_de440_spody benchmarks/gen_de440_spody.c)
target_link_libraries(gen_de440_spody PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(gen_de440_spody PRIVATE m)
endif()

# ---- Tests ----
add_executable(test_ephemeris_moon_sun tests/test_ephemeris_moon_sun.c)
target_link_libraries(test_ephemeris_moon_sun PRIVATE spody_core)

add_executable(test_integrator_kepler tests/test_integrator_kepler.c)
target_link_libraries(test_integrator_kepler PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(test_integrator_kepler PRIVATE m)
endif()

add_executable(test_harmonics_hpc tests/test_harmonics_hpc.c)
target_link_libraries(test_harmonics_hpc PRIVATE spody_core)
if(UNIX OR APPLE)
    target_link_libraries(test_harmonics_hpc PRIVATE m)
endif()
if(SPODY_GRGM1200A_PATH)
    add_test(NAME harmonics_hpc COMMAND test_harmonics_hpc "${SPODY_GRGM1200A_PATH}")
endif()

# Register ctest entries.
# test_integrator_kepler is self-contained (no external data), always runs.
add_test(NAME integrator_kepler COMMAND test_integrator_kepler)

if(SPODY_DE440_PATH)
    add_test(NAME ephemeris_moon_sun
             COMMAND test_ephemeris_moon_sun "${SPODY_DE440_PATH}")
else()
    message(STATUS "SPODY_DE440_PATH not set - ephemeris ctest entries skipped "
                   "(executables still built, run them manually)")
endif()
