# Tests, validations and benchmarks

This directory is **OFF by default**. It only builds when CMake is configured
with `-DSPODY_BUILD_TVB=ON`, so consumers of the library never pay for it.

## Layout

```
tvb/
├── benchmarks/
│   ├── bench_ephemeris_batch.c        single-call vs batch vs no-cache ephemeris query
│   ├── bench_harmonics_hpc.c          ref vs hpc harmonic gravity at varying degree
│   ├── bench_integrators.c            RK4 (multiple h) vs RKDP45 (multiple tol) on Kepler
│   ├── bench_integrator_physics.c     RHS cost breakdown: 2-body, 3rd body, SRP, harmonics 10..200
│   └── gen_de440_spody.c              one-shot tool: builds raw_data/DE440/de440.spody from ASCII
├── tests/
│   ├── test_ephemeris_moon_sun.c      Moon→Sun regression + triangular closure invariants
│   ├── test_harmonics_hpc.c           bit-equivalence regression: hpc vs reference (rel < 1e-12)
│   └── test_integrator_kepler.c       RK4 + RKDP45 on a circular LEO: closure & energy conservation
└── validation/                        (reserved, not populated yet)
```

## Build

### Multi-config generators (Visual Studio, Xcode)

One build tree contains all configs. Pick the config at build time.

```bash
cmake -B build -DSPODY_BUILD_TVB=ON
cmake --build build --config Release
```

All TVB executables land in a flat `Release/` sub-folder under `build/tvb/`:

```
build/tvb/Release/bench_ephemeris_batch.exe
build/tvb/Release/bench_harmonics_hpc.exe
build/tvb/Release/bench_integrators.exe
build/tvb/Release/bench_integrator_physics.exe
build/tvb/Release/gen_de440_spody.exe
build/tvb/Release/test_ephemeris_moon_sun.exe
build/tvb/Release/test_harmonics_hpc.exe
build/tvb/Release/test_integrator_kepler.exe
```

For Debug builds replace `Release` with `Debug` everywhere.

### Single-config generators (Ninja, Unix Makefiles, MinGW Makefiles)

One build tree = one config. Pick the config at configure time.

```bash
cmake -B build-gcc -G Ninja -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release \
                  -DSPODY_BUILD_TVB=ON -DSPODY_ENABLE_OMP_SIMD=ON
cmake --build build-gcc
```

No config sub-folder; executables land directly in `build-gcc/tvb/`. Use a separate
build tree per generator (e.g. `build/` for MSVC, `build-gcc/` for GCC) to avoid
clobbering CMake cache files.

Need both Debug and Release? Use two separate build trees
(e.g. `build-release/` and `build-debug/`).

## Running the benchmarks

### Generate the DE440 binary (one-shot, first time only)

If `raw_data/DE440/` contains the JPL ASCII files (`header.440`, `ascpXXXXX.440`)
but no `de440.spody`, build and run the helper from the repository root:

```powershell
# from spody-core/
build\tvb\Release\gen_de440_spody.exe
```

This produces `raw_data/DE440/de440.spody` (~100 MB) covering all the supplied
ASCII chunks. Run it only once.

### Ephemeris batch (no external data needed beyond `de440.spody`)

```powershell
build\tvb\Release\bench_ephemeris_batch.exe raw_data\DE440\de440.spody
```

Compares single-call vs batch vs no-cache ephemeris queries.

### Harmonic gravity hpc

```powershell
build\tvb\Release\bench_harmonics_hpc.exe raw_data\GRGM1200A\gggrx_1200a_sha.tab
```

Reports ns/call for the reference and `_hpc` variants of `spody_get_hgaccbodyfixed`
at degree 10/50/100/200. With GCC + `SPODY_ENABLE_OMP_SIMD=ON` expect ~1.5× at N=100.

### Integrators (Kepler 2-body)

```powershell
build\tvb\Release\bench_integrators.exe
```

Sweeps RK4 fixed step `h ∈ {30, 10, 5, 1}` s and RKDP45 adaptive
`rel_tol ∈ {1e-6, 1e-9, 1e-12}` over 10 LEO orbits. Reports steps, RHS calls,
wall time, and final energy drift. Self-contained, no external data.

### Integrator physics breakdown (RHS cost analysis)

```powershell
build\tvb\Release\bench_integrator_physics.exe raw_data\DE440\de440.spody raw_data\GRGM1200A\gggrx_1200a_sha.tab
```

Propagates a 100 km LLO for 5 orbits with progressively heavier physics
(2-body → +3rd body → +SRP → spherical harmonics N=10/50/100/200 → full).
Each harmonics scenario runs twice (`ref` + `hpc`) for a side-by-side comparison
of the reference and high-performance implementations. Useful to quantify how
much wall time the RHS dominates in real propagation.

## Running the tests

Manually (works on both generator types):

```powershell
build\tvb\Release\test_ephemeris_moon_sun.exe raw_data\DE440\de440.spody
build\tvb\Release\test_harmonics_hpc.exe       raw_data\GRGM1200A\gggrx_1200a_sha.tab
build\tvb\Release\test_integrator_kepler.exe   # self-contained, no data file
```

Exit code `0` means PASS, non-zero means FAIL and the failing checks are
printed to stdout.

Via `ctest` (auto-registers if you pass the data paths at configure time):

```bash
cmake -B build -DSPODY_BUILD_TVB=ON \
               -DSPODY_DE440_PATH=$(pwd)/raw_data/DE440/de440.spody \
               -DSPODY_GRGM1200A_PATH=$(pwd)/raw_data/GRGM1200A/gggrx_1200a_sha.tab
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Tests that don't need external data (`integrator_kepler`) are always registered.

## Locating the executables without remembering the generator

If you forget which generator you used, let the OS find the binary:

```powershell
# PowerShell
Get-ChildItem build -Recurse -Filter bench_ephemeris_batch*.exe
```

```bash
# Linux / macOS / MSYS2
find build -name 'bench_ephemeris_batch*'
```

## Writing a new test or benchmark

1. Add a self-contained C file with a `main()` under `tests/` or `benchmarks/`
   (or `validation/` if it requires an external reference dataset).
2. In `tvb/CMakeLists.txt`, add:

   ```cmake
   add_executable(my_new_thing tests/my_new_thing.c)
   target_link_libraries(my_new_thing PRIVATE spody_core)
   ```

3. If it's a test, optionally register it with ctest:

   ```cmake
   if(SPODY_DE440_PATH)
       add_test(NAME my_new_thing
                COMMAND my_new_thing "${SPODY_DE440_PATH}")
   endif()
   ```

4. Re-run cmake configure + build.

> **No framework dependency**: each test is just a `main()` that returns `0`
> on success and non-zero on failure.
