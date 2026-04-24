# Tests, validations and benchmarks

This directory is **OFF by default**. It only builds when CMake is configured
with `-DSPODY_BUILD_TVB=ON`, so consumers of the library never pay for it.

## Layout

```
tvb/
├── benchmarks/
│   └── bench_ephemeris_batch.c      single-call vs batch vs no-cache
├── tests/
│   └── test_ephemeris_moon_sun.c    Moon→Sun regression + invariants
└── validation/                      (reserved, not populated yet)
```

## Build

### Multi-config generators (Visual Studio, Xcode)

One build tree contains all configs. Pick the config at build time.

```bash
cmake -B build -DSPODY_BUILD_TVB=ON
cmake --build build --config Release
```

Executables land in a config sub-folder:

```
build/tvb/benchmarks/Release/bench_ephemeris_batch.exe
build/tvb/tests/Release/test_ephemeris_moon_sun.exe
```

For Debug builds replace `Release` with `Debug` everywhere.

### Single-config generators (Ninja, Unix Makefiles, MinGW Makefiles)

One build tree = one config. Pick the config at configure time.

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODY_BUILD_TVB=ON
cmake --build build
```

No config sub-folder; executables land directly in:

```
build/tvb/benchmarks/bench_ephemeris_batch
build/tvb/tests/test_ephemeris_moon_sun
```

Need both Debug and Release? Use two separate build trees
(e.g. `build-release/` and `build-debug/`).

## Running the benchmark

It needs a compiled DE440 binary (`*.spody`) as its only argument.

```powershell
# Visual Studio / multi-config
build\tvb\benchmarks\Release\bench_ephemeris_batch.exe path\to\de440.spody
```

```bash
# Ninja / Makefile / single-config
./build/tvb/benchmarks/bench_ephemeris_batch path/to/de440.spody
```

Expected output (numbers are indicative):

```
=== Ephemeris query benchmark ===
Iterations      : 100000 dates * 4 bodies = 400000 queries

A) single calls, no cache   : 1234.56 ms  (  3086 ns / query)
B) single calls, with cache :  678.90 ms  (  1697 ns / query)
C) batch API, with cache    :  456.78 ms  (  1142 ns / query)

speedup B vs A (cache)      : 1.82x
speedup C vs A (cache+batch): 2.70x
speedup C vs B (batch only) : 1.49x
```

## Running the tests

Manually (works on both generator types):

```powershell
# multi-config
build\tvb\tests\Release\test_ephemeris_moon_sun.exe path\to\de440.spody
```

```bash
# single-config
./build/tvb/tests/test_ephemeris_moon_sun path/to/de440.spody
```

Exit code `0` means PASS, non-zero means FAIL and the failing checks are
printed to stdout.

Via `ctest` (only if you passed `-DSPODY_DE440_PATH=<abs path>` at configure):

```bash
cmake -B build -DSPODY_BUILD_TVB=ON -DSPODY_DE440_PATH=C:/data/de440.spody
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Locating the executables without remembering the generator

If you forget which generator you used, let the OS find the binary:

```powershell
# PowerShell
Get-ChildItem build -Recurse -Filter bench_ephemeris_batch*.exe
```

```cmd
:: cmd
dir build\*bench_ephemeris_batch*.exe /s /b
```

```bash
# Linux / macOS
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
