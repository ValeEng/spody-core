# 🛰️ Spody Core

> SpOdy (**S**imultaneous **P**ropagation of **O**rbital **DY**namics) Core — A high-performance C library for astrodynamics and space mechanics.

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Language: C](https://img.shields.io/badge/Language-C99-lightgrey.svg)](https://en.wikipedia.org/wiki/C99)
[![Build: CMake](https://img.shields.io/badge/Build-CMake%20%E2%89%A5%203.10-brightgreen.svg)](https://cmake.org/)

---

## Overview

**Spody Core** is a lightweight, dependency-free C library designed for high-precision space dynamics and astrodynamics calculations. It enables simultaneous propagation of multiple independent space objects with a focus on performance and numerical accuracy.

The library provides a clean, modular API covering the core pillars of orbital mechanics:

- 🌍 **Ephemeris parsing** — ASCII planetary data from JPL
- 🌑 **Eclipse detection** — Umbra and penumbra modeling
- 🌐 **Spherical harmonics gravity** — High-fidelity lunar gravity
- 🚀 **Numerical integrators** — Adaptive Dormand-Prince 5(4) and classical RK4

The runtime API is **thread-safe by construction**: shared, read-only data structures
(ephemeris, gravity coefficients) are decoupled from per-thread query handles, so a
single dataset can drive many concurrent propagations without contention.

---

## Features

| Module | Description |
|---|---|
| `ephemeris` | Parses and queries JPL ASCII ephemeris files for planetary positions and velocities (e.g. DE440). Memory-mapped binary format with thread-safe handle/data split. |
| `eclipse` | Detects solar eclipse conditions (umbra/penumbra) for orbiting objects |
| `harmonics` | Computes gravitational acceleration using spherical harmonic coefficients (e.g. GRGM1200A). Includes a `_hpc` SIMD-friendly variant for production hot paths. |
| `integrators` | ODE integrators with a generic RHS callback. RKDP45 adaptive (validated, GMAT-style step control) and RK4 fixed-step are implemented; Verlet/RK78 are reserved. |
| `math` | Shared math utilities (rotation matrices, vector operations) |
| `mapping` | Cross-platform memory-mapped file I/O |

---

## External Data Sources

Spody Core relies on standard, publicly available scientific datasets:

- **JPL DE440 Ephemeris (ASCII)**
  → [https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/)

- **GRGM1200A Gravitational Model**
  Lunar gravitational harmonics coefficients up to degree/order 1200.
  → [https://pgda.gsfc.nasa.gov/products/50](https://pgda.gsfc.nasa.gov/products/50)

---

## Project Structure

```
spody-core/
├── include/                # Public headers — include this folder in your build
│   ├── spody_core.h        # Umbrella header (include this one)
│   ├── spody_const.h       # Physical and numerical constants
│   ├── spody_eclipse.h
│   ├── spody_ephemeris.h
│   ├── spody_harmonics.h
│   ├── spody_integrators.h
│   ├── spody_mapping.h
│   └── spody_math.h
├── src/                    # Implementation files (.c)
├── tvb/                    # Tests, validations, and benchmarks
├── raw_data/               # External datasets (DE440, GRGM1200A; see raw_data/README.md)
├── CMakeLists.txt
├── LICENSE
└── README.md
```

---

## Prerequisites

- A C99-compatible compiler (GCC, Clang, or MSVC)
- **CMake** ≥ 3.10 *(only needed for Options A and B below)*
- On Unix/Linux/macOS: the standard math library (`libm`) — linked automatically

---

## Integration

There are three ways to use Spody Core in your project. Pick the one that fits your workflow.

### Option A — Build standalone and link against the compiled library

Use this if you want a precompiled `.a`/`.lib` you can link against, or if you're distributing the library separately from your app.

**1. Clone and build**

```bash
git clone https://github.com/ValeEng/spody-core.git
cd spody-core
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

> **Note on `--config`:**
> Multi-config generators (Visual Studio, Xcode) pick the build type at *build* time via `--config`. Single-config generators (Makefile, Ninja) pick it at *configure* time and default to `Release` — you can override with `cmake -DCMAKE_BUILD_TYPE=Debug ..`.

This produces:
- `build/Release/spody_core.lib` on Windows (MSVC)
- `build/libspody_core.a` on Linux / macOS

**2. (Optional) Install to a clean distribution folder**

```bash
cmake --install . --prefix ./dist
```

This creates:

```
dist/
├── include/
│   ├── spody_core.h
│   ├── spody_const.h
│   ├── spody_eclipse.h
│   ├── spody_ephemeris.h
│   ├── spody_harmonics.h
│   ├── spody_integrators.h
│   ├── spody_mapping.h
│   └── spody_math.h
└── lib/
    └── libspody_core.a    (or spody_core.lib on Windows)
```

**3. Compile your program**

```bash
gcc main.c -I./dist/include -L./dist/lib -lspody_core -lm -o my_space_app
```

---

### Option B — Use as a CMake sub-project (recommended)

Use this if your project is already CMake-based. Add `spody-core/` as a subdirectory (or a git submodule) of your project:

```cmake
add_subdirectory(external/spody-core)

target_link_libraries(my_target PRIVATE spody_core)
```

CMake handles include paths, the `libm` link, and build order automatically. No install step needed.

---

### Option C — Drop-in sources (no build system)

Use this for simple projects, Makefile-based builds, embedded targets, or when you just want to compile everything in one shot. Copy `include/` and `src/` into your project, then:

```bash
gcc main.c src/spody_*.c -Iinclude -lm -O2 -o my_space_app
```

Since all public headers live in `include/`, a single `-Iinclude` is enough.

---

## Build Options

The `CMakeLists.txt` exposes a few flags you can toggle on the `cmake` command line with `-D<n>=<VALUE>`.

| Option | Default | Description |
|---|:---:|---|
| `SPODY_FAST_MATH` | `OFF` | Enables aggressive floating-point optimizations (`/fp:fast` on MSVC, `-ffast-math` on GCC/Clang). Faster but may introduce small numerical drift over long integrations. Leave `OFF` for bit-reproducible results. |
| `SPODY_WHOLE_PROGRAM_OPT` | `ON` | Enables whole-program / link-time optimization in Release builds (`/GL + /LTCG` on MSVC, `-flto` on GCC/Clang). Slower link, faster runtime. |
| `SPODY_SILENCE_MSVC_CRT_WARNINGS` | `ON` | On MSVC, silences warnings about "unsafe" standard C functions (`strtok`, `sprintf`, ...). These are portable ISO C functions; the `_s` alternatives are non-portable Microsoft extensions. |
| `SPODY_ENABLE_OMP_SIMD` | `OFF` | Activates `#pragma omp simd` hints in the harmonic-gravity hot loops (`/openmp:experimental` on MSVC, `-fopenmp-simd` on GCC/Clang). On GCC/Clang this enables explicit SIMD reduction for ~1.5× speedup at degree N≥50. On MSVC the macros expand to no-ops (front-end ignores `reduction` on `simd`); the flag is still useful because it changes the auto-vectorizer behaviour. |
| `SPODY_VEC_REPORT` | `OFF` | Emits the compiler's auto-vectorization diagnostic at build time (`/Qvec-report:2` on MSVC, `-fopt-info-vec-all` on GCC, `-Rpass=loop-vectorize` on Clang). Useful while tuning hot loops. |
| `SPODY_BUILD_TVB` | `OFF` | Builds the tests, validations, and benchmarks under [`tvb/`](tvb/). See [tvb/README.md](tvb/README.md). |

**Example:** production build on MSVC (default):

```bash
cmake -B build
cmake --build build --config Release
```

**Example:** production build on GCC with SIMD pragmas active (recommended for hot-path harmonic gravity):

```bash
cmake -B build-gcc -G Ninja -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DSPODY_ENABLE_OMP_SIMD=ON
cmake --build build-gcc
```

**Example:** safer numerical build (no fast-math, no LTO):

```bash
cmake -B build -DSPODY_FAST_MATH=OFF -DSPODY_WHOLE_PROGRAM_OPT=OFF
cmake --build build --config Release
```

### Compiler recommendation

For the harmonic-gravity hot path the cost-model differences between compilers are
visible. On MSVC the `_hpc` variant gives roughly 5–10% speedup; on GCC (or Clang)
with `SPODY_ENABLE_OMP_SIMD=ON` it gives 1.5× at degree 100, 1.7× at degree 200 over
the reference implementation. For production high-fidelity propagation we therefore
recommend GCC or Clang with `SPODY_ENABLE_OMP_SIMD=ON`. MSVC remains fully supported
and produces identical numerical results.

---

## Usage

Regardless of the integration method, include the umbrella header:

```c
#include "spody_core.h"
```

This exposes the full public API (ephemeris, eclipse, harmonics, integrators, math, mapping).

### Quick tour

**Ephemeris** (read-only data shared across threads, per-thread query handle):

```c
MappedEphemerisData med = {0};
spody_setup_MappedEphemerisData(&med, "raw_data/DE440/de440.spody");

MappedEphemeris map = {0};                     // one per thread
spody_setup_MappedEphemeris(&map, &med);

double r[3];
spody_get_ephposition(&map, /*central=*/399, /*target=*/301, jd, r);  // Earth -> Moon

spody_free_MappedEphemeris(&map);
spody_free_MappedEphemerisData(&med);
```

**Harmonic gravity** (read-only coefficients shared, per-thread scratch buffers):

```c
HarmonicGravityData hgd = {0};
spody_load_HarmonicGravityData(&hgd, "raw_data/GRGM1200A/gggrx_1200a_sha.tab", 100);

HarmonicGravity hg = {0};                       // one per thread
spody_setup_HarmonicGravity(&hg, &hgd);

double pos_pa[3] = { /* lunar PA frame, km */ };
double acc[3];
spody_get_hgaccbodyfixed_hpc(&hg, pos_pa, acc); // production hot path

spody_free_HarmonicGravity(&hg);
spody_free_HarmonicGravityData(&hgd);
```

The unsuffixed `spody_get_hgaccbodyfixed` is the algorithmic reference used as
regression baseline; `_hpc` produces bit-equivalent results (rel error < 1e-12)
and is the variant to call from a propagator's RHS.

**Numerical integrator** (RKDP45 adaptive on a 6-state two-body propagation):

```c
typedef struct { double mu; } Params;

int rhs(double t, const double *y, double *dy, void *user) {
    (void)t;
    Params *p = (Params*)user;
    double r2 = y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
    double k  = -p->mu / (r2 * sqrt(r2));
    dy[0] = y[3]; dy[1] = y[4]; dy[2] = y[5];
    dy[3] = k*y[0]; dy[4] = k*y[1]; dy[5] = k*y[2];
    return 0;
}

Params params = { .mu = EARTH_MU };
IntegratorAllData integ;
IntegratorOptions opt;
spody_default_integrator_options(SPODY_INTEG_RK45, &opt);

spody_setup_integrator(&integ, SPODY_INTEG_RK45, &opt, /*dim=*/6, rhs, &params);

double y0[6] = { 7000.0, 0,0,  0, 7.546, 0 };
spody_set_integrator_state(&integ, /*t0=*/0.0, y0);
spody_propagate_untilend(&integ, /*t_end=*/5400.0);   // one LEO orbit

spody_free_integrator(&integ);
```

The integrator carries no global state; the per-thread `IntegratorAllData` holds
the state vector, the RHS callback, an opaque user payload (typically containing
shared pointers to ephemeris + gravity + per-thread handles), and the scratch
buffers for the Runge-Kutta stages.

---

## License

Distributed under the **Apache License 2.0**.
See [`LICENSE`](./LICENSE) for full details.

---

## Author

**ValeEng** — [github.com/ValeEng](https://github.com/ValeEng)

---

*Built for precision. Designed for space.*
