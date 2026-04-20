# 🛰️ Spody Core

> **Sp**ace **O**rbital **Dy**namics — A high-performance C library for astrodynamics and space mechanics.

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

---

## Features

| Module | Description |
|---|---|
| `ephemeris` | Parses and queries JPL ASCII ephemeris files for planetary positions and velocities (e.g. DE440) |
| `eclipse` | Detects solar eclipse conditions (umbra/penumbra) for orbiting objects |
| `harmonics` | Computes gravitational acceleration using spherical harmonic coefficients (e.g. GRGM1200A) |
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
│   ├── spody_eclipse.h
│   ├── spody_ephemeris.h
│   ├── spody_harmonics.h
│   ├── spody_mapping.h
│   └── spody_math.h
├── src/                    # Implementation files (.c)
├── tvb/                    # Test / validation benchmarks
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
│   ├── spody_eclipse.h
│   ├── spody_ephemeris.h
│   ├── spody_harmonics.h
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

**Example:** maximum-speed build with fast math enabled:

```bash
cmake .. -DSPODY_FAST_MATH=ON
cmake --build . --config Release
```

**Example:** safer numerical build with whole-program optimization disabled:

```bash
cmake .. -DSPODY_WHOLE_PROGRAM_OPT=OFF
cmake --build . --config Release
```

---

## Usage

Regardless of the integration method, include the umbrella header:

```c
#include "spody_core.h"
```

This exposes the full public API (ephemeris, eclipse, harmonics, math, mapping).

---

## License

Distributed under the **Apache License 2.0**.
See [`LICENSE`](./LICENSE) for full details.

---

## Author

**ValeEng** — [github.com/ValeEng](https://github.com/ValeEng)

---

*Built for precision. Designed for space.*
