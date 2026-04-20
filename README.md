# 🛰️ Spody Core

> **Sp**ace **O**rbital **Dy**namics — A high-performance C library for astrodynamics and space mechanics.

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Language: C](https://img.shields.io/badge/Language-C99-lightgrey.svg)](https://en.wikipedia.org/wiki/C99)
[![Build: CMake](https://img.shields.io/badge/Build-CMake%20%E2%89%A5%203.10-brightgreen.svg)](https://cmake.org/)

---

## Overview

**Spody Core** is a lightweight, dependency-free C library designed for high-precision space dynamics and astrodynamics calculations. It enables simultaneous propagation of multiple independent space objects with a focus on performance and numerical accuracy.

The library provides a clean, modular API covering the core pillars of orbital mechanics:

- 🌍 **Ephemeris parsing** — DE440 ASCII planetary data from JPL
- 🌑 **Eclipse detection** — Umbra and penumbra modeling
- 🌐 **Spherical harmonics gravity** — High-fidelity lunar gravity via GRGM1200A

---

## Features

| Module | Description |
|---|---|
| `ephemeris` | Parses and queries JPL DE440 ASCII ephemeris files for planetary positions and velocities |
| `eclipse` | Detects solar eclipse conditions (umbra/penumbra) for orbiting objects |
| `harmonics` | Computes gravitational acceleration using spherical harmonic coefficients (GRGM1200A) |
| `spody_mapping` | Internal coordinate and state-vector mapping utilities |

---

## External Data Sources

Spody Core relies on standard, publicly available scientific datasets:

- **JPL DE440 Ephemeris (ASCII)**
  Planetary positions and velocities.
  → [https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/)

- **GRGM1200A Gravitational Model**
  Lunar gravitational harmonics coefficients up to degree/order 1200.
  → [https://pgda.gsfc.nasa.gov/products/50](https://pgda.gsfc.nasa.gov/products/50)

---

## Project Structure

```
spody-core/
├── include/
│   └── spody_core.h        # Public umbrella header — include this in your project
├── src/
│   ├── eclipse.c / .h      # Eclipse detection logic
│   ├── ephemeris.c / .h    # Ephemeris parsing and interpolation
│   ├── harmonics.c / .h    # Spherical harmonics gravity model
│   └── spody_mapping.c / .h # Coordinate mapping utilities
├── tvb/                    # Test/validation benchmarks
├── CMakeLists.txt
├── LICENSE
└── README.md
```

---

## Prerequisites

- A C99-compatible compiler (GCC, Clang, or MSVC)
- **CMake** ≥ 3.10
- On Unix/Linux/macOS: the standard math library (`libm`) — linked automatically

---

## Build Instructions

### 1. Clone the repository

```bash
git clone https://github.com/ValeEng/spody-core.git
cd spody-core
```

### 2. Configure and compile

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

This produces `libspody_core.a` (Unix) or `spody_core.lib` (Windows) inside the `build/` directory.

### 3. Install (optional)

To generate a clean distribution package with headers and the compiled library:

```bash
cmake --install . --prefix ./dist
```

This creates:
```
dist/
├── include/
│   ├── spody_core.h
│   ├── eclipse.h
│   ├── ephemeris.h
│   ├── harmonics.h
│   └── spody_mapping.h
└── lib/
    └── libspody_core.a
```

---

## Usage

Include the umbrella header in your project:

```c
#include "spody_core.h"
```

Compile and link against the library:

```bash
gcc main.c -I./dist/include -L./dist/lib -lspody_core -lm -o my_space_app
```

Or, if integrating via CMake, add to your `CMakeLists.txt`:

```cmake
add_subdirectory(spody-core)
target_link_libraries(your_target PRIVATE spody_core)
```

---

## License

Distributed under the **Apache License 2.0**.
See [`LICENSE`](./LICENSE) for full details.

---

## Author

**ValeEng** — [github.com/ValeEng](https://github.com/ValeEng)

---

*Built for precision. Designed for space.*
