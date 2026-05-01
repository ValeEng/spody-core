# raw_data

This folder contains external scientific datasets required by the Spody Core library.
The data files are **not included in the repository** due to their size. Follow the instructions below to populate each subfolder.

---

## Folder structure

```
raw_data/
├── DE440/          ← JPL planetary ephemeris (ASCII source files + generated binary)
└── GRGM1200A/      ← GRAIL lunar gravity model coefficients
```

---

## DE440 — JPL Planetary Ephemeris

Used by: `spody_ephemeris` module (`spody_setup_MappedEphemerisData`, `spody_setup_MappedEphemeris`, `spody_get_ephposition`, `spody_get_lunarlibrationangles`)

### 1. Download the ASCII files

Go to the JPL FTP server:

```
https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/
```

Download the following files into `raw_data/DE440/`:

| File | Coverage |
|---|---|
| `header.440` | Header (mandatory) |
| `ascp01550.440` | 1550 – 1650 |
| `ascp01650.440` | 1650 – 1750 |
| `ascp01750.440` | 1750 – 1850 |
| `ascp01850.440` | 1850 – 1950 |
| `ascp01950.440` | 1950 – 2050 |
| `ascp02050.440` | 2050 – 2150 |
| `ascp02150.440` | 2150 – 2250 |
| `ascp02250.440` | 2250 – 2350 |
| `ascp02350.440` | 2350 – 2450 |
| `ascp02450.440` | 2450 – 2550 |
| `ascp02550.440` | 2550 – 2650 |

> You only need the files covering the time range of your simulation. `ascp01950.440` is enough for most modern-era scenarios.

### 2. Generate the Spody binary file

The ASCII files must be converted once into the internal binary format (`.spody`) used by the memory-mapped ephemeris system. The TVB suite ships a one-shot helper that does it for you:

```powershell
# from spody-core/, after building with -DSPODY_BUILD_TVB=ON
build\tvb\Release\gen_de440_spody.exe
```

This produces `raw_data/DE440/de440.spody` (~100 MB) starting from all the JPL ASCII chunks present under `raw_data/DE440/`. Run it once.

If you prefer to do it from your own code, the underlying API is:

```c
#include "spody_core.h"

const char *path     = "raw_data/DE440";
const char *dates[]  = {"01550","01650","01750","01850","01950",
                         "02050","02150","02250","02350","02450","02550"};
const int   n_files  = 11;

int ret = spody_createfile_MappedEphemerisData(path, dates, n_files, "440");
```

After that, load the binary directly. The library uses a two-step setup: a shared, read-only `MappedEphemerisData` (one per process) and a per-thread `MappedEphemeris` query handle that holds a private Chebyshev cache.

```c
// Setup once (e.g. on the main thread): shared, read-only
MappedEphemerisData med = {0};
spody_setup_MappedEphemerisData(&med, "raw_data/DE440/de440.spody");

// Per worker thread: bind a private handle to the shared data
MappedEphemeris map = {0};
spody_setup_MappedEphemeris(&map, &med);

// Query Moon position wrt Earth at a given Julian Date
double position[3];
spody_get_ephposition(&map, 399, 301, jd_epoch, position);

// Free in reverse order when done
spody_free_MappedEphemeris(&map);
spody_free_MappedEphemerisData(&med);
```

> **Threading:** the same `MappedEphemerisData` can be safely shared across threads; each thread must own its own `MappedEphemeris`. Never share a `MappedEphemeris` across threads.

**Body index reference (JPL NAIF standard):**

| Index | Body |
|---|---|
| 399 | Earth |
| 301 | Moon |
| 10 | Sun |
| 499 | Mars |
| 299 | Venus |
| 599 | Jupiter |

---

## GRGM1200A — GRAIL Lunar Gravity Model

Used by: `spody_harmonics` module (`spody_load_HarmonicGravityData`, `spody_setup_HarmonicGravity`, `spody_get_hgaccbodyfixed`, `spody_get_hgaccbodyfixed_hpc`)

### Download

Go to the NASA PGDA page:

```
https://pgda.gsfc.nasa.gov/products/50
```

Download and place in `raw_data/GRGM1200A/`:

| File | Description |
|---|---|
| `gggrx_1200a_sha.tab` | Spherical harmonic coefficients (C, S) up to degree/order 1200 |
| `gggrx_1200a_sha.lbl` | PDS label file (metadata) |

### Usage

```c
#include "spody_core.h"

// Shared, read-only: load once per process
HarmonicGravityData hgd = {0};
spody_load_HarmonicGravityData(&hgd, "raw_data/GRGM1200A/gggrx_1200a_sha.tab", 100);

// Per-thread scratch buffers
HarmonicGravity hg = {0};
spody_setup_HarmonicGravity(&hg, &hgd);

// Compute gravitational acceleration in the Moon body-fixed (PA) frame [km, km/s^2]
double pos[3] = { /* lunar PA frame, km */ };
double acc[3] = {0};
spody_get_hgaccbodyfixed_hpc(&hg, pos, acc);   // production hot path

// or, for regression / audit:
// spody_get_hgaccbodyfixed(&hg, pos, acc);    // reference variant

spody_free_HarmonicGravity(&hg);
spody_free_HarmonicGravityData(&hgd);
```

> **Threading:** like the ephemeris module, the same `HarmonicGravityData` is safely shared across threads; each thread must own its own `HarmonicGravity` handle (because of the rolling-row scratch buffers).

> **Degree selection:** N=80–100 is a good balance between accuracy and speed for most mission scenarios. Use N≥200 for high-fidelity low-lunar-orbit propagation. Cost scales as ~N² (verified empirically: N=100 → ~10 µs/call, N=200 → ~34 µs/call on GCC/Ryzen-class CPU with `SPODY_ENABLE_OMP_SIMD=ON`).

> **`_hpc` vs reference:** `spody_get_hgaccbodyfixed_hpc` is the production-grade variant (branch-free + peeled inner loops + restrict-tagged pointers + `#pragma omp simd reduction` on GCC/Clang). It produces results bit-equivalent to the reference (rel error < 1e-12) and is ~1.5–1.7× faster on GCC/Clang at N≥50. Use it inside the RHS of a propagator; keep the reference for regression tests and audit.

---

## Summary

| Folder | File(s) needed | Source |
|---|---|---|
| `DE440/` | `header.440` + `ascp0XXXX.440` + generated `de440.spody` | [JPL FTP](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/) |
| `GRGM1200A/` | `gggrx_1200a_sha.tab` + `gggrx_1200a_sha.lbl` | [NASA PGDA](https://pgda.gsfc.nasa.gov/products/50) |