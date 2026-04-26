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

The ASCII files must be converted once into the internal binary format (`.spody`) used by the memory-mapped ephemeris system. Call this function once in your code before any simulation:

```c
#include "spody_core.h"

const char *path     = "raw_data/DE440";
const char *dates[]  = {"01550","01650","01750","01850","01950",
                         "02050","02150","02250","02350","02450","02550"};
const int   n_files  = 11;

// Generates: raw_data/DE440/de440.spody
int ret = spody_createfile_MappedEphemerisData(path, dates, n_files, "440");
```

This writes `de440.spody` into `raw_data/DE440/`. Run it **once** — after that, load the binary directly. The library uses a two-step setup: a shared, read-only `MappedEphemerisData` (one per process) and a per-thread `MappedEphemeris` query handle that holds a private Chebyshev cache.

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

Used by: `spody_harmonics` module (`spody_load_HarmonicGravityData`, `spody_setup_HarmonicGravity`)

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

HarmonicGravityData hgd;
HarmonicGravity     hg;

// Load up to degree N (max 1200, higher = more accurate, slower)
spody_load_HarmonicGravityData(&hgd, "raw_data/GRGM1200A/gggrx_1200a_sha.tab", 80);
spody_setup_HarmonicGravity(&hg, &hgd);

// Compute gravitational acceleration at a position in the Moon body-fixed frame [km, km/s^2]
double pos[3] = { ... };   // position in lunar PA frame [km]
double acc[3] = {0};
compute_harmonic_lunar_gravity_hpc(&hg, pos, acc);

// Free when done
spody_free_HarmonicGravity(&hg);
spody_free_HarmonicGravityData(&hgd);
```

> **Degree selection:** degree 80 is a good balance between accuracy and speed for most mission scenarios. Use degree ≥ 200 for high-fidelity low-lunar-orbit propagation.

---

## Summary

| Folder | File(s) needed | Source |
|---|---|---|
| `DE440/` | `header.440` + `ascp0XXXX.440` + generated `de440.spody` | [JPL FTP](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/) |
| `GRGM1200A/` | `gggrx_1200a_sha.tab` + `gggrx_1200a_sha.lbl` | [NASA PGDA](https://pgda.gsfc.nasa.gov/products/50) |