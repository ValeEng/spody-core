# raw_data

This folder contains external scientific datasets required by the Spody Core library.
The data files are **not included in the repository** due to their size. Follow the instructions below to populate each subfolder.

---

## Folder structure

```
raw_data/
├── DE440/          ← JPL planetary ephemeris (ASCII source files + generated binary)
├── GRGM1200A/      ← GRAIL lunar gravity model, original release
└── GRGM1200B/      ← GRAIL lunar gravity model, updated release (recommended)
```

A **Moon-centred** propagation needs only the two above. An
**Earth-centred** one additionally needs a terrestrial gravity field,
the Earth-orientation inputs, and — with drag on — space weather; see
[Earth-centred datasets](#earth-centred-datasets) below. Those are not
lunar-mission assets and are usually kept in the consuming
application's own data directory rather than here (the SpOdy app's
Setup wizard downloads them into its portable `data/`), so this folder
carries no subdirectory for them by default.

---

## DE440 — JPL Planetary Ephemeris

Used by: `spody_ephemeris` module (`spody_setup_MappedEphemerisData`, `spody_setup_MappedEphemeris`, `spody_get_ephposition`, `spody_get_ephvelocity`, `spody_get_ephstate`, `spody_get_lunarlibrationangles`)

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

The ASCII files must be converted once into the internal binary format (`.spody`) used by the memory-mapped ephemeris system. On-disk layout:

```
| magic "SPDYEPET" (8 B) | format_version (4 B) | reserved (4 B) | header ... | records ... |
```

Epochs in the file are in **Ephemeris Time (ET) seconds past J2000** — the same scale used at the API boundary. This is ~250× more precise near today's epoch than the legacy JD-days format (ULP ~183 ns vs ~47 µs at JD ≈ 2.46e6).

Conversion is done programmatically via `spody_createfile_MappedEphemerisData`:

```c
#include "spody_core.h"

const char *path     = "raw_data/DE440";
const char *dates[]  = {"01550","01650","01750","01850","01950",
                         "02050","02150","02250","02350","02450","02550"};
const int   n_files  = 11;

int ret = spody_createfile_MappedEphemerisData(path, dates, n_files, "440");
```

This produces `raw_data/DE440/de440.spody` (~100 MB) from the JPL ASCII chunks. Run it once after downloading the source files; the resulting binary is portable across machines of the same endianness.

After that, load the binary directly. The library uses a two-step setup: a shared, read-only `MappedEphemerisData` (one per process) and a per-thread `MappedEphemeris` query handle that holds a private Chebyshev cache.

```c
// Setup once (e.g. on the main thread): shared, read-only
MappedEphemerisData med = {0};
spody_setup_MappedEphemerisData(&med, "raw_data/DE440/de440.spody");

// Per worker thread: bind a private handle to the shared data
MappedEphemeris map = {0};
spody_setup_MappedEphemeris(&map, &med);

// Query Moon position wrt Earth at a given ET epoch (seconds past J2000).
// If you have a Julian Date, convert with ET_FROM_JD(jd) from spody_const.h.
double position[3];
spody_get_ephposition(&map, 399, 301, et, position);

// Velocity (km/s) and full state [x,y,z,vx,vy,vz] come from the analytic
// derivative of the same Chebyshev series -- exact, no finite differences.
double velocity[3], state[6];
spody_get_ephvelocity(&map, 399, 301, et, velocity);
spody_get_ephstate(&map, 399, 301, et, state);

// Free in reverse order when done
spody_free_MappedEphemeris(&map);
spody_free_MappedEphemerisData(&med);
```

> **Threading:** the same `MappedEphemerisData` can be safely shared across threads; each thread must own its own `MappedEphemeris`. Never share a `MappedEphemeris` across threads.

**Body index reference (JPL NAIF standard):**

DE440 stores positions relative to the **Solar System Barycenter (SSB)**, plus the Moon relative to Earth-Moon Barycenter (EMB). `spody_get_ephposition(map, observer, target, et, r)` returns the position of `target` as seen from `observer`, handling the EMB/EMRAT bookkeeping internally; `spody_get_ephvelocity` / `spody_get_ephstate` apply the same bookkeeping to the km/s rates.

| NAIF ID | Body |
|:-:|---|
| 0 | Solar System Barycenter (SSB) |
| 3 | Earth-Moon Barycenter (EMB) |
| 10 | Sun |
| 199 / 1 | Mercury |
| 299 / 2 | Venus |
| 399 | Earth |
| 301 | Moon |
| 499 / 4 | Mars |
| 599 / 5 | Jupiter |
| 699 / 6 | Saturn |
| 799 / 7 | Uranus |
| 899 / 8 | Neptune |

---

## GRGM1200 — GRAIL Lunar Gravity Model

Used by: `spody_harmonics` module (`spody_load_HarmonicGravityData`, `spody_setup_HarmonicGravity`, `spody_get_hgaccbodyfixed`, `spody_get_hgaccbodyfixed_hpc`)

Two releases of the same GRAIL-derived model are supported:

| Release | Folder | Notes |
|---|---|---|
| **GRGM1200A** | `raw_data/GRGM1200A/` | Original release, paired in literature with many published validations. |
| **GRGM1200B** | `raw_data/GRGM1200B/` | Updated release with refined coefficients at high degree. Recommended for new work, used by the LRO validation bench. |

Both share the same on-disk format (PDS spherical-harmonic `.tab`), API, and degree range — pick one per simulation.

### Download

Go to the NASA PGDA page:

```
https://pgda.gsfc.nasa.gov/products/50
```

Download and place into the matching folder:

| File | Folder | Description |
|---|---|---|
| `gggrx_1200a_sha.tab` | `GRGM1200A/` | Spherical harmonic coefficients (C, S) up to degree/order 1200 |
| `gggrx_1200a_sha.lbl` | `GRGM1200A/` | PDS label file (metadata) |
| `gggrx_1200b_sha.tab` | `GRGM1200B/` | Updated coefficients (recommended) |
| `gggrx_1200b_sha.lbl` | `GRGM1200B/` | PDS label file (metadata) |

### Usage

```c
#include "spody_core.h"

// Shared, read-only: load once per process
HarmonicGravityData hgd = {0};
spody_load_HarmonicGravityData(&hgd, "raw_data/GRGM1200B/gggrx_1200b_sha.tab", 100);

// Per-thread scratch buffers
HarmonicGravity hg = {0};
spody_setup_HarmonicGravity(&hg, &hgd);

// Compute the disturbing acceleration in the Moon body-fixed (PA) frame [km, km/s^2].
// acc_out IS the acceleration (-grad V_pert, n>=2 only); the two-body central
// term is provided separately by spody_force_twobody().
double pos[3] = { /* lunar PA frame, km */ };
double acc[3] = {0};
spody_get_hgaccbodyfixed_hpc(&hg, pos, acc);   // production hot path

// or, for regression / audit:
// spody_get_hgaccbodyfixed(&hg, pos, acc);    // reference variant

spody_free_HarmonicGravity(&hg);
spody_free_HarmonicGravityData(&hgd);
```

> **Threading:** like the ephemeris module, the same `HarmonicGravityData` is safely shared across threads; each thread must own its own `HarmonicGravity` handle (because of the rolling-row scratch buffers).

> **Degree selection:** picked empirically from the LRO 6-day bench against an external reference.
>
> | N | use case |
> |:-:|---|
> | 30–50 | quick sanity propagation, low-fidelity orbit averaging |
> | 80 | reasonable default (sub-km vs SPICE LRO POD over 6 days) |
> | **150** | **sweet spot for LLO/LRO** — ~95% of the residual reduction of N=200 at half the cost |
> | 200 | high-fidelity floor; beyond ~200 the GRGM1200B coefficients become weakly observed and adding terms can slightly *increase* the mean drift |
>
> Cost scales as ~N². Indicative timings for a 6-day LRO propagation on a desktop x86-64 in Release: N=80 → ~3 s, N=150 → ~12 s, N=200 → ~23 s. Absolute numbers depend on CPU and compiler; the relative scaling holds.

> **`_hpc` vs reference:** `spody_get_hgaccbodyfixed_hpc` is the production-grade variant (branch-free + peeled inner loops + restrict-tagged pointers + `#pragma omp simd reduction` on GCC/Clang). It produces results bit-equivalent to the reference (rel error < 1e-12) and is ~1.5–1.7× faster on GCC/Clang at N≥50. Use it inside the RHS of a propagator; keep the reference for regression tests and audit.

---

## Earth-centred datasets

Only needed when the central body is the Earth. None of these is a
`.tab`-style drop-in: the gravity field needs a conversion step, and
the other three are read directly from their published text formats.

### EIGEN-6C4 — terrestrial gravity field

Used by: `spody_harmonics`, after conversion by `spody_icgem`.

Download the model in ICGEM `.gfc` form from
[https://icgem.gfz-potsdam.de/tom_longtime](https://icgem.gfz-potsdam.de/tom_longtime),
then convert it once to the engine's `.tab` — from C:

```c
spody_convert_icgem_to_tab("eigen-6c4.gfc", "eigen-6c4.tab", /*max_degree=*/0);
```

or, if you have the SpOdy application built on top of this library:

```bash
spody convert harmonics_icgem eigen-6c4.gfc eigen-6c4.tab
```

`max_degree <= 0` (equivalently: omitting `--max-degree`) keeps the
full field, and that is the recommended setting even though EIGEN-6C4
runs to degree 2190. The harmonics loader truncates again at read time
according to the caller's requested degree, so storing everything
leaves every future run free to pick its own truncation without a
re-conversion. Cap it only if the on-disk size is the binding
constraint.

### IERS Earth Orientation Parameters

Used by: `spody_eop`, feeding `spody_earth_orientation`.

One file, `finals2000A.all`, from
[https://datacenter.iers.org/products/eop/rapid/standard/](https://datacenter.iers.org/products/eop/rapid/standard/).
It carries polar motion, UT1−UTC and celestial-pole offsets, with the
recent span observed and a forward span predicted.

**It goes stale.** The predicted tail degrades as it ages, so a
long-lived install should re-download periodically rather than treat
this as a one-time asset.

### IAU 2006/2000A_R06 precession-nutation tables

Used by: `spody_earth_orientation` for the `(X, Y, s)` series.

The tabulated series from
[https://iers-conventions.obspm.fr/](https://iers-conventions.obspm.fr/).
Point the consuming application at the directory holding them; unlike
the EOP file these are fixed and never need refreshing.

### CelesTrak space weather

Used by: `spody_atmosphere` → `spody_nrlmsise00`, only when drag is on.

One file, `SW-All.csv`, from
[https://celestrak.org/SpaceData/](https://celestrak.org/SpaceData/):
observed daily F10.7 plus the 3-hour Ap history NRLMSISE-00 needs for
storm-time density. Like the EOP file it has an observed span and a
predicted tail, and re-downloading is what keeps a run near the
present epoch honest.

---

## Summary

| Dataset | File(s) needed | Source |
|---|---|---|
| `DE440/` | `header.440` + `ascp0XXXX.440` + generated `de440.spody` (SPDYEPET, ET seconds) | [JPL FTP](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/) |
| `GRGM1200A/` | `gggrx_1200a_sha.tab` + `gggrx_1200a_sha.lbl` | [NASA PGDA](https://pgda.gsfc.nasa.gov/products/50) |
| `GRGM1200B/` | `gggrx_1200b_sha.tab` + `gggrx_1200b_sha.lbl` (recommended) | [NASA PGDA](https://pgda.gsfc.nasa.gov/products/50) |
| EIGEN-6C4 *(Earth)* | `eigen-6c4.gfc` → converted `.tab` | [ICGEM](https://icgem.gfz-potsdam.de/tom_longtime) |
| IERS EOP *(Earth)* | `finals2000A.all` — refresh periodically | [IERS](https://datacenter.iers.org/products/eop/rapid/standard/) |
| IAU 2006 tables *(Earth)* | the `(X, Y, s)` series directory | [IERS Conventions](https://iers-conventions.obspm.fr/) |
| Space weather *(drag)* | `SW-All.csv` — refresh periodically | [CelesTrak](https://celestrak.org/SpaceData/) |