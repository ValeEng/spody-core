# SpOdy Core
**Simultaneous Propagation of Orbital DYnamics — Core Library**

SpOdy Core is a **high-performance, high-precision orbital dynamics library**
written in **C**, designed for the **simultaneous propagation of multiple
independent space objects**. It provides the numerical integrators, orbital
dynamics models, and propagation engine that can be used in research, simulation,
and engineering applications.

This library focuses on **modularity, efficiency, and scientific rigor**, while
keeping each trajectory dynamically independent (no mutual perturbations are
currently modeled).

---

## Key Features
- **Main kernels**
  - Thread safe ephemeris kernel on DE440
  - Thread safe gravitational armonics kernel on GRGM1200A

- **High-precision numerical integrators**
  - Runge-Kutta, Adams-Bashforth, and multi-step methods
  - Step-size control for accurate propagation

- **Flexible orbital dynamics models**
  - Central-body gravity
  - Third-body perturbations (optional, applied independently)
  - Modular and extensible force-model framework

- **Propagation engine**
  - Efficient batch propagation of multiple objects
  - Dynamically independent trajectories

- **Mathematical utilities**
  - Linear algebra functions
  - Vectors, matrices, and orbital state operations

- **Modular API**
  - Clear interface for integration into external programs
  - Well-defined public functions and structures

---

### Build

```bash
gcc -Iinclude src/*.c -o example1
