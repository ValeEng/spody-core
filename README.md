# Spody Core Library

**Spody Core** is a high-performance C library designed for space dynamics and astrodynamics calculations. It provides robust tools for ephemeris parsing, eclipse detection, spherical harmonics gravity modeling.

## Project Structure

The project is organized to separate the public API from the internal implementation logic:

* **`include/`**: Contains the public "umbrella" header `spody_core.h`.
* **`src/`**: Contains the source code (`.c`) and internal module headers (`.h`).
* **`CMakeLists.txt`**: The cross-platform build configuration file.

## Key Features
- **Main kernels to use**
  - Ephemeris ascii DE440 : https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/
  - Gravitational armonics data from GRGM1200A : https://pgda.gsfc.nasa.gov/products/50

## Prerequisites

To build and use this library, you will need:
- A C compiler (e.g., GCC, Clang, or MSVC).
- **CMake** version 3.10 or higher.

## Build Instructions

To compile the library on your local machine:

1. **Clone the repository**:
   ```bash
   git clone [https://github.com/ValeEng/spody-core](https://github.com/ValeEng/spody-core.git)
   cd spody_core

2. **Generate build files and compile**:
   ```bash
  mkdir build
  cd build
  cmake ..
  cmake --build .

3. **Installation**:
  If you want to create a clean package for distribution (containing only the necessary headers and the compiled library), use the following command:

    ```bash
    cmake --install . --prefix ./dist

4. **Compiling your program**:
  ```bash
  gcc main.c -I./dist/include -L./dist/lib -lspody_core -lm -o my_space_app