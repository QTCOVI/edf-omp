# Installation

## Prerequisites

- `gfortran` ≥ 9 (supports Fortran 2008 features and `-fopenmp`)
- `make`
- For `edf` (static): no extra runtime libraries needed — bundled `liblapack.a`/`libblas.a` are included in the repo.
- For `edf-omp` (dynamic, **recommended**): system LAPACK and BLAS at `/usr/lib/x86_64-linux-gnu/lapack/` and `/usr/lib/x86_64-linux-gnu/blas/`.

On Ubuntu/Debian:

```bash
sudo apt install gfortran liblapack-dev libblas-dev
```

## Building `edf` (serial, static)

```bash
make clean
make
```

The resulting binary is `./edf`. It is compiled with `-fopenmp` but the static linking
(`-static -static-libgcc`) prevents the OpenMP runtime from loading dynamically, so all
calculations run on a single core regardless of `OMP_NUM_THREADS`.

!!! note
    The `.o` files checked into the repository are stale build artifacts. Always run
    `make clean` before a fresh build or before benchmarking.

## Building `edf-omp` (OpenMP-enabled, dynamic)

```bash
make clean
make edf-omp
```

The resulting binary is `./edf-omp`. It uses **dynamic** LAPACK/BLAS and a dynamically
loaded OpenMP runtime (`libgomp`), so `!$omp parallel` directives in the source actually
run multi-threaded.

!!! warning "System LAPACK path"
    The `edf-omp` target hard-codes the library path to
    `/usr/lib/x86_64-linux-gnu/lapack/liblapack.a` and
    `/usr/lib/x86_64-linux-gnu/blas/libblas.a`.
    On other distributions edit `Makefile` lines 104–106 accordingly,
    or replace the `.a` paths with `-llapack -lblas` if dynamic `.so` files are preferred.

## Cleaning

```bash
make clean
```

Removes all `.o` files, `.mod` files, the `edf` binary, and the `edf-omp` binary.

## Compiler flags

The default flags (Makefile line 10) are:

```
gfortran -fopenmp -fbounds-check -O
```

`-fbounds-check` enables array-bounds checking at a small runtime cost. Remove it for
production runs on large systems by editing the `FCOMPL` line. A stricter warning set
is commented out on Makefile line 9:

```
gfortran -fopenmp -Wall -Wunused-label -Wmaybe-uninitialized -fbounds-check
```

Enable this when hunting uninitialized-variable bugs.
