# Performance — edf-omp

## edf vs. edf-omp

Both binaries are compiled from the same source with the same flags
(`gfortran -fopenmp -fbounds-check -O`). The only difference is in the link step:

| | `edf` | `edf-omp` |
|-|-------|-----------|
| Link | `-static -static-libgcc` + bundled `.a` files | Dynamic system LAPACK/BLAS |
| OpenMP | Compiled in but **non-functional** | **Active** |
| Portability | Self-contained, no runtime deps | Requires system LAPACK/BLAS + libgomp |
| Threading | Always single-core | Multi-core via `OMP_NUM_THREADS` |

The Makefile comment (lines 4–5) explains the reason: with `-static -static-libgcc` the
OpenMP runtime library (`libgomp`) is linked statically in a way that prevents the
user-level `!$omp` directives from running multi-threaded. Only `edf-omp` (dynamic link)
activates threading.

## Building edf-omp

```bash
make clean && make edf-omp
```

Requires system LAPACK and BLAS at `/usr/lib/x86_64-linux-gnu/{lapack,blas}/liblapack.a`
and `libblas.a`. See [Installation](installation.md) for prerequisites.

## Setting the thread count

```bash
OMP_NUM_THREADS=8 ./edf-omp < input.edfinp > output.edfout
```

Set `OMP_NUM_THREADS` to the number of **physical cores** (not hyperthreads) for
compute-bound workloads. Hyperthreading rarely helps for floating-point-heavy loops.

To check the thread count being used, look at the timer section at the end of the output:
the `_linearsys` and `_splitedf` timers reflect wall time for the parallel regions.

## What is parallelized

The `!$omp parallel` directives in `binedf.f` (the 2-group binomial recurrence engine)
create three parallel regions:

| Region | binedf.f line | Schedule | Parallelizes |
|--------|--------------|----------|-------------|
| Alpha `_linearsys` | 649 | `dynamic,4` | Loop over alpha spin-configs (triangular; `npa` iterations) |
| Beta `_linearsys` | 804 | `dynamic,4` | Loop over beta spin-configs (triangular; `npb` iterations) |
| `_splitedf` | 938 | `dynamic,1` | Outer loop over alpha config pairs (`m1=1,npa`) |

!!! info "When does this help?"
    The `binedf.f` engine is used automatically when `ngroup = 2` (the two-group case)
    unless the `RECUR` keyword forces it also for larger systems. For `ngroup > 2`, the
    other engine files (`calcedfd.f`, `rcalcedf.f`, `xcalcedf.f`, etc.) handle the
    computation; these are **not yet parallelized**.

**Expected speedup** (2-group case):

- Near-linear for large CASSCF wavefunctions with many determinants (large `npa`, `npb`).
- Modest for small single-determinant systems where overhead dominates.
- The `_linearsys` regions use `schedule(dynamic,4)` to handle the triangular loop
  imbalance. The `_splitedf` region uses `schedule(dynamic,1)` because inner loop costs
  vary with configuration pair.

## Over-subscription warning

If the system BLAS (OpenBLAS, MKL) is compiled with internal threading, each of the N
OpenMP threads may spawn M BLAS threads, giving N×M system threads on a machine with N
cores. This severely degrades performance.

To prevent over-subscription when using `edf-omp`:

```bash
export OPENBLAS_NUM_THREADS=1   # for OpenBLAS
export MKL_NUM_THREADS=1        # for Intel MKL
export BLIS_NUM_THREADS=1       # for BLIS
```

Check which BLAS is installed:

```bash
ls -la /usr/lib/x86_64-linux-gnu/blas/libblas.a
ldconfig -p | grep -E "openblas|lapack"
```

## Practical recommendations

1. **Default**: `OMP_NUM_THREADS` = number of physical cores. Check with `nproc --all`.
2. **Large CASSCF runs** (`ndets` > 1000, `ngroup = 2`): full speedup from parallelism.
3. **Single-determinant small systems**: speedup may be modest; `edf-omp` still correct.
4. **Many-group systems** (`ngroup > 2`): no speedup from threading — use `edf` (serial
   is faster due to no OpenMP overhead).
5. **Cluster nodes**: request one full node (all cores), set `OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE`.

## Example: timed run on 1 vs. 4 threads

```bash
cd test/
time OMP_NUM_THREADS=1 ../edf-omp < cas-ch4.edfinp > /dev/null
time OMP_NUM_THREADS=4 ../edf-omp < cas-ch4.edfinp > /dev/null
```

Look at the `_linearsys` and `_splitedf` lines in the timer section to isolate the
contribution of the parallel regions.
