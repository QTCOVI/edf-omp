# FAQ

## Build and setup

### `make edf-omp` fails with "cannot find -llapack"

The `edf-omp` target links against system LAPACK/BLAS. Install them:

```bash
sudo apt install liblapack-dev libblas-dev    # Ubuntu/Debian
```

If they are installed at a non-standard path, edit `Makefile` lines 104–106 to reflect the correct paths.

### OpenMP doesn't activate — `edf-omp` still runs single-threaded

Check that you are actually running `edf-omp`, not `edf`. The static `edf` binary cannot activate OpenMP regardless of `OMP_NUM_THREADS`.

```bash
OMP_NUM_THREADS=4 ./edf-omp < input.edfinp > output.edfout
```

Verify from the timer section: if `_linearsys` wall time scales with thread count, OpenMP is working.

---

## AOM normalization errors

### `edf.f: AOM does not satisfy the sum rule`

The AOM fails the check \(\sum_A S^A_{ij} \approx \delta_{ij}\) for some orbital pair. Causes:

1. **Incomplete integration**: PROMOLDEN did not converge the atomic basin integrals — increase angular/radial grid in PROMOLDEN.
2. **Wrong AOM file**: The `.aom` file does not match the `.wfn` file (different geometry, basis, or orbital ordering).
3. **Large tolerance needed**: For fuzzy partitions, small numerical errors can accumulate. Loosen `TOLAOM` temporarily to diagnose.

**Quick fix**: add `AOMNORM n` (where `n` is the last atom) to adjust the last atom's AOM row so the sum rule holds exactly.

### All probabilities are zero

This usually means the AOM and WFN are inconsistent — the AOM integrals correspond to a different wavefunction. Verify that:
- The `.aom` file was generated from the same geometry as the `.wfn` file.
- The basis set in PROMOLDEN matches the one used in GAMESS/Gaussian.
- `ioverlap` is correct for the source of the `.aom` file.

---

## CASSCF wavefunction issues

### "The WFN file does not contain CASSCF determinants"

The `.wfn` file must come from a patched version of GAMESS10 that writes determinant coefficients and active-space occupation numbers in the extended aimpac format. Standard GAMESS10 outputs do not include this information.

### Probability sum is significantly less than 1

With `PROBCUT > 0`, only configurations above the threshold are printed. The `TOTAL SUM` line (not the `SUM` line) gives the true total. If the `TOTAL SUM` itself is much less than 1:
1. Check `TOLAOM` — the AOM normalization may be poor.
2. Check `EPSDET` — if set too aggressively, many determinant pairs are skipped.

### Calculation is very slow for large CASSCF

The cost scales as \(O(N_\text{det}^2 \times N_\alpha \times N_\beta)\). To speed up:
- Use `EPSDET 1e-8` to skip negligible determinant pairs.
- Use `EPSWFN 1e-5` to drop small-coefficient determinants before the EDF loop.
- Use `NDETS n` to limit the expansion.
- Build `edf-omp` and set `OMP_NUM_THREADS` to the number of physical cores (only effective for `ngroup = 2` or with `RECUR`).

---

## Threading and performance

### Performance degrades with more threads

Most likely cause: BLAS over-subscription. If the system BLAS (OpenBLAS, MKL) uses internal threading, each of the N OpenMP threads may spawn M BLAS threads, giving N×M total threads on a machine with N cores.

Fix:
```bash
export OPENBLAS_NUM_THREADS=1   # for OpenBLAS
export MKL_NUM_THREADS=1        # for Intel MKL
```

### No speedup for `ngroup > 2`

Correct. The parallelized `binedf.f` engine is used only for `ngroup = 2` (or when `RECUR` forces it). For more groups, the other engines (`calcedfd.f`, `rcalcedf.f`, etc.) are not yet parallelized.

---

## Output interpretation

### Printed probability sum is well below 1 but `TOTAL SUM` is 1

This is expected when `PROBCUT > 0`. The `SUM` line counts only printed configurations; the `TOTAL SUM` includes all of them. Decrease `PROBCUT` (or set `PROBCUT -1`) to see more configurations.

### Delocalization indices are negative or very large

Negative two-center DIs can appear for fragments that contain many electrons and are not physically bonded. Very large DIs between fragments that formally share a double or triple bond (e.g., N₂) are correct — \(\delta(\text{N,N}) \approx 3\) for a triple bond.

Self-DI (localization index) \(\delta(i,i)\) is always positive by construction.

### The N5+ output shows different orbital assignments than the reference

The `dsyev` LAPACK eigensolver (used since the Tier 2A upgrade) may return degenerate eigenvectors in a different order than the old `jacobi` routine. Both orderings are mathematically valid. This is an acceptable difference and does not affect EDF probabilities or delocalization indices.

---

## Input file questions

### Can I put comments in the input file?

Yes. Any line whose first non-blank character is `#` is treated as a comment and ignored.

### Does keyword order matter?

Keywords after Record 3 are mostly order-independent. However:
- `NGROUP` must appear before the fragment lines (Records 5.i).
- `MAXPOP`/`MINPOP` must appear after the `PRSRS` keyword they modify.
- `BONDING`/`ENDBONDING` must form a complete block.

### Can I omit `END`?

Yes. edf stops at end-of-file if `END` is not given.
