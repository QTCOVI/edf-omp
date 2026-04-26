# Basic Control Keywords

These keywords control thresholds, truncation of the wavefunction, and the core computation path.

---

## `AOMNORM [atom]`

Renormalizes the Atomic Overlap Matrix so that the sum rule

$$\sum_A S^A_{ij} = \delta_{ij}$$

is satisfied exactly. The correction is applied by adjusting the AOM elements of atom `atom` (1-indexed). If `atom` is omitted, the last atom in the WFN file is used.

After renormalization, a file `<filedat>.normalized` is written containing the corrected AOM.

!!! tip
    Use `AOMNORM` when `TOLAOM` failures occur on molecules where numerical integration left one atom's AOM slightly off.

---

## `PROBCUT real`

Probabilities smaller than `real` are not written to the output. Default is `0.0` (all computed probabilities are printed).

Set `PROBCUT -1.0` to guarantee that every non-zero probability appears, even if round-off makes some of them formally zero.

**Example**:
```
probcut 0.001
```
Prints only configurations with probability > 0.001.

---

## `EPSDET real`

In a CASSCF calculation, skip the computation of the contribution from a determinant pair \((i, j)\) when \(|C_i \cdot C_j| < \text{real}\). Default is `0.0` (all pairs are included).

Setting `EPSDET` to a small positive number (e.g. `1e-8`) can greatly speed up large CASSCF runs at the cost of a small loss of accuracy.

---

## `EPSPROBA real`

Skip the beta-electron probability loop for an alpha configuration when the highest alpha probability for that configuration is smaller than `real`. Avoids computing beta probabilities for negligible alpha states. Default is `0.0`.

---

## `TOLAOM real`

Tolerance for the AOM sum-rule check. If

$$\left|\sum_A S^A_{ij} - \delta_{ij}\right| > \text{real}$$

for any \((i,j)\), edf stops with an error. Default is `0.01`.

Tighten this to `1e-4` for high-accuracy work; loosen it (e.g. `0.05`) only if the AOM source is known to be slightly imprecise.

---

## `NDETS n`

In a CASSCF calculation, use only the first `n` determinants from the WFN file. The wavefunction is renormalized after truncation. If absent, all determinants are used.

Useful for testing convergence with respect to the number of determinants or for speeding up exploratory runs.

---

## `EPSWFN real`

Drop all determinants whose coefficient \(|C_i| < \text{real}\) and renormalize. Applied before `NDETS`. Default is `0.0` (all determinants kept).

---

## `RECUR`

Force the binomial recurrence path (`binedf.f`, Cançes et al.) even when there are more than 2 groups. By default, the recurrence is used only for `ngroup = 2`.

!!! note
    The recurrence path is the one parallelized in `edf-omp`. For `ngroup > 2`, forcing `RECUR` enables threading but may not be faster than the direct path.

---

## `NOMEM`

By default, for multi-determinant wavefunctions, intermediate alpha and beta probability arrays are kept in memory (`memedf = .true.`, handled by the in-memory engine). `NOMEM` switches to disk-based storage.

Use `NOMEM` on machines with limited RAM and large active spaces. The calculation becomes I/O-bound instead of memory-bound.

---

## `RANDOM minrandom maxrandom`

Sets the minimum and maximum values of the random number used to perturb the linear system in `binedf.f`. This is only relevant for the two-group recurrence path.

The default values are chosen to give a well-conditioned system; this keyword is rarely needed.
