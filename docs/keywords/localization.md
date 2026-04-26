# Localization Keywords

These keywords request an isopycnic (Cioslowski) localization of molecular orbitals. All localization methods maximize a localization criterion constructed from the AOM and produce a set of MOs that are normalized but **not** necessarily orthogonal.

---

## Localization strategies

### `FULLOC [ALLCENTER]`

Localize the **natural orbitals** from the input WFN using isopycnic localization. The localization function is built from the fragments defined by `NGROUP`.

If `ALLCENTER` is given, the localization function is built using every atomic center as an individual fragment (instead of the user-defined groups).

### `CANLOC [ALLCENTER]`

Localize the **canonical orbitals** from the input WFN. Otherwise identical to `FULLOC`.

### `CHEMLOC`

Chemical localization variant. Uses the same isopycnic algorithm but with different weighting — the criterion is more directly linked to chemical bond regions.

### `AOMLOC`

A second chemical localization variant, closer to the AOM-based formulation.

### `DAFHLOC`

Localization via progressive diagonalization of the **Domain-Averaged Fermi Hole (DAFH)** in 1-, 2-, 3-, … center domains up to `MXCENT` centers. Produces MOs that are localized on a growing number of atomic centers, yielding 1c-2e (lone pairs), 2c-2e (bonds), and higher-center contributions.

### `OPENLOC`

Localization by diagonalizing the **Open Quantum System (OQS)** 1-RDM, i.e., \(S_A^{1/2} \cdot \text{1-RDM} \cdot S_A^{1/2}\), for each fragment \(A\). Produces Fragment Natural Orbitals (FNOs) localized on each domain.

---

## Block keywords

### `DODAFH … END`

Provides full manual control over the DAFH localization sequence. Within the `DODAFH` … `END` block, each line specifies one localization step:

```
ncen  cutoff  atom1  [atom2 …]  [DEPLET]
```

| Field | Meaning |
|-------|---------|
| `ncen` | Number of centers to localize on |
| `cutoff` | Acceptance threshold (close to but below 1.0) |
| `atom1 …` | Indices of the `ncen` atoms |
| `DEPLET` | After this step, deplete (remove found MO from) the full DAFH |

**Example** — C₂H₆, manual DAFH sequence:
```
dodafh
  1 0.95 1
  1 0.95 2   deplet
  2 0.95 1 2 deplet
  2 0.95 1 5
  2 0.95 1 6
  2 0.95 1 7 deplet
  2 0.95 2 3
  2 0.95 2 4
  2 0.95 2 8 deplet
end
```

### `OTHEROPEN … END`

Equivalent to `DODAFH … END` but uses OQS 1-RDM diagonalization instead of DAFH. Same syntax.

---

## Localization control parameters

### `DAFHSELECT n n1 n2 … nn`

Restrict the DAFH progressive localization to specific atom tuples. For `n`-center steps, only the tuples `(n1, n2, …, nn)` are explored instead of all combinations.

Multiple `DAFHSELECT` lines with the same `n` accumulate — all listed tuples are explored. If no `DAFHSELECT` line is given for a value of `n`, all n-tuples are explored.

Currently only works for single-determinant wavefunctions in `DAFHLOC` mode.

### `CRITICS real`

When `CANLOC` or `CHEMLOC` is active, an MO is considered to have partial localization on a fragment if its self-overlap \(S^A_{ii} > \text{real}\) in that fragment. Default is `0.02`.

### `COVX real`

Bond detection scaling factor. Atoms A and B are considered bonded if their distance is less than \((R_A + R_B) \times \text{COVX}\) where \(R_A\) and \(R_B\) are stored covalent radii. Default is `1.2`.

### `DAMPS real`

Damping factor for `CHEMLOC`. In each localization cycle, all `CRITOV` values are multiplied by `DAMPS` to progressively relax the localization criterion. Must be less than 1.0. Default is `0.95`.

### `MXCENT n [skip1 skip2 …]`

Maximum number of centers in the progressive DAFH localization (1c, 2c, …, up to `n`c). Default is `6`.

Optional integers after `n` specify center counts to skip. For example:

```
MXCENT 6 3 4
```

analyzes 1-, 2-, 5-, and 6-center localizations, skipping 3 and 4.

### `CRITOV [n critval]`

Overlap threshold for accepting an MO as localized on an n-center domain. Default for all `CRITOV(1..MXCENT)` is `0.95`.

- `CRITOV` alone — no change.
- `CRITOV n critval` — set `CRITOV(n) = critval`.
- `CRITOV critval` (real number, no integer first) — set all `CRITOV(1..MXCENT) = critval`.

Shorthand forms: `CRITOV1 critval` sets `CRITOV(1)`, `CRITOV2 critval` sets `CRITOV(2)`.

### `MAXBOND maxbond`

In the progressive DAFH partitioning for atom pairs, the maximum number of bonds separating atom1 from atom2 that is actually explored. Default is `1` (only directly bonded pairs). Maximum is `4`.

### `SKIPH`

When true, hydrogen atoms are skipped in the single-atom DAFH analysis. Default is false.

### `EPSNEG [epsneg]`

Largest allowed (in absolute value) negative eigenvalue in `OPENLOC`/`OTHEROPEN` diagonalizations of matrices that should be positive definite. Default is `-1e-6`.

### `EPSEIGEN [epseigen]`

Eigenvalues smaller than `epseigen` in the OQS 1-RDM diagonalization are considered to correspond to MOs localized outside the fragment and are discarded.

### `SMALLEIGEN [smalleigen]`

Minimum eigenvalue floor when computing \(S_A^{-1/2}\). Eigenvalues of \(S_A\) below this value are replaced by `smalleigen` to avoid division by zero.
