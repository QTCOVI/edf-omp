# AOM Formats

The **Atomic Overlap Matrix (AOM)** stores the integrals \(S^A_{ij} = \langle \phi_i | \hat{w}_A | \phi_j \rangle\) of each molecular orbital pair \(i,j\) over the basin of atom \(A\). These integrals encode the partitioning of space into molecular fragments.

The source and format of the AOM is selected by `ioverlap` (Record 1).

---

## `ioverlap = 0` — PROMOLDEN format or fuzzy partition

When `filedat` is a filename, edf reads an AOM file produced by **PROMOLDEN** (the standard case for QTAIM basins). The PROMOLDEN format is a plain-text table with one block per atom.

When `filedat` is a keyword (`mulliken`, `lowdin`, `mindef`, etc.), edf computes the AOM internally using a fuzzy-atom weighting scheme. See [Input Format — Record 2](input-format.md#record-2-filedat) for the full list of keywords.

---

## `ioverlap = 1` — PROAIM format

AOM from the **PROAIM** program (QTAIM basins). PROAIM is an older topological analysis code; its AOM format differs from PROMOLDEN in the header and ordering of matrix elements. edf detects and parses this format automatically when `ioverlap = 1`.

---

## `ioverlap = 2` — TOPMOD, ELF basins

AOM from **TOPMOD** using **ELF** (Electron Localization Function) basins. ELF basins divide space according to the topology of the electron localization function rather than the electron density.

---

## `ioverlap = 3` — TOPMOD, QTAIM basins

AOM from **TOPMOD** using **QTAIM** (Quantum Theory of Atoms in Molecules) basins, i.e., the gradient-field basins of the electron density \(\rho\). Same program as `ioverlap = 2` but with different basin topology.

---

## `ioverlap = 4` — AIMALL format

AOM from **AIMALL** (QTAIM basins). AIMALL produces a dedicated AOM output file; select this value when the `.aom` file was generated with AIMALL rather than PROMOLDEN or PROAIM.

---

## `ioverlap = -1` — Complex AOM (edfx)

A **complex** Atomic Overlap Matrix. Control passes immediately to `edfx.f`, which has its own input grammar. Consult the `edfx.f` source file for the format of additional input records in this mode.

---

## `ioverlap = -2` — Complex AOM (edfcrit)

A **complex** Atomic Overlap Matrix sourced from **critic2**. Control passes to `edfcrit.f`. Consult `edfcrit.f` for the additional input expected in this mode.

---

## `ioverlap = -11` — Single-determinant fast path

For **single-determinant** wavefunctions (RHF, ROHF, UHF). After reading the AOM, control passes to `edfstd.f`, which uses an optimized closed-shell algorithm. Additional keywords specific to `edfstd.f` may follow. Consult `edfstd.f` for its input grammar.

---

## AOM normalization check

For all formats, edf verifies that the AOM satisfies the idempotency sum rule:

$$\sum_A S^A_{ij} = \delta_{ij}$$

If any diagonal deviation exceeds `TOLAOM` (default 0.01), edf stops with an error. Use `AOMNORM` to fix a poorly normalized AOM (see [Keywords — Basic](keywords/basic.md#aomnorm-atom)).

---

## Summary table

| `ioverlap` | Source | Basin type |
|-----------|--------|-----------|
| `0` | PROMOLDEN file **or** fuzzy keyword | QTAIM or fuzzy |
| `1` | PROAIM | QTAIM |
| `2` | TOPMOD | ELF |
| `3` | TOPMOD | QTAIM |
| `4` | AIMALL | QTAIM |
| `-1` | Complex AOM → `edfx.f` | — |
| `-2` | Complex AOM from critic2 → `edfcrit.f` | — |
| `-11` | Any AOM → `edfstd.f` fast path | Single-det only |
