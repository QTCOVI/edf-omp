# edf / edf-omp — Electron Number Distribution Functions

**edf** computes the probability distribution *p*(n₁, n₂, …, nₘ) of finding exactly n₁ electrons
in fragment 1, n₂ in fragment 2, … and nₘ in fragment m of a molecule. It reads an aimpac `.wfn`
wavefunction (from GAMESS or Gaussian) together with an Atomic Overlap Matrix (AOM) and produces
the complete electron-number distribution for any user-defined partitioning of molecular space.

Developed by **E. Francisco & A. Martín Pendás** (University of Oviedo).

## What it computes

- **Full EDF** — probability of every electron configuration (n₁, …, nₘ).
- **Spin-resolved and spinless EDFs** — alpha/beta and total-electron distributions.
- **Average populations** ⟨nᵢ⟩ and **delocalization indices** δ(i,j).
- **Covariance matrices** and mutual-information entropy.
- **Localized MOs** via Cioslowski isopycnic localization.
- **DAFH**, **EOS/Lewis structures**, and **real-space resonance structures**.

## Basin types

| Type | Source |
|------|--------|
| QTAIM atomic basins | PROMOLDEN, PROAIM, TOPMOD, AIMALL |
| ELF basins | TOPMOD |
| Fuzzy basins (Mulliken, Löwdin, MinDef, Becke) | Computed internally |

## Wavefunction types

Single-determinant (RHF, ROHF, UHF) and multi-determinant CASSCF.

## Two binaries

| Binary | Build | OpenMP |
|--------|-------|--------|
| `edf` | `make` | compiled in, not active (static link) |
| `edf-omp` | `make edf-omp` | active — use for parallel runs |

## Quick start

```
make edf-omp
./edf-omp < test/h6.edfinp
```

## Documentation

Full user manual, keyword reference, worked examples, and theory:
**https://qtcovi.github.io/edf-omp/**
