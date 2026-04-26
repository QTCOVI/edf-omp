# edf — Electron Number Distribution Functions

**edf** computes the probability distribution \(p(n_1, n_2, \ldots, n_m)\) of finding exactly
\(n_1\) electrons in fragment 1, \(n_2\) electrons in fragment 2, …, and \(n_m\) electrons in
fragment \(m\) of a molecule. It reads an aimpac `.wfn` wavefunction file (from GAMESS or
Gaussian) together with an Atomic Overlap Matrix (AOM) and produces a complete electron-number
distribution for any user-defined partitioning of molecular space.

Developed by **E. Francisco & A. Martín Pendás** (University of Oviedo), © 2022.

## What edf computes

Given a wavefunction and a set of \(m\) non-overlapping, space-filling fragments:

- **The full EDF**: \(p(n_1, \ldots, n_m)\) — probability of each electron configuration.
- **Spin-resolved EDF**: separate alpha and beta electron distributions.
- **Spinless EDF**: sum over spin to get total-electron distributions.
- **Average populations** \(\langle n_i \rangle\) and **delocalization indices** \(\delta(i,j)\).
- **Covariance matrices** and mutual information entropy (optional).
- **Localized MOs** via Cioslowski isopycnic localization (optional).
- **DAFH**, **EOS/Lewis structures**, **real-space resonance structures** (optional analyses).

## Basin types supported

| Type | Source program |
|------|---------------|
| QTAIM atomic basins | PROMOLDEN, PROAIM, TOPMOD, AIMALL |
| ELF basins | TOPMOD |
| Fuzzy atomic basins (Mulliken, Löwdin, MinDef, Becke, promolecular) | Computed internally |

## Wavefunction types

- **Single-determinant**: RHF, ROHF, UHF (closed-shell fast path via `ioverlap = -11`).
- **Multi-determinant**: CASSCF. Determinant coefficients are read from the `.wfn` file
  (GAMESS format).

## Two binaries

| Binary | Build command | OpenMP | BLAS |
|--------|--------------|--------|------|
| `edf` | `make` | Compiled in, **not active** (static link) | Bundled static |
| `edf-omp` | `make edf-omp` | **Active** | System dynamic |

See [Performance (edf-omp)](performance.md) for details on when and how to use `edf-omp`.

## Basic usage

```
./edf-omp < input.edfinp > output.edfout
```

Input is read from **stdin**; output goes to **stdout**. See [Quick Start](quickstart.md)
for a complete worked example.
