# Input Format

`edf` reads its directives from **stdin**. The first three records are **positional** and
mandatory; everything after them is keyword-driven.

```
./edf-omp < input.edfinp > output.edfout
```

Comments: any line whose first non-blank character is `#` is ignored.

Terminator: `END` (case-insensitive) ends the input.

---

## Record 1 — `ioverlap`

An integer on its own line that selects the source and format of the Atomic Overlap Matrix.

| Value | AOM source |
|-------|-----------|
| `0` | PROMOLDEN format **or** `filedat` is a fuzzy-partition keyword (see Record 2) |
| `1` | PROAIM format |
| `2` | TOPMOD (ELF basins) |
| `3` | TOPMOD (QTAIM basins) |
| `4` | AIMALL (QTAIM basins) |
| `-1` | Complex AOM — control passed to `edfx.f` |
| `-2` | Complex AOM — control passed to `edfcrit.f` |
| `-11` | Single-determinant fast path — control passed to `edfstd.f` |

See [AOM Formats](aom-formats.md) for details on each format.

---

## Record 2 — `filedat`

Either the **name of the AOM file** or one of the built-in fuzzy-partition keywords.
Requires `ioverlap = 0` when a keyword is used.

| Value | Meaning |
|-------|---------|
| *filename* | Read AOM from this file |
| `mulliken` | Compute AOM using Mulliken population analysis |
| `lowdin` | Compute AOM using Löwdin population analysis |
| `mindef` | Compute AOM using Minimally Deformed Atoms (MinDef) of Fernández Rico et al. |
| `mindefrho` | MinDef with weights \(w_A = \rho_A/\rho\) (always gives \(S_{ii}^A > 0\)) |
| `netrho` | Net-density weights (eliminates mixed-center cross terms) |
| `promrho` | Promolecular density weights: \(w_A = \rho_A^{\rm atom}/\rho\) |
| `heselmann` | Heselmann chemical-localization weights (JCTC 12, 2720, 2016) |
| `becke` | Becke fuzzy-cell partition |

For `mindefrho`, `netrho`, and `becke`, optional integration parameters (`NANG`, `NRAD`, `IRMESH`,
and for Becke also `POW`) may follow the keyword on the same line. See
[Integration Keywords](keywords/integration.md).

---

## Record 3 — `wfnfile`

The name of the `.wfn` wavefunction file. Must be in the GAMESS aimpac format.

- **Single-determinant**: standard `.wfn` from GAMESS or Gaussian.
- **Multi-determinant (CASSCF)**: `.wfn` with determinant coefficients and active-space
  occupation numbers, written by a domestic version of GAMESS10.

The AOM integrals in `filedat` must correspond to the wavefunction in `wfnfile`.

---

## Record 4 — `NGROUP ngroup`

The number of fragments into which the molecule is divided.

- If `ngroup ≤ 0`, or `NGROUP` is absent, each atomic basin becomes its own independent
  fragment with full electron-count range (0 to total alpha/beta). In this case Records 5.i
  must be **skipped**.
- Negative `ngroup`: the first `|ngroup|−1` groups are specified in Records 5.i, and the
  last group automatically contains all remaining atoms.

A variant keyword `NFRAG nfrag` sets the fragment count used exclusively by the EOS/Lewis
routines (`nfnoeos.f`, `lewisoqs.f`). If absent, `NFRAG` defaults to `NGROUP`.

---

## Records 5.i — Fragment specification (i = 1 … NGROUP)

Two formats are accepted:

**Format 1** (with electron count limits):
```
nfugrp  atom1  atom2  ...  minelecA  maxelecA  minelecB  maxelecB
```

**Format 2** (full range assumed):
```
nfugrp  atom1  atom2  ...
```

| Field | Meaning |
|-------|---------|
| `nfugrp` | Number of atoms in this group |
| `atom1 atom2 …` | Indices of the atoms belonging to this group (1-indexed) |
| `minelecA`, `maxelecA` | Min/max alpha electrons in this group |
| `minelecB`, `maxelecB` | Min/max beta electrons in this group |

When Format 2 is used, the program sets `minelec = 0` and `maxelec = NALPHA` (or `NBETA`).

!!! warning "Core electrons"
    The `minelec`/`maxelec` values **include** core electrons (as defined by `COREMO`).
    They are total electron counts, not valence-only counts.

**Example** — H₆ ring, 6 groups of 1 atom each:
```
ngroup 6
1 1
1 2
1 3
1 4
1 5
1 6
```

**Example** — N₅⁺, 4 groups with atoms 2 and 4 together:
```
ngroup 4
1 1
1 3
1 5
2 2 4
```

---

## Records 5bis.i — EOS/Lewis fragment specification

These are used only when `NFRAG` differs from `NGROUP` and EOS/Lewis analysis is requested.
Format: `nwithin atom1 atom2 …` (similar to Records 5.i, Format 2, but for `NFRAG` fragments).

---

After the fragment specification, all remaining input is **keyword-driven** and
optional (except `END`). See the [Keywords](keywords/index.md) section.
