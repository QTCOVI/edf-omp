# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`edf` — Electron Number Distribution Functions. A Fortran 77/90 quantum-chemistry program by E. Francisco & A. Martín Pendás (University of Oviedo) that reads an aimpac `.wfn` wavefunction (from GAMESS/GAUSSIAN) plus an Atomic Overlap Matrix (AOM), and computes the probability distribution of finding exactly (n_1, n_2, …, n_m) electrons in m user-defined molecular fragments. Supports QTAIM, ELF, and fuzzy (Mulliken/Löwdin/MinDef/Becke) basin partitions; single-determinant and multi-determinant (CASSCF) wavefunctions; optional Cioslowski isopycnic localization; DAFH, EOS, and Lewis-structure analyses.

The full input-file grammar is documented as comments at the top of `edf.f` (records 1–N, starting around line 89). That file is authoritative — always read it before reasoning about input keywords.

## Build

```
make              # build ./edf binary (default, static)
make clean        # remove objects, *.mod, and the edf binary
```

Compiler settings live in the top of `Makefile`. Two configurations are toggleable by commenting/uncommenting blocks:
- **Static** (default, lines 7–15): `-static`, links bundled `libblas.a` / `liblapack.a` shipped in the repo root. OpenMP is compiled in but effectively unused — the note in the Makefile says "openmp CAN NOT be used" in the static build.
- **Non-static** (lines 21–28): uses system LAPACK/BLAS at `/usr/lib/x86_64-linux-gnu/{lapack,blas}/`. Use this if you want OpenMP to actually run.

Flags default to `-O -fopenmp -fbounds-check`. A stricter warning line is commented out on Makefile:9 — enable it when hunting uninitialized-variable bugs.

No autoconf/cmake. The `Makefile` hand-lists every object in `OBJECTS` and every include-file dependency after it. **When adding a new `.f`/`.f90` source file, you must append it to `OBJECTS` and add its `*.inc` dependencies** — there is no auto-dependency generation.

## Running

`edf` reads its directives from **stdin** and writes to **stdout**. Typical usage:

```
./edf < test/h6.edfinp > test/h6.edfout
```

By convention, input files use the `.edfinp` suffix and outputs use `.edfout`. The `test/` directory holds reference inputs (`h6`, `n5+`, `cas-ch4`, `cas-n2h2`, `n2casb`) alongside their `.wfn`, `.aom`, and expected `.edfout` files — use these as regression fixtures. There is no automated test runner; validate changes by diffing new output against the committed `.edfout`.

The first three input records are positional (`ioverlap`, AOM filename, WFN filename); everything after is keyword-driven (`NGROUP`, `AOMNORM`, `PROBCUT`, `FULLOC`, `COREMO`, `DELETE`, `SUPPRESS`, `PRSRS`, `DAFH`, `OQSEOS`, `LEWISOQS`, …). See the `edf.f` header for the full list.

## Architecture

### Entry point and driver
- `edf.f` — `program edf`. Large (~3700 lines) monolithic driver: parses stdin, reads WFN + AOM, then dispatches to the EDF engines based on `ioverlap` and the keyword flags it accumulated. Standard units: `stdin=5`, `stdout=6`. When `ioverlap` has special values, control is handed off whole to `edfx` (complex AOM), `edfcrit` (complex AOM from critic2), or `edfstd` (single-determinant fast path).

### EDF computation engines
Multiple specialized engines compute the EDF under different assumptions — they share structure but are kept as separate files:
- `calcedf.f` / `mcalcedf.f` — generic multi-determinant EDF (in-file vs. in-memory storage; see `NOMEM` keyword).
- `rcalcedf.f` — real-AOM variant. `xcalcedf.f` / `xcalcedfd.f` — complex-AOM variants.
- `calcedfd.f` — direct-access-file variant. `edfstd.f` — single-determinant closed-shell fast path. `edfx.f` / `edfcrit.f` — complex-AOM drivers.
- `binedf.f` — binomial-recurrence path (Cançes et al.) used for 2-group cases unless `RECUR` forces it.

### Data partitioning (AOM construction)
When `filedat` is a keyword rather than a filename, the AOM is computed in-process:
- `aomulliken.f`, `aomlowdin.f`, `aomindef.f`, `aomindefrho.f`, `aombecke.f`, `netrhow.f`, `promrhow.f` — one file per fuzzy-atom weighting scheme.
- `readaom.f`, `readaomstd.f`, `readcmplx.f` — parsers for PROMOLDEN / PROAIM / TOPMOD / AIMALL / complex AOM formats.

### Localization and analysis
- Isopycnic (Cioslowski) localization: `lmowfn.f`, `itaomloc.f`, `uhfitaomloc.f`, `locsome.f`, `coremos.f`, `delmo.f`.
- DAFH (Domain-Averaged Fermi Hole): `dafh.f`, `dafhdo.f`, `dafhdrv.f`, `dafhstd.f`, `dafhcrit.f`, `gendafh.f`, `itdafh.f`, `itdafhc.f`.
- EOS / Lewis structures: `neweos.f`, `sdweos.f`, `fnoeos.f`, `nfnoeos.f`, `lewisoqs.f`, `ulewisoqs.f`, `sdwlewis.f`, `oqsloc.f`.
- Probabilities of real-space resonance structures (`PRSRS`): `prsrs.f`, `wrsrs.f`, `sdwrsrs.f`, `probres.f`.
- Symmetry: `sym.f`, `newsym.f`, `pntgrp.f`, `symgr.f`, `typeop.f`, `mod_sym.f90` (+ `sym.inc`).

### State: modules and include files
State is shared two ways — both are pervasive; grep first, then edit.

1. **Fortran modules (`space_for_*`)** — one per file, each a thin wrapper around allocatable arrays that are filled once and read everywhere. Used for larger/dynamic data:
   - `wfnbasis.f` (`space_for_wfnbasis`), `wfncoef.f` (`space_for_wfncoef`), `primgto.f` (`space_for_primgto`), `cidetarea.f` (`space_for_cidet`), `rdm1space.f` (`space_for_rdm1`), `rsspace.f` (`space_for_rsrs`), `bondarea.f` (`space_for_bonds`), `dafhspace.f` (`space_for_dafh`), `sgarea.f` (`space_for_sgarea`), `spaceconf.f` (`space_for_conf`), `symarea.f` (`space_for_sym`), `aomspin.f` (`space_for_aomspin`), `wfxedf.f` (`space_for_wfxedf`).
   - Larger F90 modules: `mod_sym.f90`, `mod_periodic.f90`, `mod_param.f90`.

2. **Classic COMMON blocks via `.inc` files** — `INCLUDE`-d at the top of each subroutine that touches them. The key ones:
   - `implicit.inc` → `implicit real(kind=8) (a-h,o-z)` (project-wide implicit typing — assume it everywhere).
   - `param.inc` → `maxtype`, `nstack`, `maxndet`, `pi`, and the `/logic/` common (`icorr, rhf, rohf, uhf, cciqa, onlynat, thereisEDF, woccup`) that flags the wavefunction type.
   - `wfn.inc` → `/iatdat/ nprims, nmo, ncent` (wavefunction dimensions).
   - `corr.inc` → `/icor/ ndets, ncore, nact, nelact, nel, nalpha, nbeta, mocore, moval` (electron / orbital counts).
   - `constants.inc`, `lengrec.inc`, `fact.inc`, `integ.inc`, `mline.inc`, `stderr.inc`, `error.inc`, `opt.inc`, `sym.inc`, `datatm.inc`, `primgto.inc`, `maxrdm.inc`.

If you change a struct in a `.inc` file, every `.f` that `INCLUDE`s it must be rebuilt — `make clean && make` is the safe move.

### I/O helpers and utilities
- WFN reading/writing: `rdwfn.f`, `wrtwfn.f`, `wfxedf.f`, `cutwfn.f`, `wfnmax.f`.
- Sorting: `qcksort.f`, `iqcksort.f`, `qqsort.f`, `chsort.f`.
- Linear algebra: bundled `eispack.f`, `jacobi.f`, `svdcmp.f`, `ludcmp.f` / `lubksb.f`, plus `mylapack.f` / `lulapack.f90` / `detlapack.f90` wrappers over the shipped LAPACK/BLAS `.a` files.
- Lebedev angular grids: `sphere_lebedev_rule.f` (single 237k-line file — don't try to refactor casually).
- Misc parsers: `atoi.f`, `isdigit.f`, `uppcase.f`, `lower.f`, `leng.f`, `setint.f`, `setdble.f`, `setword.f`, `searchstring.f90`.

## Working notes

- Fortran fixed-form (`.f`) vs. free-form (`.f90`) is respected by the Makefile via suffix rules — don't rename files across the boundary without updating the source.
- `implicit real(kind=8)` is active project-wide via `implicit.inc`. Variables starting `i–n` are integer, everything else is `real*8`, unless explicitly typed otherwise.
- The `.o` files checked into the repo root (~100 of them) are stale build artifacts. Prefer `make clean` before benchmarking or diffing.
- Archive `edf-06-04-2026.tar` and `test/edf-examples.tar` are snapshots; ignore unless the user specifically references them.
- Shell is `tcsh` on this machine — redirection syntax differs slightly from bash (`>&` combines stdout+stderr; no `&>`).
