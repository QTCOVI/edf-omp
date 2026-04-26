# Examples

The `test/` directory contains five reference cases that cover the main usage modes of edf. All input (`.edfinp`), wavefunction (`.wfn`), AOM (`.aom`), and reference output (`.edfout`) files are included.

Run any example from the `test/` directory:

```bash
cd test/
OMP_NUM_THREADS=4 ../edf-omp < <name>.edfinp > /tmp/<name>.out
```

| Example | System | WFN type | Groups | Key features |
|---------|--------|----------|--------|--------------|
| [H6 ring](h6.md) | H₆ | RHF | 6 (one H each) | Basic QTAIM EDF |
| [N5+ cation](n5plus.md) | N₅⁺ | RHF | 4 (asymmetric) | `PROBCUT`, `AOMNORM` |
| [CH4 CASSCF](cas-ch4.md) | CH₄ | CASSCF | default (each basin) | Multi-det, spin-split |
| [N2H2 CASSCF](cas-n2h2.md) | N₂H₂ | CASSCF | 3 | `PROBCUT -1`, `AOMNORM` |
| [N2 CASSCF+rhfmode](n2casb.md) | N₂ | CASSCF | default | Natural orbitals, RHF-like path |

---

## Regression testing

There is no automated test runner. To validate a change, diff the new output against the committed reference:

```bash
cd test/
OMP_NUM_THREADS=1 ../edf-omp < h6.edfinp > /tmp/h6_new.out
diff <(grep -v 'Calculation starts\|Calculation ends' h6.edfout) \
     <(grep -v 'Calculation starts\|Calculation ends' /tmp/h6_new.out)
```

Acceptable differences:
- Timestamp lines (`Calculation starts/ends`).
- Last-digit (1-ULP) floating-point variations from different LAPACK routine ordering.
- Different ordering of degenerate eigenvectors in EOS output (both valid).
