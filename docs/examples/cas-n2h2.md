# N₂H₂ CASSCF

Diazene (N₂H₂) with a CASSCF wavefunction (208 determinants, 16 electrons, 12 active). This example shows `PROBCUT -1` to print all computed probabilities including near-zero ones, combined with `AOMNORM` for AOM correction.

## Files

```
test/cas-n2h2.edfinp      — input
test/cas-n2h2.wfn         — CASSCF wavefunction
test/cas-n2h2.wfn.aom     — AOM from PROMOLDEN
test/cas-n2h2.edfout      — reference output
```

## Input

```
0
cas-n2h2.wfn.aom
cas-n2h2.wfn
probcut -1.0
aomnorm 1
ngroup 3
1 1
1 2
2 3 4
end
```

## Fragment setup

Three groups:
- **Group 1**: atom 1 (one nitrogen)
- **Group 2**: atom 2 (other nitrogen)
- **Group 3**: atoms 3 and 4 (both hydrogens grouped together)

This choice probes the N–N bond and the two N–H bonds collectively.

## Keywords used

### `PROBCUT -1.0`

A negative `PROBCUT` prints **all probabilities**, even those that are formally zero or tiny due to round-off. The output will include the complete 2025-entry spin-split probability table:

```
#     0.9999999999902743  <- SUM (2025 PROBS with P_{a} & P_{b} > -1.0E+00)
#     0.9999999999902743  <- TOTAL SUM
```

The sum reaches 0.99999999 (essentially 1.0), confirming that the AOM + wavefunction are consistent. All 45 spinless alpha and 45 spinless beta probabilities are printed.

### `AOMNORM 1`

Renormalize the AOM using atom 1 as the reference.

## Key output

### Wavefunction info

```
# Number of electrons               =   16
# Number of determinants in WFN     =  208
# Number of active electrons in WFN =   12
```

With 208 determinants and only 3 groups, the calculation is fast. The active-space electrons are the 12 electrons in the CASSCF active space; 4 are core.

### Average populations

```
<n(  1)_alpha>  =  3.687780     <- N atom 1
<n(  2)_alpha>  =  3.687653     <- N atom 2
<n(  3)_alpha>  =  0.624567     <- both H atoms
```

The two nitrogen atoms carry similar populations (~3.69 alpha electrons each ≈ 7.38 total), and the two hydrogens together carry ~1.25 alpha electrons ≈ 2.5 total. The nearly equal N populations confirm the trans-N₂H₂ symmetry.

### Delocalization indices

The N–N delocalization index (between groups 1 and 2) and the N–H index (between group 3 and each N) characterize the bonding in the CASSCF picture. With a multi-determinant wavefunction, the N–N DI reflects the bond order including correlation effects.

## Running

```bash
cd test/
OMP_NUM_THREADS=4 ../edf-omp < cas-n2h2.edfinp > cas-n2h2_new.edfout
```

With 208 determinants, this runs in under a second on any modern machine.
