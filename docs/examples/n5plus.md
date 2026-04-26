# N5⁺ Cation (RHF)

N₅⁺ is a cyclic pentanitrogen cation with 34 electrons. This example demonstrates asymmetric fragment grouping, `PROBCUT` to filter low-probability configurations, and `AOMNORM` to correct the AOM normalization.

## Files

```
test/n5+.edfinp       — input
test/n5+.wfn          — RHF wavefunction
test/n5+.wfn.aom      — AOM from PROMOLDEN
test/n5+.edfout       — reference output
```

## Input

```
0
n5+.wfn.aom
n5+.wfn
probcut 0.001
aomnorm 1
ngroup 4
1 1
1 3
1 5
2 2 4
end
```

## Fragment setup

The molecule is divided into 4 groups with an asymmetric grouping: atoms 1, 3, and 5 each form their own fragment, while atoms 2 and 4 are grouped together into a single two-atom fragment.

```
Group 1: atom 1       (1 atom)
Group 2: atom 3       (1 atom)
Group 3: atom 5       (1 atom)
Group 4: atoms 2,4    (2 atoms)
```

This reflects the topology of the ring where atoms 2 and 4 play a symmetric role.

## Keywords used

### `PROBCUT 0.001`

Configurations with probability below 0.001 are not printed. For N₅⁺ this reduces the output from ~1600 lines (at `PROBCUT 0.0`) to the most chemically relevant configurations.

The sum line in the output shows what fraction of probability is captured:
```
#     0.9813550159  <-- SUM, 40 PROBABILITIES > 0.001
#     1.0000000001  <--- TOTAL SUM
```
The 40 printed configurations account for 98.1% of the total; the remaining 1.9% is distributed among very small probabilities. The `TOTAL SUM` reaching 1.0 confirms AOM normalization is correct.

### `AOMNORM 1`

Renormalizes the AOM by adjusting atom 1's elements so the sum rule is exactly satisfied. The corrected AOM is written to `n5+.wfn.aom.normalized`.

## Key output

### Average populations

```
<n(  1)_alpha>  =  3.422784     <- N atom 1
<n(  2)_alpha>  =  3.656765     <- N atom 3
<n(  3)_alpha>  =  3.131843     <- N atom 5
<n(  4)_alpha>  =  6.788608     <- N atoms 2+4 (two-atom group)
```

Group 4 (two N atoms) carries ~13.6 electrons total, consistent with 2 × 6.8 ≈ 13.6 electrons per N atom in the ring. The asymmetry between single-atom groups reflects the N₅⁺ electronic structure.

## Running

```bash
cd test/
OMP_NUM_THREADS=4 ../edf-omp < "n5+.edfinp" > n5plus_new.edfout
```

!!! note
    Shell quoting may be needed for the `n5+` filename depending on the shell.
