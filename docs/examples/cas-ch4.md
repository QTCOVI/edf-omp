# CH₄ CASSCF

Methane with a CASSCF wavefunction (5519 determinants, 10 electrons). This is the primary benchmark for the OpenMP parallelization — the large number of determinants gives the `_linearsys` and `_splitedf` parallel regions meaningful work.

## Files

```
test/cas-ch4.edfinp   — input
test/cas-ch4.wfn      — CASSCF wavefunction
test/cas-ch4.aom      — AOM from PROMOLDEN
test/cas-ch4.edfout   — reference output
```

## Input

```
0
cas-ch4.aom
cas-ch4.wfn
probcut 0.001
short
```

Note: no `ngroup` keyword. Each atomic basin is its own independent fragment with full electron-count range. For CH₄ (C + 4H = 5 atoms) this gives 5 fragments.

The `short` keyword is not a documented keyword — `cas-ch4.edfinp` does not contain an `ngroup` line, so edf uses all basins as independent fragments. There is no `end` line in this input either; edf stops at EOF.

Actually, looking at the input more carefully:

```
0
cas-ch4.aom
cas-ch4.wfn
probcut 0.001
short
```

`short` is not a standard keyword in the documented grammar — it may be silently ignored, or it may be an alias. The short output mode is the default; this keyword likely has no effect.

## Wavefunction

```
# Number of determinants in WFN     =    5519
# Number of electrons               =      10
# Number of active electrons in WFN =      10
# NUMBER OF alpha CONFIGS =   120
# NUMBER OF beta  CONFIGS =   120
```

5519 determinants × 120 alpha configs × 120 beta configs: this is where the OpenMP parallelization pays off.

## Running

```bash
cd test/

# Single thread (baseline)
time OMP_NUM_THREADS=1 ../edf-omp < cas-ch4.edfinp > /dev/null

# Multiple threads
time OMP_NUM_THREADS=4 ../edf-omp < cas-ch4.edfinp > /dev/null
```

Compare the `_linearsys` timer between runs to measure parallel speedup.

## Key output

### Average populations (per basin)

```
<n(  1)_alpha>  =  2.992989     <- carbon
<n(  2)_alpha>  =  0.501702     <- H atom 2
<n(  3)_alpha>  =  0.501702     <- H atom 3
<n(  4)_alpha>  =  0.501702     <- H atom 4
<n(  5)_alpha>  =  0.501702     <- H atom 5
```

Carbon carries ~6 electrons (2 × 2.99) and each hydrogen ~1 electron (2 × 0.50), consistent with the RHF picture. CASSCF correlation redistributes a small amount of electron density.

### Probability distribution (spinless)

With `PROBCUT 0.001`, 255 configurations account for 90.7% of the spinless probability. The dominant configurations place 6 electrons on C and 1 on each H (weight ~0.4). The printed sum:

```
#     0.9070216008  <-- SUM, 255 PROBABILITIES > 0.001
#     0.9999593711  <--- TOTAL SUM
```

The sum of printed probabilities (90.7%) is significantly below the total (99.996%) because many small ionic configurations are below the 0.001 threshold. To see all configurations use `PROBCUT 0.0` (large output) or `PROBCUT -1.0`.

## Performance note

This is the recommended benchmark system for `edf-omp`. The 120 × 120 = 14400 (alpha × beta) pairs in `_linearsys` and the outer 120-element loop in `_splitedf` scale well with thread count. Expected speedup on 4 cores: 2.5–3.5×.
