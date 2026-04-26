# Output Guide

edf writes all output to **stdout**. A typical run produces five sections:

1. Header and input echo
2. Wavefunction summary
3. EDF probability tables
4. Populations and delocalization indices
5. Timer

---

## 1. Header and input echo

The output begins with an ASCII banner, then echoes the input file verbatim between `# The Input file is` and `# End of the Input file` markers. Verify this against your input to catch silent truncation or encoding issues.

---

## 2. Wavefunction summary

```
# AOM file                          = h6.aom
# Wave Function file                = h6.wfn
# Input number of Primitives        48 reduced to       48
# Description of the Primitive Basis Set
# Total number of Primitive Gaussians:     48
# CENTER   1
# S Shell (Z=0.33865E+02) :    1
...
```

This section lists:
- The AOM and WFN filenames.
- The basis set (number of Gaussian primitives per shell per center).
- For CASSCF: the number of determinants, core/active/virtual orbital counts, and the CI coefficients.

If the primitive count after `reduced to` differs from `Input number`, some primitives were eliminated by the wavefunction reader — check `TOLAOM` if the AOM normalization is suspect.

---

## 3. EDF probability tables

### Spin-split (full) EDF

```
# Spin-Splitted probabilities
# M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION INCLUDING SPIN
# NUMBER OF GROUPS               =                  6
# TOTAL NUMBER OF PROBABILITIES  =               3136
# Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS IN GROUP i
# -------------------------------------------------------------------------------------
#     Probability           G1(a) G1(b) G2(a) G2(b) ...
#     0.0000000075662346      3     3     0     0     0     0     0     0     0     0     0     0
#     0.0000001158377170      3     2     0     1     0     0     0     0     0     0     0     0
```

Each row: probability followed by the alpha (`Gi(a)`) and beta (`Gi(b)`) electron counts for each group. The total number of rows is `(Nα+1)^G × (Nβ+1)^G` before truncation by `PROBCUT`.

The section ends with:
```
#     0.9994921681758134  <- SUM (     3136 PROBS with P_{a} & P_{b} >  0.0E+00)
#     0.9994921681758134  <- TOTAL SUM
```

The sum should be close to 1.0; the deviation reflects the probability carried by configurations below the `PROBCUT` threshold.

### Alpha spinless EDF

```
# M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION FOR ALPHA ELECTRONS
# FOR EACH VALUE, A SUM OVER ALL BETA  RESONANT STRUCTURES HAS BEEN DONE
# NUMBER OF GROUPS                                  =        6
# TOTAL NUMBER OF PROBABILITIES FOR ALPHA ELECTRONS =       56
#     Probability            n1    n2    n3 ...
#     0.0000869620157188      3     0     0     0     0     0
...
#     0.9994921681758127  <-- SUM,      56 PROBABILITIES >  0.000000000E+00
```

The alpha EDF is obtained by summing the spin-split EDF over all beta configurations for each alpha configuration. Similarly, the beta EDF (printed next) sums over alpha configurations.

### Spinless (total) EDF

After the alpha and beta tables, a **total-electron spinless EDF** sums both spins:

```
#     Probability            n1    n2    n3 ...
#     0.0002607553   1  1  1  1  1  1    <- dominant: 1 electron per H
#     0.0001234281   2  1  1  1  1  0
```

For a closed-shell molecule with equal alpha and beta, the spinless EDF is the convolution of alpha and beta EDFs.

---

## 4. Populations and delocalization indices

### Average populations

```
<n(  1)_alpha>                  =      0.499746084
<n(  1)_beta>                   =      0.499746084
<n(  2)_alpha>                  =      0.499760074
...
```

\(\langle n_i \rangle_\alpha\) is the average number of alpha electrons in fragment \(i\), and similarly for beta. For a closed-shell molecule these are equal. The total population is \(\langle n_i \rangle = \langle n_i \rangle_\alpha + \langle n_i \rangle_\beta\).

### Two-body expectations

```
<n(  2)_alpha n(  1)_alpha>     =      0.139042602
<n(  2)_alpha n(  1)_beta>      =      0.249880037
```

Mixed two-body populations \(\langle n_i n_j \rangle\) enter the delocalization index:

$$\delta(i,j) = 2\left[\langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle\right]$$

### Delocalization indices

```
 Delocalization indices, Eq. (28) J. Chem. Phys.  126, 094102 (2007)
 # delta_(  2  1)         =      0.443349226
 # delta_(  3  1)         =      0.060120455
 # delta_(  3  2)         =      0.443399770
```

Two-center DIs measure electron sharing: \(\delta(i,j) \approx 1\) for a covalent bond, \(\delta(i,j) \ll 1\) for non-bonded pairs. Three-center DIs \(\delta(i,j,k)\) quantify three-center bonding.

The localization index for fragment \(i\) is:

```
# delta_(  1  1)         =      0.439111188  % Localization =  43.9334
```

\(\delta(i,i) = 2 \left[\langle n_i^2 \rangle - \langle n_i \rangle^2\right]\) measures how "localized" the electrons of fragment \(i\) are; the percentage is \(\delta(i,i)/\langle n_i \rangle \times 100\).

### Covariance matrix

```
# Covariance Matrix
      0.56038098     -0.22167461 ...
# Covariance Eigenvalues
      0.00050886      0.42511057 ...
# Covariance Eigenvectors
      0.40824829      0.40827125 ...
```

The covariance matrix \(\sigma_{ij} = \langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle\) (diagonal entries are variances). Its eigenvalues and eigenvectors encode the collective electron fluctuation modes across fragments.

---

## 5. Timer

```
#    timer:
#    -pid----name--------------cumtime-------pcalls--popen-
#     1     _edf               0.017662          1     T
#     4     _binedf            0.014454          1     F
#     6     _linearsys         0.000534          1     F
#     7     _splitedf          0.000009          1     F
```

| Field | Meaning |
|-------|---------|
| `cumtime` | Cumulative wall time in seconds |
| `pcalls` | Number of times this timer was entered |
| `popen` | Whether this timer is still open at exit |

`_linearsys` and `_splitedf` are the parallelized regions in `binedf.f`. Compare their wall times between 1-thread and N-thread runs to measure parallel speedup. See [Performance](performance.md).
