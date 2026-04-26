# Keywords Overview

After the three positional records and the fragment specification, all remaining input is keyword-driven. Keywords are case-insensitive. Lines whose first non-blank character is `#` are comments. `END` terminates the input.

---

## Categories

### [Basic Control](basic.md)

Control thresholds, truncation, and the core computation options.

| Keyword | Brief description |
|---------|------------------|
| `AOMNORM [atom]` | Renormalize the AOM for the given atom |
| `PROBCUT real` | Omit probabilities below this threshold from output |
| `EPSDET real` | Skip determinant pairs \|Ci·Cj\| below this threshold |
| `EPSPROBA real` | Skip beta computation when max alpha probability is below this |
| `TOLAOM real` | Tolerance for AOM sum-rule check (default 0.01) |
| `NDETS n` | Use only the first `n` determinants |
| `RECUR` | Force binomial recurrence path (default for ngroup = 2) |
| `NOMEM` | Use disk instead of memory for intermediate probabilities |
| `RANDOM min max` | Set random bounds for the `binedf.f` linear system |

### [Wavefunction](wavefunction.md)

Modify or trim the wavefunction before the EDF is computed.

| Keyword | Brief description |
|---------|------------------|
| `EPSWFN real` | Drop determinants with \|coeff\| < epswfn |
| `EPSWFNLOC real` | Same but applied to the localized-MO expansion |
| `COREMO real` | Treat MOs with self-overlap > real as frozen core |
| `DELETE n1 n2 …` | Remove specific core MOs from the wavefunction |
| `SUPPRESS n1 n2 …` | Write a trimmed AOM+WFN pair and stop |
| `LOCORE` | Localize core MOs of a CASSCF wavefunction in place |
| `LOCSOME n1 n2 …` | Isopycnic-localize a subset of canonical MOs |

### [Localization](localization.md)

Perform an isopycnic (Cioslowski) localization of molecular orbitals.

| Keyword | Brief description |
|---------|------------------|
| `FULLOC [ALLCENTER]` | Localize natural orbitals |
| `CANLOC [ALLCENTER]` | Localize canonical orbitals |
| `CHEMLOC` | Chemical localization variant |
| `AOMLOC` | AOM-based chemical localization |
| `DAFHLOC` | DAFH-based localization |
| `OPENLOC` | OQS 1-RDM diagonalization localization |
| `DODAFH … END` | Full control DAFH localization block |
| `OTHEROPEN … END` | OQS 1-RDM localization block |
| `DAFHSELECT n n1 n2 …` | Restrict DAFH search to specific atom tuples |
| `CRITICS real` | Threshold for partial localization detection |
| `COVX real` | Covalent radius scaling for bond detection |
| `DAMPS real` | Damping factor for CHEMLOC localization cycles |
| `MXCENT n [n1 n2 …]` | Max centers in progressive DAFH localization |
| `CRITOV [n critval]` | Overlap threshold for MO localization acceptance |
| `ROHFLOC` | Localize ROHF spin-orbitals separately and stop |
| `UHFLOC` | Localize UHF spin-orbitals separately and stop |

### [Analysis](analysis.md)

Optional population analyses and structure searches.

| Keyword | Brief description |
|---------|------------------|
| `PRSRS n1 n2 …` | Probability of a specific resonance structure |
| `MAXPOP m1 m2 …` | Max α/β electrons per group for PRSRS filter |
| `MINPOP m1 m2 …` | Min α/β electrons per group for PRSRS filter |
| `TWOCENDI [eps]` | Compute all two-center delocalization indices |
| `OQSEOS` | Open quantum system EOS analysis and stop |
| `OQSLEWIS` | Lewis structure search via OQS and stop |
| `FNO atom1 atom2 …` | Fragment Natural Orbitals and stop |
| `DIRPROB` | Brute-force probability computation (SDW only) |
| `DOENTROPY` | Compute mutual information entropy |
| `BONDING … ENDBONDING` | Fit EDF to a bond model |

### [Integration](integration.md)

Grid parameters for fuzzy AOM integration (used with `mindefrho`, `netrho`, `becke`).

| Keyword | Brief description |
|---------|------------------|
| `LEBEDEV nang` | Angular Lebedev grid points (default 434) |
| `NRAD nrad` | Radial grid points (default 100) |
| `IRMESH irmesh` | Radial mapping scheme (default 1) |
| `KPOWER pow` | Becke iterative K parameter (default 3) |
| `RHOPOW rhopow` | Exponent for promolecular density weights |

### [Output Control](output-control.md)

Suppress or expand output sections.

| Keyword | Brief description |
|---------|------------------|
| `LARGE` / `VERBOSE` / `BIGOUTPUT` | Enable extended output |
| `NOORDER` | Do not sort spinless probabilities by decreasing value |
| `NOPROB` | Do not compute the EDF for the canonical wavefunction |
| `PROBLOC` | Compute EDF for the localized-MO wavefunction |
| `NOPROBX` | Suppress the approximate localized-MO EDF |
| `DOENTROPY` | Print mutual information entropy |
