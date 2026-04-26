# Analysis Keywords

These keywords request population analyses, resonance structure probabilities, and structure searches beyond the standard EDF.

---

## `PRSRS n1 n2 … nG`

Compute the probability of finding exactly `n1` electrons in group 1, `n2` in group 2, …, `nG` in group G — a specific **Real-Space Resonance Structure (RSRS)**.

Multiple `PRSRS` lines may appear in the input to request several structures at once.

**Example** — probability of the (2,2,2) configuration in a 3-group molecule:
```
PRSRS 2 2 2
```

---

## `MAXPOP m1 m2 … mG`

When computing PRSRS, restrict the spin-resolved search to configurations where group `i` has at most `mi` alpha or beta electrons. Must appear after a `PRSRS` line.

---

## `MINPOP m1 m2 … mG`

When computing PRSRS, restrict the spin-resolved search to configurations where group `i` has at least `mi` alpha or beta electrons. Must appear after a `PRSRS` line.

---

## `TWOCENDI [epsbond]`

Compute all **two-center delocalization indices** \(\delta(A,B)\) between every pair of atoms. `epsbond` is a threshold: atoms A and B are reported as bonded if \(\delta(A,B) > \text{epsbond}\). Default is `1.0`.

The calculation is run ten times with thresholds `epsbond`, `epsbond/2`, …, `epsbond/2⁹`.

!!! warning
    Only works with **single-determinant wavefunctions**. Not available for CASSCF.

---

## `OQSEOS`

Perform an **Open Quantum System Effective Oxidation State (OQS-EOS)** analysis of each fragment and stop. Computes the EOS of each fragment from the 1-RDM of the open system.

Requires a prior alpha and beta 1-RDM computation. Only closed-shell coupled-cluster wavefunctions are supported for CCWFN inputs; CASSCF wavefunctions are also supported.

---

## `OQSLEWIS`

Search for **Lewis structures** based on the user-defined fragments using the OQS approach. Control passes to `lewisoqs.f`, and edf stops after the search.

Use `NFRAG` to define a fragment count different from `NGROUP` for the Lewis search.

---

## `FNO atom1 atom2 …`

Compute the **Fragment Natural Orbitals (FNOs)** for the fragment formed by atoms `atom1`, `atom2`, …. Multiple `FNO` lines may appear. After all FNO computations, edf stops.

FNOs are the natural orbitals of the reduced density matrix of the fragment, obtained by diagonalizing the fragment's one-body density matrix within the AOM framework.

---

## `DIRPROB`

Compute probabilities by a **direct brute-force** enumeration over all real-space resonance structures whose electron counts fall within `[minpopul(i), maxpopul(i)]` for each fragment. Stopping criterion after all structures are found.

!!! warning
    Only works for **single-determinant wavefunctions**. Scales exponentially with molecule size — use only for small systems with few fragments.

---

## `DOENTROPY`

Compute the **Mutual Information Entropy (MEI)** after the EDF is computed. MEI quantifies the total electron correlation between the defined fragments.

Default is to skip this — MEI output can be very large for many fragments. See also [Output Control — DOENTROPY](output-control.md#doentropy).

---

## `BONDING … ENDBONDING`

Fit the molecular EDF to a model composed of 2-center-2-electron (2c,2e) bonds and optionally 3-center-2-electron (3c,2e) bonds. The block ends with `ENDBONDING`.

This fitting minimizes a weighted objective:

$$\Delta = W_\text{EDF} \cdot A + W_\text{POP} \cdot B + W_\text{DIS} \cdot C$$

where \(A\) is the mean squared difference between model and exact probabilities, \(B\) is the population error, and \(C\) is the delocalization-index error.

**Keywords within the block:**

| Keyword | Syntax | Meaning |
|---------|--------|---------|
| `ENDBONDING` | — | End the block (mandatory) |
| `TYPE` | `type qq ifixq ffty ifixf` | Define a 2c,2e bond type |
| `TYPE3C` | `type p200 p020 p002 p110 p101 f200 f020 f002 f110 f101` | Define a 3c,2e bond type |
| `PAIR` | `frag1 frag2 type [rtype]` | Place a 2c,2e bond between fragments |
| `TRIO` | `frag1 frag2 frag3 type [rtype]` | Place a 3c,2e bond between fragments |
| `EPSBOND` | `real` | Convergence threshold (default 1e-4) |
| `INICVAR` | `real` | Initial step parameter (default 0.1) |
| `PRIN` | `int` | Print level (default 0) |
| `WEDF` | `real` | Weight of EDF term (default 1.0) |
| `WPOP` | `real` | Weight of population term (default 0.0) |
| `WDIS` | `real` | Weight of DI term (default 0.0) |

For `TYPE`: `qq` is the bond polarity (optimized if `ifixq ≠ 0`); `ffty` is the correlation factor (optimized if `ifixf ≠ 0`). Default starting values of `qq = 0.5`, `ffty = 0.0` are typical.

**Example** — two H–H bond types, one bond:
```
BONDING
  TYPE  1  0.5  1  0.0  1
  PAIR  1  2  1
ENDBONDING
```
