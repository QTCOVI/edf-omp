# Wavefunction Keywords

These keywords modify or trim the wavefunction before the EDF is computed.

---

## `EPSWFN real`

Discard all determinants with \(|C_i| < \text{real}\) and renormalize the wavefunction. Default is `0.0` (all determinants are kept).

Useful for large CASSCF expansions where many small coefficients contribute negligibly. Combine with `NDETS` to control the effective wavefunction size.

---

## `EPSWFNLOC real`

Same as `EPSWFN` but applied to the wavefunction expressed in terms of **localized** (isopycnic) MOs. Only relevant when `FULLOC` or another localization keyword is active.

---

## `COREMO overlap`

Designates any MO whose self-overlap \(S^A_{ii}\) in any fragment exceeds `overlap` as a **core** orbital for that fragment. Core orbitals are excluded from the EDF computation.

If `overlap` is omitted, the default value is `0.99`.

!!! warning "Core electrons in minelec/maxelec"
    When `COREMO` is active, the `minelec`/`maxelec` limits in Records 5.i must include core electrons — they are total, not valence, electron counts. See [Input Format — Records 5.i](../input-format.md#records-5i-fragment-specification-i-1-ngroup).

---

## `DELETE n1 n2 …`

Remove specific MOs from a multi-determinant wavefunction. The listed indices must all be ≤ `NCORE`. After deletion, the actual number of MOs and core orbitals (alpha and beta) is reduced by the number of valid indices given.

Use `DELETE` to manually excise a known core orbital before the EDF calculation.

---

## `SUPPRESS n1 n2 …`

Write a trimmed AOM file (without the listed MOs) and a trimmed WFN file (without the listed MOs from the WFN), then stop. Does not run the EDF.

`SUPPRESS` is a utility operation for preprocessing: it produces a reduced-dimension AOM+WFN pair for a subsequent edf run that omits the suppressed orbitals.

---

## `LOCORE`

Localize the core MOs of a multi-determinant wavefunction using isopycnic localization, and replace the canonical core MOs with the localized ones in memory. The original canonical MOs are discarded.

Applied before the EDF computation so that the core-localized MOs are used throughout.

---

## `LOCSOME n1 n2 …`

Perform isopycnic localization on the subset of canonical MOs `n1, n2, …` (typically core orbitals in a CASSCF wavefunction). The AOM is updated to reflect the localized MOs.

`LOCSOME` is more selective than `LOCORE`: it localizes only the listed orbitals and leaves the rest unchanged.

---

## `ROHFLOC`

Localize the alpha and beta spin-orbitals of a **ROHF** wavefunction separately using isopycnic localization. Writes a WFN file with the localized MOs and an AOM file with the overlaps between them, then stops.

---

## `UHFLOC`

Same as `ROHFLOC` but for a **UHF** wavefunction.
