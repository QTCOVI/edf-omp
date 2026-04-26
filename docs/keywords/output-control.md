# Output Control Keywords

These keywords suppress or expand sections of the edf output.

---

## `LARGE` / `VERBOSE` / `BIGOUTPUT`

All three are equivalent: they enable **large output mode**, which prints additional detail throughout the calculation. Default is short output.

Large output includes expanded wavefunction summaries, per-determinant contributions, and extended probability tables.

---

## `NOORDER`

By default, spinless EDF probabilities are sorted in decreasing order before being printed (provided the total count is below the `NSTACK` parameter, currently 4000). `NOORDER` disables this sorting — probabilities are printed in the order they are computed.

Useful when comparing outputs between runs to avoid reordering artifacts.

---

## `NOPROB`

Do not compute the exact EDF for the canonical (input) wavefunction. Useful when only a localized-MO EDF or a population analysis is wanted.

---

## `PROBLOC`

Compute the EDF using the wavefunction expressed in terms of **localized MOs**. By default this is not computed.

Since isopycnic localized MOs are normalized but not orthogonal, the sum of all probabilities in this EDF is generally not 1.0.

Only meaningful when a localization keyword (`FULLOC`, `CANLOC`, etc.) has been given.

---

## `NOPROBX`

Suppress the **approximate** localized-MO EDF (the one computed by assuming all localized MOs are orthogonal within their fragments). By default this approximate EDF is computed alongside `PROBLOC`.

`NOPROBX` is relevant only when a localization has been performed.

---

## `DOENTROPY`

Compute and print the **Mutual Information Entropy (MEI)** after the EDF. The MEI quantifies total inter-fragment electron correlation via:

$$I = \sum_{n_1, \ldots, n_m} p(n_1,\ldots,n_m) \ln \frac{p(n_1,\ldots,n_m)}{\prod_i p_i(n_i)}$$

MEI output is suppressed by default because it can be very large for molecules with many fragments (the full joint probability table is printed).
