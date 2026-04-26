# Integration Keywords

These keywords control the numerical integration grid used when the AOM is computed internally (i.e., when `filedat` is a fuzzy-partition keyword: `mindefrho`, `netrho`, `becke`, `promrho`).

They may also be given directly on the same line as the fuzzy keyword in Record 2. See [Input Format — Record 2](../input-format.md#record-2-filedat).

---

## `LEBEDEV nang`

Number of angular grid points in the Lebedev quadrature used for the angular part of the fuzzy AOM integration. Default is `434`.

Available Lebedev grid sizes are the standard Lebedev orders (e.g., 110, 194, 302, 434, 590, 770, …). Higher values increase accuracy at the cost of compute time.

---

## `NRAD nrad`

Number of radial grid points in the Becke-type integration. Default is `100`.

Higher values improve accuracy for atoms with steep density variations near the nucleus (heavy elements). Values of 200–400 are common for production calculations.

---

## `IRMESH irmesh`

Selects the radial mapping function that transforms the \([0, \infty)\) range to the finite integration interval. Default is `1`.

| `irmesh` | Mapping |
|----------|---------|
| `1` | Becke mapping (default) |
| `2` | Gauss-Chebyshev |

---

## `KPOWER pow`

The iterative K parameter in Becke's fuzzy-cell partitioning scheme, controlling how sharply the weight function \(w_A\) transitions between atoms. Default is `3` (the standard value from Becke's original paper).

Higher values make the partition more step-like (sharper atom boundaries); lower values make it smoother.

Only relevant when `filedat = becke`.

---

## `RHOPOW rhopow`

Exponent applied to each promolecular atomic density \(\rho_A\) when computing the weight:

$$w_A = \frac{\rho_A^{\text{rhopow}}}{\sum_B \rho_B^{\text{rhopow}}}$$

Default is `1.0` (standard Hirshfeld / promolecular weights). Only relevant when `filedat = promrho`.

---

## Inline integration parameters

For `mindefrho`, `netrho`, and `becke`, the integration parameters may be given directly on the `filedat` line in Record 2:

```
mindefrho  NANG  NRAD  IRMESH
netrho     NANG  NRAD  IRMESH
becke      NANG  NRAD  IRMESH  POW
```

These are equivalent to using the `LEBEDEV`, `NRAD`, `IRMESH`, and `KPOWER` keywords later in the input.
