# Theory

## Electron Number Distribution Function

Given a molecule partitioned into \(m\) non-overlapping, space-filling fragments \(\Omega_1, \ldots, \Omega_m\), the **Electron Number Distribution Function (EDF)** is the joint probability of finding exactly \(n_1\) electrons in \(\Omega_1\), \(n_2\) electrons in \(\Omega_2\), …, and \(n_m\) electrons in \(\Omega_m\):

$$p(n_1, n_2, \ldots, n_m) = \binom{N}{n_1, n_2, \ldots, n_m} \int_{\Omega_1^{n_1} \otimes \cdots \otimes \Omega_m^{n_m}} |\Psi|^2 \, d\mathbf{r}_1 \cdots d\mathbf{r}_N$$

where the multinomial coefficient accounts for the indistinguishability of electrons and the integration runs over all assignments of \(n_i\) electrons to each region \(\Omega_i\). The EDF satisfies:

$$\sum_{n_1 + \cdots + n_m = N} p(n_1, \ldots, n_m) = 1$$

---

## Atomic Overlap Matrix

The central computational object is the **Atomic Overlap Matrix (AOM)**:

$$S^A_{ij} = \int_{\Omega_A} \phi_i(\mathbf{r}) \, \phi_j(\mathbf{r}) \, d\mathbf{r}$$

where \(\phi_i\) are molecular orbitals (natural orbitals for CASSCF) and \(\Omega_A\) is the basin of atom \(A\). The AOM satisfies the sum rule:

$$\sum_A S^A_{ij} = \delta_{ij}$$

(which edf verifies via `TOLAOM`).

---

## Single-determinant case

For a single Slater determinant of \(N\) spin-orbitals \(\{\phi_i\}\), the EDF probability of configuration \((n_1, \ldots, n_m)\) is a **permanent** of the AOM:

$$p(n_1, \ldots, n_m) = \frac{1}{N!} \sum_{\sigma \in S_N} \prod_{i=1}^N S^{A_{\sigma(i)}}_{ii}$$

where the sum is over all permutations \(\sigma\) that assign \(n_k\) orbitals to fragment \(k\). In practice this is evaluated as a sum of products of matrix permanents.

For the two-group case (\(m = 2\)), a binomial recurrence relation (Cançes, Keriven, Lodier, Savin, *Theor. Chem. Acc.* 111, 373, 2004) allows the full EDF to be computed in \(O(N^2)\) time instead of \(O(N!)\).

---

## Multi-determinant (CASSCF) case

For a CASSCF wavefunction

$$|\Psi\rangle = \sum_I C_I |\Phi_I\rangle$$

the EDF is a double sum over determinant pairs:

$$p(n_1, \ldots, n_m) = \sum_{I,J} C_I C_J \, p_{IJ}(n_1, \ldots, n_m)$$

where \(p_{IJ}\) is the cross-determinant contribution. The total computation scales as \(O(N_\text{det}^2)\) in the number of determinants. The `EPSDET` keyword truncates pairs with \(|C_I C_J| < \epsilon\), and `NDETS`/`EPSWFN` trim the determinant expansion.

---

## Average populations and delocalization indices

**Average population** of fragment \(i\):

$$\langle n_i \rangle = \sum_{\mathbf{n}} n_i \, p(\mathbf{n}) = \sum_j S^i_{jj}$$

(the trace of the fragment AOM over occupied orbitals).

**Two-fragment delocalization index** \(\delta(i,j)\):

$$\delta(i,j) = 2\left[\langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle\right] = -2 \sum_{k,l} S^i_{kl} S^j_{lk}$$

where the second expression holds for single-determinant wavefunctions (Eq. 28 of *J. Chem. Phys.* 126, 094102, 2007). The delocalization index measures electron sharing: \(\delta \approx 1\) for a single covalent bond, \(\delta \approx 2\) for a double bond.

**Localization index**: \(\delta(i,i) = 2[\langle n_i^2 \rangle - \langle n_i \rangle^2]\) is the variance of the electron population in fragment \(i\). The percentage localization is \(\delta(i,i)/\langle n_i \rangle \times 100\).

---

## Fuzzy atomic partitions

When QTAIM or ELF basins are unavailable, edf computes the AOM internally using a **fuzzy partition** in which the sharp basin indicator function is replaced by a smooth weight \(w_A(\mathbf{r}) \geq 0\) with \(\sum_A w_A(\mathbf{r}) = 1\):

$$S^A_{ij} = \int w_A(\mathbf{r}) \, \phi_i(\mathbf{r}) \, \phi_j(\mathbf{r}) \, d\mathbf{r}$$

Available weight functions:

| Keyword | Weight \(w_A\) |
|---------|---------------|
| `mulliken` | Mulliken (step at AO center) |
| `lowdin` | Löwdin (symmetric orthogonalization) |
| `mindef` | Minimally Deformed Atoms (Fernández Rico et al.) |
| `mindefrho` | MinDef with \(w_A = \rho_A/\rho\) |
| `netrho` | Net density (eliminates cross-center terms) |
| `promrho` | Promolecular: \(w_A = \rho_A^{\rm atom}/\rho\) |
| `heselmann` | Heselmann chemical localization weights |
| `becke` | Becke fuzzy cells |

The numerical integrals use Lebedev angular quadrature (controlled by `LEBEDEV`) and a radial Becke-type grid (`NRAD`, `IRMESH`).

---

## Key references

- Francisco, E.; Martín Pendás, A.; Blanco, M. A. "EDF: Computing electron number probability distribution functions in real space from molecular wave functions." *Comput. Phys. Commun.* 178, 621 (2008).
- Cançes, E.; Keriven, R.; Lodier, F.; Savin, A. "How electrons guard the space." *Theor. Chem. Acc.* 111, 373 (2004).
- Bader, R. F. W. *Atoms in Molecules: A Quantum Theory.* Oxford University Press, 1990.
- Becke, A. D. "A multicenter numerical integration scheme for polyatomic molecules." *J. Chem. Phys.* 88, 2547 (1988).
- Heselmann, A. "Improved description of chemical bonding from a localized molecular orbital perspective." *J. Chem. Theory Comput.* 12, 2720 (2016).
- Matito, E.; Solà, M.; Salvador, P.; Duran, M. "Electron sharing indexes at the correlated level." *Faraday Discuss.* 135, 325 (2007). [Eq. 28 for DI]
