# Critical three-state Potts levels

`bethe-potts-pbc` computes the finite-ring vacuum and a selected charged
one-hole branch at every lattice momentum. These are solutions of finite-size
Bethe equations, not energies obtained by inserting conformal dimensions into
an asymptotic formula. The branch supplies a useful discrete-symmetry benchmark
for MPS excitation calculations.

```sh
bethe-potts-pbc 16 --json potts.json --csv potts.csv
bethe-potts-pbc 64 --momentum-index 1 --precision long-double
bethe-potts-pbc 32 --branch charged --charge -1 --tsv charged.tsv
bethe-potts-pbc 128 --branch ground --precision fp128
bethe-potts-pbc --references
```

fp64 is the default; native long-double and optional MPLAPACK-enabled fp128
arithmetic are supported. [Build instructions](building.md) and
[output conventions](output.md) apply.

## Hamiltonian and physical labels

For a ring of L clock sites, use

```math
H=-\sum_{j=1}^{L}\left(X_j+X_j^\dagger+
 Z_jZ_{j+1}^\dagger+Z_j^\dagger Z_{j+1}\right),\qquad Z_{L+1}=Z_1.
```

Here $`X|a\rangle=|a+1\pmod3\rangle`$ and
$`Z|a\rangle=\omega^a|a\rangle`$, with
$`\omega=e^{2\pi i/3}`$. The exchange and critical transverse field are both
one. The two oriented bonds are both included when L=2. Off-critical, chiral,
antiferromagnetic and nonperiodic Potts chains are outside this frontend.

`charge` is q=0,+1,-1, the exponent in the eigenvalue
$`\omega^q`$ of $`\prod_j X_j`$. `k` labels physical translation by one
clock site, with momentum $`p=2\pi k/L`$. Charge conjugation gives equal
energies at q=+1 and q=-1, at the **same** momentum. Spatial reflection relates
k and L-k. These are distinct symmetries; do not conflate their multiplicities.

The Hamiltonian is $`\sqrt3/2`$ times the normalization of
[Dasmahapatra et al.](../CITATIONS.md#dasmahapatra-kedem-mccoy-melzer-1994).
In our units the exact thermodynamic energy density and velocity are

```math
e_\infty=-\frac43-\frac{2\sqrt3}{\pi},\qquad
v=\frac{3\sqrt3}{2}.
```

## Why use an auxiliary XXZ problem?

The Potts and twisted-XXZ Hamiltonians represent the same Temperley–Lieb
generators, with loop weight $`\sqrt3`$. The auxiliary XXZ ring has **2L**
spin-half sites and $`\Delta=\sqrt3/2`$. With the usual spin-half XXZ energy
including its ferromagnetic reference constant, the map is

```math
E_{\mathrm{Potts}}=2\sqrt3\,E_{\mathrm{XXZ}}+\frac L2.
```

This follows by writing the clock Hamiltonian as
$`-\sqrt3\sum_{j=1}^{2L}e_j+2L`$ and using the XXZ realization of the same
generators. One clock translation advances two auxiliary sites. The map and
its representation-theory qualifications are discussed in
[Fukai et al., Sections 2.2–2.3](https://arxiv.org/html/2309.07472v2) and
[Nichols, Section 6](https://arxiv.org/html/hep-th/0509069v1).

Let $`\gamma=\pi/6`$ and store the L real auxiliary roots in the bounded
coordinate $`z_j=\tanh\lambda_j/\tan(\gamma/2)`$. The shared XXZ solver
solves

```math
\begin{aligned}
F_i&=4L\arctan z_i-2\pi I_i-\phi-2\sum_{j\ne i}\alpha_{ij}=0,\\{}
\alpha_{ij}&=\arctan\frac{\Delta(z_i-z_j)}
 {1+\Delta-(1-\Delta)z_i z_j},\\{}
E_{\mathrm{XXZ}}&=-\frac{L\Delta}{2}
 +\sum_{j=1}^{L}\frac{z_j^2-1}{z_j^2+1}.
\end{aligned}
```

Roots stay strictly inside the finite-real branch. No finite-size string
deviations are discarded. The internal twist parameter extends the shared
iteration; it does not widen the public XXZ explicit-label solver's supported
window or turn it into an arbitrary-twist spectrum enumerator.

## Physical-root selection and counting audit

**An arbitrary converged twisted-XXZ solution is not a physical Potts level.**
At this root of unity, auxiliary representations contain submodules which
must be excluded or combined differently in the Potts representation. Even
the small-ring decomposition is not simply the entire zero-magnetization XXZ
spectrum. A regression explicitly constructs an excluded level in the vacuum
twist sector by adding a root at negative infinite rapidity to a zero-twist
XXZ state with one fewer root. Its Bethe equations hold and its energy is
unchanged, but that energy is absent from every charge/momentum block of the
four-site clock Hamiltonian. Root convergence is not a physical-state test.

The public selection is fixed as follows:

| Family | Auxiliary twist | Occupied Bethe labels | Physical output |
| --- | --- | --- | --- |
| Vacuum | $`\phi=\pi/3`$ | $`I_j=j-(L-1)/2,\ j=0,\ldots,L-1`$ | One q=0, k=0 state |
| Charged one-hole | $`\phi=\pi`$ | $`I_j=j-(L+1)/2,\ j=0,\ldots,L`$, omit j=h | One state at k=h for each q=+1,-1 |

Use **0≤h<L**, not 0≤h≤L: the second endpoint has the same physical momentum
and energy as h=0 and is not another level of either charge. The retained
family contains **2L charged levels**, or **2L+1 rows including the vacuum**.
At even L, k=L/2 is its own reflection partner; it still has two charge copies,
not four. Summing the occupied labels gives $`\sum I_j=-h`$, consistent
with the physical momentum map above.

This is a **selected family**, not a complete decomposition of the Potts
Hilbert space of dimension $`3^L`$. The charged roots at L=4 reproduce the
real-root examples and their Potts identifications in Fukai et al., Tables 3
and 11. Independent clock-Hamiltonian tests resolve both charges and all
momenta for L=2 through 6. The family is the lowest charged level in each of
those tested blocks; the API does **not** claim a proof of a global sector
minimum for arbitrary L. Convergence alone is never used to classify an
arbitrary auxiliary solution as physical.

The direct Potts Bethe equations use a different set of complex roots,
including a paired-root ground sea and auxiliary “ghost” species involved in
counting. Their bookkeeping is not interchangeable with these L auxiliary
XXZ roots. The audit of those species in Dasmahapatra et al., Section 2 and
Appendix A, is the reason this frontend does not offer generic direct-root
enumeration. **Neutral excited levels**, other charged families, arbitrary
root labels and full-spectrum `--excitations all` remain unimplemented.

## Output and finite-size checks

The `levels` table has `branch`, `charge`, `k`, `momentum`, `energy`, `gap`,
`x_scaled`, `residual`, `iterations`, and `status`. All gaps subtract the
separately calculated finite-ring vacuum, even with `--branch charged`.
`--branch all` means all **supported families**, not all physical excitations.
`--momentum-index` and `--charge` filter only the charged rows; the vacuum
remains the reference and is output unless `--branch charged` is requested.

The scaled gap is an estimator,

```math
x_L=\frac{L(E-E_0)}{2\pi v},\qquad
c_L=-\frac{6L(E_0-Le_\infty)}{\pi v}.
```

The vacuum approaches central charge c=4/5. The charged k=0 state approaches
scaling dimension x=2/15; k=1,L-1 approach its first descendants at x=17/15.
Finite-size corrections are retained, not fitted away. At fixed physical
momentum the branch tends to $`3\sqrt3\sin(p/2)`$ for
$`0\le p\le2\pi`$. Energies alone provide no spectral weights or
guarantee that a particular observable couples to a level.

The normalized residual is $`\max_i|F_i|/(2L)`$, in radians. The default
tolerance is 32 machine epsilons; `--max-iterations` limits sweeps **per solve**.
A zero budget evaluates the initial roots but does not manufacture a level.
Unconverged energies are null/empty; gaps also require a converged vacuum.
Incomplete numerical output exits 2, invalid arguments exit 1. Invalid inputs
are checked before opening or overwriting exports. The default site budget is
512, adjustable with `--max-sites`; each sweep costs O(L²), and a full momentum
scan requires L charged solves. Charge conjugation reuses a solve.

Validation includes exact L=2 expressions, independent 75-digit clock-matrix
references at L=3, all momentum/charge blocks through L=6, reflected roots,
exponentiated Bethe equations, conformal scaling at L=16,32,64, and explicit
failure/serialization contracts. Separate sparse-oracle audits at L=8 and 9
match all selected levels to within 1.4e-13 in double precision.
The independent oracle can be rerun without invoking the Bethe code:

```sh
python3 scripts/reference_potts.py --sites 3 --digits 75
OPENBLAS_NUM_THREADS=1 python3 scripts/reference_potts.py --sites 8 --backend scipy
```

The C++ entry points are `potts::ground_state<Real>(L)` and
`potts::charged_one_hole<Real>(L, charge, k)` in
[`bethe/potts.hpp`](../include/bethe/potts.hpp). The result keeps auxiliary roots
and half-integer labels for inspection. Its energy is a last-iterate diagnostic
unless `converged` is true; the frontend withholds unconverged observables.
