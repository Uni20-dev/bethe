# Sutherland gas on a circle

`bethe-sutherland-pbc` evaluates the finite-volume spectrum of the scalar,
bosonic trigonometric inverse-square gas. Here “Sutherland” means the continuum
Calogero–Sutherland family, not the Uimin–Lai–Sutherland spin chain. Energies
follow explicit spectral rules: there are no iterative root solves or
convergence tolerances. fp64, long-double and optional fp128 are supported.

```sh
build/bethe-sutherland-pbc 8 --length 8 --lambda 2
build/bethe-sutherland-pbc 4 --length 4 --lambda 0.5 --labels -1,0,0,2
build/bethe-sutherland-pbc 4 --length 4 --lambda 2 --levels 10 --window 2
build/bethe-sutherland-pbc 4 --length 4 --lambda 2 --levels all --window 2 --format csv
build/bethe-sutherland-pbc 4 --length 4 --lambda 2 --pseudomomenta --precision fp128
```

## Hamiltonian and collision branch

With $\hbar =2m =1$, circumference $L >0$ and $N \ge 1$, our convention is

```math
H=-\sum_i\partial_i^2
+2\lambda(\lambda-1)\left(\frac\pi L\right)^2
\sum_{i<j}\frac{1}{\sin^2[\pi(x_i-x_j)/L]}.
```

The input $\lambda \ge 0$ selects the bosonic Jastrow domain, with wavefunctions
behaving as $\lvert x_{i} -x_{j} \rvert ^\lambda$ at a collision. It is **not** the raw potential
coefficient. In particular:

- $\lambda =0$: ordinary free bosons, allowing particles to coincide.
- $\lambda =1$: hard-core bosons; the collision boundary condition remains
  different even though the explicit potential coefficient is also zero.
- $0<\lambda <1$: negative inverse-square coefficient on this specified branch,
  not an arbitrary attractive potential with an unspecified boundary condition.

Replacing $\lambda$ by $1-\lambda$ leaves the potential unchanged but generally
changes the domain and spectrum. We do not choose a branch by inverting that
coefficient. At $\lambda =1/2$, the Jastrow branch excludes the logarithmic
collision solution. Other self-adjoint extensions are outside this interface.

The normalization and partition spectrum are given in
[Gurappa–Panigrahi, Eqs. (18)–(23)](https://arxiv.org/html/hep-th/9908127).
Their coupling $\beta$ is our $\lambda$; their descending partition entries
become our ascending labels after reversal and an optional integer boost.
The original solution is [Sutherland (1971)](../CITATIONS.md#sutherland-1971).

## Labels, pseudomomenta and gaps

Specify $N$ nondecreasing integers `n_0 <= ... <= n_(N-1)`. Repeated and negative
labels are allowed. The ground state has every label zero. With $q =2\,\pi /L$,

```math
\begin{aligned}
k_j&=q\left[n_j+\lambda\left(j-\frac{N-1}{2}\right)\right],
\quad j=0,\ldots,N-1,\\
E&=\sum_j k_j^2,\qquad \mathrm{momentum\_index}=\sum_jn_j,\\
P&=q\,\mathrm{momentum\_index},\\
E_0&=\left(\frac\pi L\right)^2\frac{\lambda^2N(N^2-1)}3.
\end{aligned}
```

Continuum momentum is **not reduced modulo $2\,\pi$**. The $k_{j}$ are
pseudomomenta, not separately occupied free-particle orbitals except in the
appropriate free limits. At $\lambda =1$, the free-fermion dual has integer
momenta for odd $N$ and half-integer, antiperiodic momenta for even $N$, while
the physical bosons always have periodic boundary conditions.

A uniform integer boost $n_{j} \to n_{j} +b$ changes $P$ by $N \,q \,b$ and $E$ by
$2\,q \,b \,P +N \,q ^{2}\,b ^{2}$. Reflecting and reversing the labels reverses $P$ without
changing the energy. Changing only $L$ scales momentum as $1/L$ and energy
as $1/L ^{2}$.

The implementation computes the gap directly from nonnegative terms,

```math
\mathrm{gap}=q^2\left[\sum_jn_j^2+
\lambda\sum_{j=1}^{N-1}j(N-j)(n_j-n_{j-1})\right].
```

This avoids subtracting nearly equal total and ground energies. A small boost
gap can remain accurate even when `E0+gap` rounds back to `E0`. Integer
coefficients are accumulated before the scalar prefactors, preserving ordinary
small-label degeneracies. Very large labels and nearly equal distinct gaps
remain subject to the chosen floating-point precision.

## What an excitation scan covers

`--levels COUNT --window W` returns the lowest COUNT states **within**
$-W \le n_{j} \le W$, including the ground state. `--levels all` returns every ordered
label list in that window. One row is one bosonic eigenstate; equal energies
are separate rows, and COUNT can cut through a degeneracy. Ordering uses the
computed gap, then total momentum index, then lexicographic labels.

There are `binomial(N+2*W,N)` states in the window. Before allocating them,
the tool checks this against `--max-states` (default 100,000). Exceeding that
budget returns exit status 2 with no state rows. It does not publish a partial
search as the lowest levels. Every accepted scan visits the whole window before
truncation; memory grows as $N$ times the number of states, even for a small
requested COUNT. Increase budgets with that cost in mind.

The continuum spectrum is infinite. No finite-window result claims global
low-energy completeness. Increase the window and compare if states outside it
could matter. Ground-state and explicit-label evaluation need no scan.

The implementation bounds $N$ at 1,000,000 and $W$ at 1,000,000. Explicit
labels and their total momentum index must fit signed 64-bit integers; a wide
internal sum allows cancellation of prefixes outside that range. Invalid
inputs, nonfinite computed results and positive energies/gaps underflowing to
zero are errors (CLI exit 1), not valid zero energies. As with other numerical
tools, extreme intermediate arithmetic can exceed the selected scalar range.

Output formats are `auto`, `pretty`, `plain`, `csv`, `tsv`, and `json`, with
independent file exports through `--csv FILE`, `--tsv FILE`, and `--json FILE`.
The primary `states` table has columns `state_id,labels,energy,gap,momentum_index,p`;
integer label lists remain space-separated. `--pseudomomenta` adds a separate
typed table `state_id,index,label,k`, joined to states by zero-based `state_id`.
Use `--table pseudomomenta` to select it for CSV/TSV, or export both tables as
JSON. See [output conventions](output.md) for metadata, streaming and table selection.
Use `--references` for centralized literature, also listed in [CITATIONS.md](../CITATIONS.md).

## Library and validation

```cpp
#include <bethe/sutherland.hpp>
namespace cs = bethe::sutherland;
auto ground = cs::ground_state(8, 8.0, 2.0);
auto excited = cs::evaluate(4, 4.0L, 0.5L, {-1, 0, 0, 2});
auto window = cs::spectrum(4, 4.0, 2.0, 2, 10); // Check window.complete.
```

Tests apply the coordinate-space differential Hamiltonian independently to
many-particle Jastrow ground states, degree-two Jack excitations and a family
of two-particle Gegenbauer excitations. Further checks cover free-boson and
hard-core limits, boosts, reflection, length scaling, window counts, range
errors, and native-precision CLI parsing/output. These run for every enabled
scalar type, without narrowing comparisons to double.

This first slice does not expose wavefunctions, matrix elements, spinful or
fermionic variants, twists, trapped/rational Calogero systems, thermodynamic
response functions, or arbitrary collision domains.
