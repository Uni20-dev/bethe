# Repulsive SU(n) fermions on a ring

[Overview](../README.md) · [Gaudin–Yang](gaudin-yang.md) · [Model catalogue](models.md)

`bethe-sun-fermions-pbc` generalizes the two-component Gaudin–Yang gas to
an arbitrary number of internal components. To distinguish the component
count from the particle count N, this guide calls it κ. All particles have
the same mass and every pair of different components has the same contact
repulsion:

```math
\begin{aligned}
H&=-\sum_j\frac{\partial^2}{\partial x_j^2}+2c\sum_{i<j}\delta(x_i-x_j),\\
\frac{\hbar^2}{2m}&=1,\quad c\ge0,\quad \ell>0\ \text{(circumference)},\\
E&=\sum_jk_j^2,\qquad P=\sum_jk_j.
\end{aligned}
```

Identical-component contact interactions vanish by antisymmetry. There
are no component-dependent fields, traps or energy shifts. This is a
continuum gas, not the [SU(3) permutation spin chain](su3.md).

## First calculations and supported sectors

```sh
build/bethe-sun-fermions-pbc --populations 3,3,3 --length 9 --c 1
build/bethe-sun-fermions-pbc --populations 5,3,1 --length 9 --c 0.01 --roots
build/bethe-sun-fermions-pbc --populations 3,1,1,1 --length 6 --c 2 --precision long-double
# Requires a binary128-enabled build:
build/bethe-sun-fermions-pbc --populations 1,1,1 --length 1 --c 1 --precision fp128
```

The interacting ground-state branch currently requires **every occupied
component to contain an odd number of particles**. Both balanced and
imbalanced populations are supported. For example `3,3,3`, `5,3,1` and
`3,1,1,1` work; `2,2,2` and `2,1,1` do not. This is a finite-ring
shell/Bethe-number restriction, not a failure of the model's integrability.
Other shell branches need their own state selection and validation.

Empty components are allowed: `1,0,3,1` has the same energy as `3,1,1`.
The physical component order is preserved in the result. For nesting,
occupied components are stably sorted by decreasing population; the report
shows the mapping, with zero-based indices and rank 0 as the reference
component. Empty components have no Bethe sea.

At **c=0**, or with **at most one occupied component**, arbitrary
nonnegative counts are exact free limits. A vacuum such as `0,0,0` has
E=P=0. Each free component fills integer modes
`j-floor((N_a-1)/2)`, `j=0,...,N_a-1`. An even free population selects the
positive-current representative of the degenerate ground state. Free
occupations are not interacting Bethe labels; their tables are separate.

The existing `bethe-gaudin-yang-pbc` remains the convenient two-spin
frontend with `--sz`. This frontend instead takes an explicit population
list. It does not implement Bose–Fermi mixtures, attraction, excitations,
open ends, unequal masses or unequal pair couplings.

## Nested Bethe equations

The original multicomponent solution is due to
[Sutherland](../CITATIONS.md#sutherland-1968). Our finite-size equations use
[Lee–Guan–Batchelor, Eqs. (2)–(3)](https://arxiv.org/pdf/1011.0128), with
`c'=c/2`; their logarithmic equations (13)–(19) fix the nesting signs.
We define c by the **positive `2c delta` interaction above**, without an
additional scattering-length convention. No string hypothesis or
thermodynamic approximation is made in the finite-size solver.

After sorting the occupied populations `N_1>=...>=N_κ`, let

```math
\begin{aligned}
M_0&=N,\quad M_a=N_{a+1}+\cdots+N_\kappa,\quad M_\kappa=0,\\
x^{(0)}&=q=k\ell,\quad x^{(a)}=\lambda^{(a)}\ell\quad(a\ge1),\quad g=c\ell,\\
\theta_w(d)&=2\arctan(d/w).
\end{aligned}
```

The charge and spin equations are

```math
\begin{aligned}
q_j+\sum_b\theta_{g/2}(q_j-x_b^{(1)})&=2\pi I_j,\\
\sum_b\theta_{g/2}(x_j^{(a)}-x_b^{(a-1)})
+\sum_b\theta_{g/2}(x_j^{(a)}-x_b^{(a+1)})
-\sum_{b\ne j}\theta_g(x_j^{(a)}-x_b^{(a)})&=2\pi J_j^{(a)}.
\end{aligned}
```

There is no next-level sum at a=κ-1. Energy depends only on the charge
momenta; auxiliary roots describe the internal component state.

Every sea uses centered consecutive labels $j -(M_{a} -1)/2$, stored as
Uni20 `half_int`. The rational equations require charge-label parity
`2I=M_1 mod 2`, and spin-label parity
`2J^(a)=M_(a-1)+M_(a+1)-M_a+1 mod 2`. Centering the charge sea therefore
requires N_1 odd; centering all the spin seas propagates that requirement
to every occupied N_a. This parity derivation is part of our branch
selection, not a claim that centered labels cover other sectors.

Summing the equations gives `P=2*pi*(sum I+sum_a sum J^(a))/ell` with these
signs. All label sums vanish on the interacting branch, hence P=0. Free
states instead use the sum of their occupied integer modes. Momentum is
signed, with no Brillouin-zone reduction.

## Solving and reporting incomplete continuations

Reflection symmetry leaves `d=sum_a floor(M_a/2)` independent positive
roots. Odd-sized seas have an exactly zero central root. The solver uses
an analytic reduced Jacobian, native-precision dense pivoting, and damped
Newton corrections which preserve the order and positivity of each sea.
Storage is O(d²), with O(d³) linear solves. There is no compiled-in maximum
κ, but this is a dense finite-system tool, not a large-κ thermodynamics code.

Continuation starts at $g =\max (c \,\ell,64N)$ and halves g down to the requested
value. Charge roots start from a centered impenetrable-gas sea; spin roots
use the inverse balanced SU(κ) spin-chain density only as an initial guess.
All finite-size equations are solved, also for imbalanced populations.
Spin coordinates are scaled by `max(1,g)`; changes of scale preserve the
physical roots. At weak coupling, exact integer multiples of pi are
subtracted before the small complementary phases are summed.

Residual normalization is N for g>=1. For g<1 it is
$\max (\sqrt{g },\lvert x_{j} ^{a }\rvert)$ **at every nesting level**. Scaling the spin
equations too is essential when an entire small root cluster shrinks as
sqrt(g). The zero-root equations follow exactly from reflection symmetry;
tests separately check the original rational equations for every root.

`--tolerance` defaults to 32 machine epsilons and is not an energy-error
bound. All calculation and formatting uses the selected fp64, long-double
or optional fp128 arithmetic. Very weak coupling can bring roots near a
nonzero free momentum too close to resolve: use higher precision rather
than treating an incomplete continuation as converged.

`--max-iterations` counts attempted Newton corrections across all stages;
`--max-stages` counts attempted coupling stages. Both default to 10000.
With zero budget no interacting stage is solved, so **no energy or roots
are returned**. If continuation stops after a successful stage, only that
stage's roots and energy are retained; the report explicitly names its
coupling and gives both reached- and requested-coupling residuals. Exit
status 2 means incomplete, 1 means invalid input/numerical exception, and
0 means converged. Exact free cases need no iteration budget.

## Validation and library interface

The tests check the rational equations independently of logarithmic phases,
finite-difference Jacobians, population permutation, unit scaling, honest
budget exhaustion and all three arithmetic types. Physical checks include:

- Reduction to the existing two-component Gaudin–Yang solver.
- One particle per component: the fully antisymmetric spin singlet has
  symmetric spatial ground state, with exactly the Lieb–Liniger energy.
- Free occupations and the first-order shift
  `E-E_free=2c*sum_(a<b) N_a*N_b/ell+O(c²)` in closed-shell sectors.
- A separate plane-wave Fock-space Hamiltonian for populations `3,1,1`.
  Increasing the momentum cutoff gives decreasing variational upper
  energies approaching the Bethe result; a finite cutoff is not exact.
- The balanced SU(3) strong-coupling spin-chain limit:
  `E_infinity=pi²*N*(N²-1)/(3ell²)` and
  $E /E_{\mathrm{infinity}} =1-2\,(N -E_{\mathrm{permutation}})/(c \,\ell)+O ((c \,\ell)^{-2})$.
  The scaled nested spin roots approach those of the permutation chain.
- Weak-to-strong sweeps through six occupied components and 30 particles,
  plus a tiny-c three-particle test that resolves the leading interaction
  energy rather than merely returning a small absolute residual.

```cpp
#include <bethe/su_fermions.hpp>
std::vector<std::size_t> populations{5,3,1};
auto s = bethe::su_fermions::ground_state<long double>(populations, 9.L, 1.L);
if (s.converged) {
    auto energy = *s.energy;
    auto const& momenta = s.rapidities[0];
}
```

`State<Real>` keeps physical populations, the occupied-component nesting
permutation, quantum numbers at all levels, physical rapidities, optional
energy/reached coupling, residuals, counters and status. Before any
successful interacting stage, the energy/coupling optionals and rapidity
array are empty; residual diagnostics are then unavailable (their scalar
storage is zero-initialized). Free modes remain in the original component
order. `ground_quantum_numbers(populations)` is available only for the
supported interacting branch. The model header has no presentation or
CLI dependency, leaving future Python bindings independent of the frontend.
