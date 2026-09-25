# Repulsive Gaudin–Yang gas on a ring

[Overview](../README.md) · [Command-line controls](command-line.md) · [Model catalogue](models.md)

This is the spin-1/2 fermionic counterpart of the [Lieb–Liniger gas](lieb-liniger.md).
The extra spin degree of freedom requires a second, auxiliary set of Bethe
roots. Charge momenta determine the energy; spin rapidities determine how
the two spin populations scatter.
For more components, [SU(n) fermions](su-fermions.md) supplies a separate
population-list API/frontend with arbitrary nesting depth and the analogous
odd-population restriction. The two-component frontend here remains available.

## Units and first calculations

We use equal masses, a ring of physical circumference $`\ell`$, and

```math
\begin{aligned}
H&=-\sum_j\frac{\partial^2}{\partial x_j^2}+2c\sum_{i<j}\delta(x_i-x_j),\\{}
\frac{\hbar^2}{2m}&=1,\qquad E=\sum_jk_j^2,\qquad P=\sum_jk_j.
\end{aligned}
```

The particle count N is **not** a lattice length. The coupling c has units
of inverse length, momentum has units of inverse length, and energy has
units of inverse length squared. Same-spin contact interactions vanish by
fermionic antisymmetry. The original solutions are
[Gaudin (1967)](../CITATIONS.md#gaudin-1967) and
[Yang (1967)](../CITATIONS.md#yang-1967); our explicit finite-ring equations
follow [Oelkers et al.](../CITATIONS.md#oelkers-2006), Eqs. (1), (2), (7), (8), (25).

```sh
build/bethe-gaudin-yang-pbc 6 --length 6 --c 1
build/bethe-gaudin-yang-pbc 8 --length 8 --c 1 --sz 1 --roots
build/bethe-gaudin-yang-pbc 18 --length 18 --c 0.001 --precision long-double
# Requires a binary128-enabled build:
build/bethe-gaudin-yang-pbc 6 --length 6 --c 1 --precision fp128
```

## Supported sectors

$`S^z =(N_{\mathrm{up}} -N_{\mathrm{down}})/2`$. The default is 0 for even N and 1/2 for odd N;
`--sz` accepts Uni20's exact integer/half-integer notation, including negative
values. Invalid population parity or $`\lvert 2S^z \rvert \gt N`$ is rejected.

The first interacting implementation requires **both N_up and N_down odd**.
Equivalently, N is even and the minority count M is odd. This includes
balanced N=2,6,10,... and polarized sectors such as N=8, Sz=1. The centered
filled-sea branch is the one selected in Sec. 5 of Oelkers et al.; the same
sector restriction is used in their weak/strong-coupling analysis.
Other mixed-spin periodic shell branches are **not implemented**. In
particular, N=4, Sz=0, c>0 is rejected. A balanced default does not imply
that every N has an implemented interacting ground state.

Spin reversal leaves the energy unchanged. The solver chooses the majority
spin as reference vacuum, and M spin roots describe the minority particles.
The report identifies this reference; it does not relabel physical N_up/down.

Two exact exceptions allow arbitrary nonnegative populations:

- **c=0:** two independent filled free-fermion seas.
- **Full polarization:** a single free-fermion sea at any finite c>=0.

These include N=0. For a species containing K particles, the occupied
integer modes are `j-floor((K-1)/2)`, `j=0,...,K-1`. If K is even this chooses
the positive-current representative of the degenerate free ground state;
the reflected occupation has the same energy. Consequently the free path
can report nonzero P. It returns **occupations, not interacting Bethe labels**,
and has no auxiliary spin roots. `--roots` prints a separate mode table for
each physical spin component.

Attraction, excitations, twists, hard walls, thermodynamic integral equations,
and unimplemented shell branches are not inferred from the cited literature.
In particular, the Hubbard Shiba mapping does not supply a continuum c<0 solver.

## Two nested equations

Let $`q =k \,\ell`$, $`l =\lambda \,\ell`$, and $`g =c \,\ell`$. The implementation solves

```math
\begin{aligned}
q_j+\sum_a2\arctan\!\left(\frac{2(q_j-l_a)}g\right)&=2\pi I_j,\\{}
\sum_j2\arctan\!\left(\frac{2(l_a-q_j)}g\right)
-\sum_{b\ne a}2\arctan\!\left(\frac{l_a-l_b}g\right)&=2\pi J_a,\\{}
I_j&=j-\frac{N-1}{2},\quad j=0,\ldots,N-1,\\{}
J_a&=a-\frac{M-1}{2},\quad a=0,\ldots,M-1.
\end{aligned}
```

Here the I labels are half-odd integers and the J labels are integers.
They are stored as `uni20::half_int`. Physical k and lambda, not q and l,
are returned and displayed. Unlike the XXX tool's z coordinate, lambda
has no additional factor of two.

Summing these equations cancels pairwise spin phases and gives
`P=2*pi*(sum I + sum J)/ell` with **these signs**. Both sums vanish in the
supported interacting branch. This is continuum momentum: no modulo-2pi
or Brillouin-zone reduction is made. The free path instead obtains the
integer momentum index from the sum of its occupied integer modes.

## Numerical method and diagnostics

The numerical choices below are our implementation, not a promise of the
literature's full coverage. Reflection symmetry leaves N/2 positive charge
roots and (M-1)/2 positive spin roots; the central spin root is exactly zero.
The solver builds an analytic reduced Jacobian and uses Uni20's dense linear
solve with backtracking to preserve ordered positive roots. For reduced
order $`d =(N +M -1)/2`$, storage is O(d²) and a dense Newton solve is O(d³).

Continuation begins at $`g =\max (c \,\ell,64N)`$ and halves g towards the requested
value after each converged stage. Spin coordinates are stored as
$`l /\max (1,g)`$ to accommodate their large-c growth. Changing that scale preserves
the physical spin roots; it does not move them towards zero at weak coupling.
At g<1 the code subtracts integer multiples of pi analytically before summing
complementary phases, resolving the central charge pair of size sqrt(g).

The residual is the maximum over the independent positive-root equations:

- For g>=1, both families' absolute residuals are divided by N.
- For g<1, each charge residual is divided by $`\max (\sqrt{g },\lvert q_{j} \rvert)`$;
  spin residuals remain divided by N.

Thus residuals on opposite sides of g=1 have different normalizations.
`--tolerance` defaults to 32 times the selected type's epsilon. It is **not an
energy-error bound**. Charge, spin, and maximum residuals are reported separately.

`--max-iterations` counts accepted Newton updates across **all** continuation
stages, not independently at each coupling. On exhaustion or a stalled line
search, the tool returns exit status 2 and labels the result unconverged.
Even then, residuals are recomputed at the **requested** c; `Root coupling
reached` records the intermediate continuation coupling. An unfinished root
set must not be used as an eigenstate at the requested interaction.

fp64, platform-native long double, and optional fp128 run the same equations
without narrowing parameters or output. Finite positive ell, finite c>=0,
and a finite positive tolerance are required; interacting $`c \,\ell /2`$ must be
finite and nonzero. Particle-label and dense-matrix allocation sizes are
checked before allocation. Nonfinite physical roots/observables and complete
underflow of a nonzero root or individual kinetic term raise errors: change
units or precision rather than interpreting a rounded zero as a physical result.
These checks do not guarantee relative accuracy in the subnormal range.

## Validation and useful limits

The precision-typed regression suite checks the original complex rational
equations independently of the solver's logarithmic residual and compares
its analytic Jacobian with central finite differences. It also checks:

- Exact free-fermion occupation energies, polarization, vacuum, and spin reversal.
- The two-body contact jump condition: at ell=1 the positive root obeys
  $`q \,\tan (q /2)=c /2`$. At c=pi, $`q =\pi /2`$ and `E=pi²/2`, providing an
  independent native-precision oracle. The two-body energies also agree
  with the bosonic contact problem.
- The first-order weak-coupling shift `E-E_free = 2c*N_up*N_down/ell + O(c²)`.
- The strong-coupling charge energy `E_infinity=pi²*N*(N²-1)/(3ell²)` and
  its leading correction `E/E_infinity=1-2*Omega/(c*ell)+O(c^-2)`, where
  `Omega=N/2-2E_XXX` for the N-site XXX sector with M overturned spins.
  Spin roots satisfy $`\lambda /c \to z_{\mathrm{XXX}} /2`$. The weak/strong expansions are
  checked against Eqs. (12), (15)-(16) of Oelkers et al.
- A controlled dilute-Hubbard discretization: with lattice spacing a,
  $`U_{\mathrm{lattice}} =2c \,a`$ and hopping set to 1, the continuum energy is approached
  by `(E_Hubbard+2N)/a²`. Doubling the lattice resolution checks convergence;
  this is not a claim that the two finite models have identical energies.
- Unit rescaling, continuation budgets, weak-to-strong sweeps through N=32,
  invalid inputs, and precision-preserving command-line output.

The strong-coupling E_infinity is the limit of this interacting branch,
not the canonical even-N fully polarized free ground state with its different
periodic shell occupation. No small residual alone establishes state selection
or completeness of an excitation spectrum.

## Library interface

```cpp
#include <bethe/gaudin_yang.hpp>

auto state = bethe::gaudin_yang::ground_state<long double>(
    5, 3, 8.0L, 1.0L); // N_up, N_down, physical length, c
if (!state.converged) {
    // Handle an incomplete continuation; do not treat energy as converged.
}
auto labels = bethe::gaudin_yang::ground_quantum_numbers(5, 3);
```

`State<Real>` owns physical populations, length, c, charge momenta, spin
rapidities, exact labels, observables, both residuals, status, and iteration
count. Its `free_modes` are populated only on the exact free path; conversely
the interacting quantum-number arrays are empty on that path.
`ground_quantum_numbers` describes only the supported interacting branch.
The header has no CLI or presentation dependency and is suitable for future
Python bindings without changing the numerical API.
