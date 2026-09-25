# Repulsive Lieb–Liniger bosons on a ring

[Overview](../README.md) · [Build guide](building.md) · [Model catalogue](models.md)

For Dirichlet boundaries, use [bethe-lieb-liniger-obc](lieb-liniger-open.md).

This is the continuum counterpart of a chain calculation: specify a physical
ring circumference $`\ell`$, a particle number $`N`$, and a repulsive coupling $`c`$.
There is no lattice spacing or site count. The executable
`bethe-lieb-liniger-pbc` supports the ground state, explicit Bethe labels, and
finite-window excitation scans, using fp64, long double, or optional fp128.

## Start with a calculation

```sh
build/bethe-lieb-liniger-pbc 4 --length 4 --c 1 --roots
build/bethe-lieb-liniger-pbc 4 --length 4 --c 1 --precision long-double
build/bethe-lieb-liniger-pbc 4 --length 4 --c 1 --excitations all --padding 1
```

Both `--length` and `--c` are required; neither is silently chosen from $`N`$.
The first command has density $`n =N /\ell =1`$ and dimensionless interaction
$`\gamma =c /n =1`$. The last solves 15 label configurations, including the ground
state. Add `--roots` to see the physical momenta and labels for each retained
state, or `--format plain` for uncolored key/value metadata and numeric tables.
The report includes solver CPU time and convergence diagnostics.

## Hamiltonian and labels

We use units $`\hbar ^{2}/(2m)=1`$ and

```math
\begin{aligned}
H&=-\sum_j\frac{\partial^2}{\partial x_j^2}+2c\sum_{i<j}\delta(x_i-x_j),\qquad c>0,\\{}
E&=\sum_j k_j^2,\qquad P=\sum_j k_j.
\end{aligned}
```

Thus $`c`$ has inverse-length units and energy has inverse-length-squared units.
The finite-ring equations are

```math
\ell k_j+2\sum_{l\ne j}\arctan\!\left(\frac{k_j-k_l}{c}\right)=2\pi I_j.
```

Labels are distinct and increasing: integers for odd $`N`$, half-odd integers
for even $`N`$. Ground-state labels are $`I_{j} =j -(N -1)/2`$, with zero-based $`j`$.
These conventions match the explicit equations (4), (24)–(27) in
[Essler–de Klerk](https://arxiv.org/html/2307.12410v1#S2.SS2), following
the original [Lieb–Liniger solution](../CITATIONS.md#lieb-liniger-1963).
Repulsive solutions have real roots; unlike a restricted spin-chain real-root
family, no complex-string sector is needed for this model at positive $`c`$.

For example, move the highest ground label up one slot:

```sh
build/bethe-lieb-liniger-pbc 4 --length 4 --c 1 \
  --quantum-numbers '-3/2,-1/2,1/2,5/2' --roots
```

Here `Q=sum(I)=1` and $`P =2\,\pi /\ell`$. **Momentum is not reduced modulo anything**:
a continuum ring has no Brillouin zone. Adding an integer $`b`$ to every label
boosts each momentum by $`v =2\,\pi \,b /\ell`$; the energy becomes
`E+2vP+Nv^2`. The vacuum ($`N =0`$, optionally `--quantum-numbers none`) and
single-particle states are included. The latter have $`k =2\,\pi \,I /\ell`$.

## What an excitation window means

There are infinitely many levels even at fixed $`N`$. A scan therefore requires
`--padding P`: append $`P`$ slots at **each** end of the ground-state label interval
and select every subset of $`N`$ labels. There are `binomial(N+2P,N)` candidates
(one empty configuration for $`N =0`$). Any number of particle/hole rearrangements
within that window is included, not just a single particle/hole pair.

`--excitations COUNT` retains the lowest `COUNT` converged levels; `all` retains
all converged candidates. Both solve the entire window. `--padding 0` means
ground state only. Distinct degenerate states are retained, and gaps are
relative to the independently converged global ground state. Exact numerical
energy ties use lexicographic labels; rounding can change the order of a
degenerate pair, not its membership.

`--max-candidates` defaults to 10000. The combinatorial count is checked before
any state is solved or roots are allocated; oversized scans fail with a useful
error. Increase this guard explicitly if the cost is acceptable. The retained
states use a bounded heap, so a small `COUNT` saves storage but not solve time.
A window does **not** certify the globally lowest `COUNT` states; enlarge it
and compare. This is also not a thermodynamic type-I/type-II dispersion API.

Failed solves are excluded from ranking. Their number and the first failed
state's diagnostics are reported, even if all requested retained levels were
found. If the ground reference failed, gaps are unavailable. Check exit status:
0 means convergence (within the window for scans), 1 means invalid inputs or
an unsupported numerical range, and 2 means an incomplete solve/scan.

## Numerical method and precision

The solver works with $`q =k \,\ell`$ and $`g =c \,\ell`$, and uses an analytic dense
Jacobian with Uni20's linear solver. Damped Newton steps preserve strict root
ordering. One step costs O(N^3) work and O(N^2) storage; scan cost multiplies
this by the number of candidates and Newton updates. There is no new Bethe
linear-algebra implementation or BLAS requirement for native long double.

At small $`g`$, direct subtraction of pi-sized phases can hide the small roots.
For $`g \lt 1`$ we instead use the algebraically equivalent ordered-root equation

```math
\begin{aligned}
q_j-2\pi J_j-2\sum_{l\ne j}\arctan\!\left(\frac g{q_j-q_l}\right)&=0,\\{}
J_j&=I_j-\left[j-\frac{N-1}{2}\right].
\end{aligned}
```

The integer subtraction in $`J_{j}`$ is performed before conversion to `Real`.
Phase sums are compensated, and kernels avoid forming unsafe squared scales.
The weak-coupling seed separates roots by order `sqrt(g)` around the free-boson
momenta $`2\,\pi \,J_{j}`$; it crosses over to $`2\,\pi \,I_{j}`$ at large $`g`$.

The reported residual is $`\max_{j} \lvert F_{j} \rvert /s_{j}`$, with
`s_j=max(min(1,sqrt(g)), |q_j|, |target_j|)` and `target_j=2*pi*J_j`
for $`g \lt 1`$, otherwise $`2\,\pi \,I_{j}`$. The vacuum's norm is zero. This measures the
equations at the **returned physical momenta**, not just the internal iterate.
The default tolerance is `32*epsilon(Real)`; `--max-iterations 0` evaluates
only the initial guess. A stalled line search or exhausted budget is reported
as an unconverged estimate. The residual is not an energy-error bound, and a
small absolute equation error need not resolve a tiny gap or a boosted cluster.

All parameters are parsed directly into the selected arithmetic. Higher
precision also applies to roots, energy sums, gap subtraction, and output.
Finite positive $`\ell`$, $`c`$, and representable nonzero $`g =c \,\ell`$ are required.
Very large labels, very close roots, or extreme physical units can exceed the
selected type's resolution or exponent range; use better-scaled units or higher
precision. Nonfinite energies and total-energy underflow are rejected rather
than reported as converged finite answers. Label sums are checked against the
signed 64-bit momentum-index range. Convergence is not guaranteed for every
input representable by the interface.

## C++ API

```cpp
#include <bethe/lieb_liniger.hpp>
namespace ll = bethe::lieb_liniger;

auto ground = ll::ground_state(4, 4.0L, 1.0L);
auto numbers = ll::ground_quantum_numbers(4);
numbers.back() += uni20::half_int{1};
auto excited = ll::solve_real(4.0L, 1.0L, numbers);
auto scan = ll::real_excitations(4, 4.0L, 1.0L, 1, {.count = 10});
// Check ground.converged, excited.converged, and scan.converged().
```

`State<Real>` carries physical `momenta`, `energy`, signed `momentum`, exact
`momentum_index`, labels, iterations, residual, and `SolveStatus`. An optional
`span<Real const>` of physical initial momenta may be supplied to `solve_real`
after its solver options. The label count determines $`N`$; an initial guess
must have matching size and finite, strictly increasing roots.

## Validation and remaining scope

[Tests](../tests/test_lieb_liniger.cpp) include an independent two-body
boundary-condition solve, analytic two-/three-body points, weak-coupling
`E ~ c*N*(N-1)/ell`, the Tonks limit
$`E \to \pi ^{2}\,N \,(N ^{2}-1)/(3\,\ell ^{2})`$, boosts, reflection, length scaling, Jacobian
finite differences, independently enumerated windows, and failure paths.
Native-precision tests distinguish extended arithmetic from fp64 narrowing.
A separate quadrature solution of the bulk root-density integral equation
checks finite-size convergence at fixed density; it is test code, not a public
thermodynamics solver. A separate [public bulk ground-state API](lieb-liniger-thermo.md)
now provides the Fermi rapidity, energy per length, and chemical potential.

Attractive interactions/bound states, exactly zero or infinite
coupling, finite-temperature TBA, and matrix
elements are not implemented. The separate [hard-wall frontend](lieb-liniger-open.md)
implements reflection equations, not a boundary toggle on the ring frontend;
negative $`c`$ must not be treated as a sign switch in this real-root solver.
For thermodynamic type-I/type-II curves, use the separate
[dispersion frontend](lieb-liniger-thermo.md), not a finite label scan.
See the [catalogue](models.md#lieb-liniger-the-simplest-new-interacting-family)
and [bibliography](../CITATIONS.md) for the next extensions and references.
