# Richardson pairing

[Overview](../README.md) | [Model catalogue](models.md) | [Precision and CLI controls](command-line.md)

`bethe-richardson` calculates the lowest energy of the reduced BCS pairing
Hamiltonian in a **specified pair-number and blocked-level sector**:

```math
\begin{aligned}
H&=\sum_i\epsilon_i(n_{i,\uparrow}+n_{i,\downarrow})-g\sum_{i,j}b_i^\dagger b_j,\\
b_i^\dagger&=c_{i,\uparrow}^\dagger c_{i,\downarrow}^\dagger,\qquad g\ge0.
\end{aligned}
```

The input epsilon_i are **single-particle energies**, in strictly increasing
order. Each level is a time-reversed doublet and accommodates at most one
pair. Distinct levels need not be equally spaced. The pairing sum **includes
i=j**: a filled pair contributes -g even when it cannot hop. There is no
chemical potential, half-filling subtraction or pseudospin vacuum offset.

A blocked level contains one fermion, contributes epsilon_i to the energy,
and cannot receive or donate a pair. Its spin orientation does not affect
the energy here. The solver does not optimize which levels are blocked.
For M pairs and B blocked levels the total fermion number is 2M+B.

There is no spatial chain or lattice momentum: PBC/OBC suffixes would have
no meaning. Repeated energies, higher orbital degeneracies, repulsive g<0,
excited states and pair-rapidity output are not implemented in this first slice.
In particular, do not simulate a higher degeneracy by repeating an energy.

## First calculations

```sh
build/bethe-richardson --levels 0,1,2,3 --pairs 2 --g 0.7
build/bethe-richardson --levels -2,-0.8,0.1,0.9,2,3.5 --pairs 2 --g 1 --blocked 1
build/bethe-richardson --levels 0,1 --pairs 1 --g 1 --precision long-double --variables
# With binary128 enabled:
build/bethe-richardson --levels 0,1,2,3 --pairs 2 --g 1 --precision fp128
```

`--levels`, `--pairs` and `--g` are required. Blocked indices are zero-based
positions in the supplied list; their input order is immaterial, but
duplicates or out-of-range indices are errors. The energy list is not sorted
silently, so the meaning of a blocked index cannot change during parsing.

The one-pair, two-level example epsilon=(0,1) has energy
$1-g -\sqrt{1+g \,g }$, hence E=-sqrt(2) at g=1. Empty, completely paired, and
zero-coupling sectors use exact formulas and need no Newton iterations.

## Avoiding pair-root collision singularities

Let e_i=2*epsilon_i for the L unblocked levels. Richardson's pair rapidities
E_alpha obey

```math
\begin{aligned}
\frac1g+\sum_i\frac1{E_\alpha-e_i}
-\sum_{\beta\ne\alpha}\frac2{E_\alpha-E_\beta}&=0,\\
E_{\mathrm{total}}&=\sum_\alpha E_\alpha+\sum_{i\ \mathrm{blocked}}\epsilon_i.
\end{aligned}
```

Real pair rapidities can meet a level pole and become a complex conjugate
pair as g changes. Following those individual roots directly is needlessly
singular. We use [Faribault et al.](../CITATIONS.md#faribault-2011)'s real,
eigenvalue-based variables

```math
\begin{aligned}
y_i&=g\Lambda(e_i)=g\sum_\alpha\frac1{e_i-E_\alpha},\\
F_i&=y_i(y_i-1)-g\sum_{j\ne i}\frac{y_i-y_j}{e_i-e_j}=0,\\
\sum_i y_i&=M.
\end{aligned}
```

At g=0, y_i is 1 on the M lowest unblocked levels and 0 elsewhere. At
nonzero g these variables are **not pair occupation probabilities**, and
the program does not label them as such. They remain regular through the
root collision, using the continuous limiting value where the displayed
root sum has individually divergent terms. `--variables` prints these y_i,
not reconstructed pair rapidities.

Summing the original equations after multiplying by E_alpha gives the
energy without root recovery:

```math
E_{\mathrm{total}}=\sum_i e_i y_i-gM(L-M+1)+\sum_{i\ \mathrm{blocked}}\epsilon_i.
```

The e_i in the referenced pseudospin equations are our pair energies,
not our input epsilon_i. Our unshifted Hamiltonian adds the pseudospin
vacuum constant `sum_(active i) epsilon_i` relative to
`sum_i e_i*S_i^z-g*S^+*S^-`, as well as blocked single-particle energies.
The [Dukelsky–Pittel–Sierra review](../CITATIONS.md#dukelsky-2004) supplies
broader pairing-model and blocking context.

## Ground-state continuation and failures

We continue the lowest zero-coupling occupation pattern towards the requested
g. For g>0 the fixed-pair Hamiltonian has connected negative off-diagonal
pair hopping, so its ground state is unique in a fixed blocked sector.
Following that branch avoids physical ground-state crossings; numerical
branch tracking still needs safeguards.

Each stage uses a tangent predictor, a small change limit in all y_i,
an analytic Jacobian and damped Newton correction. The constraint sum y=M
eliminates y_0, which is occupied at g=0. Eliminating an empty variable
instead can erase its small value by subtraction at weak coupling.
All L equations participate in a rectangular QR correction with L-1
unknowns; no equation is discarded and no normal equations are formed.
Native-precision Givens rotations provide a recoverable rank-failure path,
also for long double without a BLAS/LAPACK backend.

Stages are checked against Fermi-sea and uniform-pair variational upper
bounds, a kinetic-plus-interaction lower bound, and the concavity/tangent
bound along the ground-energy branch. These are necessary safeguards,
not an exhaustive diagonalization or a rigorous interval certificate.
Internal translation and scaling of energies improve numerical conditioning
without changing the declared physical Hamiltonian.

The residual is a polynomial backward error:

```math
r_i=\frac{|F_i|}{1+|y_i(y_i-1)|+
\displaystyle\sum_{j\ne i}\left|\frac g{e_i-e_j}\right|(|y_i|+|y_j|)},
\qquad \mathrm{residual\_norm}=\max_i r_i.
```

The denominator includes the separate polynomial terms before their
cancellation, which matters for closely spaced levels and strong g.
Default tolerance is 32 epsilon of the selected type, **not an energy-error
bound**. Neither tolerance nor arithmetic precision changes silently.
Dense storage is O(L^2), each QR correction costs O(L^3). Very small gaps
or extreme scale ratios can require higher precision or return unresolved.

`--max-iterations` counts attempted Newton corrections, including rejected
stages and failed linear solves. `--max-stages` counts all attempted stages.
Both default to 10000; zero suppresses continuation. Analytic sectors still
return their exact result with a zero budget. If continuation stops, the
result retains the **last completed stage**: its energy belongs to
`reached_coupling`, not the requested coupling. The report shows both,
plus reached and target backward residuals. CLI exits are 0 for reaching
the target, 2 for incomplete continuation, and 1 for invalid input.

## Library and validation

```cpp
#include <bethe/richardson.hpp>
std::vector<long double> levels{0,1,2,3};
std::vector<std::size_t> blocked{1};
auto state = bethe::richardson::ground_state<long double>(levels, 1, 0.7L, blocked);
if (!state.converged) { /* energy is at state.reached_coupling */ }
```

The API retains levels, sorted blocked indices, active indices and their
eigenvalue variables, all counters and diagnostics. The library and CLI
retain fp64, long double and optional fp128 arithmetic and scalar I/O.

Tests compare every pair sector of small irregular and clustered spectra,
with and without blocked levels, against independently constructed pair-space
Hamiltonians. Separate native-precision tests cover the one-pair analytic
energy, weak coupling below epsilon, exact empty/full/free limits, energy
scaling and shifts, the analytic Jacobian, rectangular QR, budgets and
larger systems through 48 levels. For four equally spaced levels with two
pairs, tests reconstruct the pair roots on both sides of the g=2/3 collision
and check the original rational equations. At the collision itself E=0,
both roots coincide at the pole, and the regular variables are exactly
`(7/9,2/3,1/3,2/9)`.

Next extensions are robust optional rapidity reconstruction, excited-state
occupation seeds, and derivative equations for repeated levels/higher
degeneracies. They are separate capabilities, not flags on this first solver.
