# Free-end Hubbard ground states

[Back to the overview](../README.md) · [Periodic Hubbard](hubbard.md)

`bethe-hubbard-obc` solves the Hubbard chain with free ends: electrons hop
between neighboring sites, but there is no bond closing the chain and no
boundary potential or magnetic field. It supports **every physical particle
number and spin projection**, either sign of U, and odd or even lengths L>=1.
It is a separate front end backed by `bethe::hubbard::open` in
`<bethe/hubbard_open.hpp>`. Hubbard excitations are not implemented yet.

## First calculations

```sh
build/bethe-hubbard-obc 6 --u 4
build/bethe-hubbard-obc 5 --u 4 --roots
build/bethe-hubbard-obc 16 --u 4 --particles 8
build/bethe-hubbard-obc 15 --u -4 --particles 7 --sz -1/2 --roots
build/bethe-hubbard-obc 32 --u 4 --precision long-double
build/bethe-hubbard-obc 6 --u 0 --max-iterations 0
# In a binary128-enabled build:
build/bethe-hubbard-obc 32 --u 4 --precision fp128
```

`--u` is required. `--particles N` defaults to L; `--sz` defaults to 0 for
even N and 1/2 for odd N. Explicit spin projections accept fractions such as
`-3/2` or decimals such as `1.5`. The physical constraints are simply

```text
N_up = N/2 + Sz,   N_down = N/2 - Sz,
0 <= N_up,N_down <= L, with both counts integral.
```

Unlike periodic rings, open chains have no shell-parity exclusions in the
interacting ground-state interface. Nor is there a conserved lattice momentum:
the returned charge wave numbers describe standing waves, not translation
eigenvalues. All three precision modes, CPU timing, plain/pretty output, and
convergence controls follow the [common CLI guide](command-line.md).

## Hamiltonian and simple checks

With hopping t=1 and sites j=1,...,L, we use

```text
H = -sum_(j=1..L-1,sigma) (c^dagger_(j,sigma) c_(j+1,sigma) + h.c.)
    + U sum_(j=1..L) n_(j,up) n_(j,down).
```

The reported energy is total and unshifted. To convert to the convention
`U*(n_up-1/2)*(n_down-1/2)`, subtract `U*N/2` and add `U*L/4`.
For L=1 there is no hopping: empty or singly occupied states have energy zero,
and the doubly occupied state has energy U.

L=2 is the ordinary **single-bond dimer**. Its two-electron Sz=0 ground energy is

```text
E = (U - sqrt(U^2+16))/2.
```

For example, U=4 gives E approximately `-0.8284271247461901`. This differs
from the periodic L=2 convention, whose site sum counts that bond twice.

At U=0, fill the lowest standing-wave levels independently for each spin:

```text
k_j = pi*j/(L+1),   epsilon_j = -2*cos(k_j),   j=1,...,L.
```

This path is exact up to arithmetic rounding and needs no Newton updates.
The returned `charge_momenta` holds occupied k values, repeated when both
spins occupy a level; Bethe-label lists and `spin_rapidities` are empty.
`free_fermion` also identifies a single-species auxiliary root sector at
nonzero U. Empty and fully filled one-spin bands have exactly zero kinetic
energy in the implementation.

## Reflected scattering and ground-state labels

First consider a repulsive sector, U>0, N<=L, and M=N_down<=N/2.
Set `u=U/4`. All charge roots obey `0<k_j<pi`, all spin rapidities are
positive, and each family is strictly increasing. We solve

```text
F_charge(j) = 2*(L+1)*k_j - 2*pi*I_j
              + 2*sum_a [atan((sin(k_j)-Lambda_a)/u)
                         + atan((sin(k_j)+Lambda_a)/u)] = 0,

F_spin(a) = 2*sum_j [atan((Lambda_a-sin(k_j))/u)
                    + atan((Lambda_a+sin(k_j))/u)]
            - 2*sum_(b!=a) [atan((Lambda_a-Lambda_b)/(2*u))
                            + atan((Lambda_a+Lambda_b)/(2*u))]
            - 2*pi*J_a = 0,

I_j = j,  j=1,...,N;    J_a = a,  a=1,...,M,
E_roots = -2*sum_j cos(k_j).
```

The plus-rapidity terms describe scattering off reflected partners. Both
spin self-scattering factors are excluded: retaining the reflected b=a
term would change the boundary problem. The charge length L+1 is also
essential; it reproduces the free standing-wave spectrum.
See [CITATIONS.md](../CITATIONS.md) for the boundary equations and provenance.

These consecutive integer labels select the sector ground state. For the
open nearest-neighbor chain, spin ordering places the lowest state at the
smallest total spin compatible with Sz; the periodic shell-parity issue
does not arise. Labels use Uni20's exact `half_int`, despite being integers
for this boundary condition. `Lambda` in returned states is the conventional
Hubbard rapidity, not the internally scaled Newton variable.

## Doping, attraction, and auxiliary roots

Every open chain is bipartite, including odd L. The same exact transformations
as in the [periodic sector guide](hubbard-sectors.md#why-attraction-can-use-a-repulsive-solver)
therefore reduce every physical request to the repulsive sector above:

1. For U<0, transform down-spin particles to holes (Shiba mapping).
2. If the transformed particle count exceeds L, transform both species to holes.
3. If necessary, reverse spins so N_up>=N_down.

For g=|U|, their unshifted-energy relations are

```text
E_-g(N_up,N_down) = E_+g(N_up,L-N_down) - g*N_up,
E_U(N_up,N_down)  = E_U(L-N_up,L-N_down) + U*(N-L).
```

There is no approximation to attractive strings here: the **energy** is
obtained by an exact symmetry, while returned roots remain those of the
auxiliary repulsive problem. We do not reconstruct physical attractive
roots or wavefunctions. The CLI labels auxiliary root tables and reports
the root sector and energy offset explicitly.

`particles`, `down_spins`, `interaction`, and `energy` describe the requested
physical sector. `root_particles`, `root_down_spins`, `root_interaction`,
the roots, quantum numbers, and residuals describe the auxiliary sector.
Always interpret them together. `energy = E_roots + energy_offset`;
mapping flags and `auxiliary_roots()` record which transformations were used.
No momentum or momentum-offset member exists in the open-chain state type.

## Numerical method and validation

The open-chain equations and analytic Jacobian are boundary-specific; the
damped Newton driver, interaction continuation, symmetry mappings, and
native-precision Uni20 linear solve are shared with the periodic solver.
There are N+M positive-root unknowns, requiring O((N+M)^2) dense storage and
O((N+M)^3) work per update. This is a moderate-chain solver, not a matrix-free
large-system implementation.

Internally we solve for `Lambda/max(1,u)`. For U<8, continuation starts at
U=8 and halves U down to the target; for U>=8 it starts directly at U.
The reported charge and spin residuals are `max|F|/[2*(L+1)]`, with default
tolerance 32 times the selected precision's epsilon. The update budget is
shared by all stages. Even a stopped solve reports residuals at the requested
root-sector interaction, not an intermediate continuation value.

An exhausted budget or stalled line search returns an unconverged estimate
(CLI exit status 2). Weak-coupling root separations can become unresolved;
strong-coupling kinetic energies suffer cancellation; attractive energy
offsets can hide small binding scales. Higher precision can help, but a small
residual is **not an energy-error bound**.

At half filling and large positive U, virtual hopping gives the open XXX chain:

```text
(U/4)*E_Hubbard -> E_XXX_open(J=1) - (L-1)/4.
```

Tests compare this limit with the existing XXX solver, independently fill
free standing waves, check the analytic Jacobian by finite differences, and
compare every spin population through six sites against fermionic exact
diagonalization for both signs of U. Dimer tests check native arithmetic in
fp64, long-double, and fp128, not merely additional printed digits.

## C++ API

```cpp
#include <bethe/hubbard_open.hpp>

// Half filling: Sz=0 for even L, Sz=1/2 for odd L.
auto state = bethe::hubbard::open::ground_state<long double>(15, 4.0L);
// An arbitrary doped, attractive physical sector.
auto paired = bethe::hubbard::open::sector_ground_state<long double>(
    15, 8, uni20::half_int{0}, -4.0L);
if (!paired.converged) {
    // Inspect status, residual_norm, and iterations before using energy.
}
auto const& k = paired.charge_momenta; // auxiliary standing-wave roots
auto const& lambda = paired.spin_rapidities;
auto labels = bethe::hubbard::open::ground_quantum_numbers(15, 8, 4);
```

`ground_quantum_numbers(L,N,M)` refers to a normalized repulsive root sector,
not an arbitrary physical sector before symmetry mapping. The one-argument
overload uses N=L and M=floor(L/2). Link `bethe::bethe` as for the other
modules. Excited states, arbitrary quantum-number solves, boundary potentials,
and boundary magnetic fields remain future extensions.
