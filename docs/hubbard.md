# Periodic Hubbard ground states

[Back to the overview](../README.md)

The Hubbard model adds mobile charge to the spin-chain problem: an electron
can hop between sites, while opposite spins pay an interaction energy when
they occupy the same site. Its nested Bethe ansatz therefore has **two root
families**, charge momenta `k_j` and spin rapidities `Lambda_a`.

`bethe-hubbard-pbc` is a separate front end, backed by `bethe::hubbard` in
`<bethe/hubbard.hpp>`. The first implementation calculates the ground state
of an **even periodic ring at half filling and Sz=0**, for repulsive `U>0`,
with an exact free-fermion path at `U=0`. Doping, polarized sectors, excited
states, attractive U, and open boundaries are not yet implemented.

## First calculations

```sh
build/bethe-hubbard-pbc 6 --u 4
build/bethe-hubbard-pbc 6 --u 4 --roots
build/bethe-hubbard-pbc 32 --u 1 --precision long-double
build/bethe-hubbard-pbc 6 --u 0 --max-iterations 0
# In a binary128-enabled build:
build/bethe-hubbard-pbc 32 --u 4 --precision fp128
```

The six-site ground energy at U=4 is approximately `-3.66870617887296`.
`--particles` and `--sz` are optional explicit sector checks: they currently
must equal L and zero respectively. Unsupported sectors are rejected, not
silently mapped to half filling. `--u` is required; zero must be requested
explicitly. The common [precision and presentation controls](command-line.md)
apply, including solver CPU time and native-precision parameter parsing.

## Hamiltonian and energy convention

We set the hopping amplitude t=1 and use periodic site indices:

```text
H = -sum_(j,sigma) (c^dagger_(j,sigma) c_(j+1,sigma) + h.c.)
    + U sum_j n_(j,up) n_(j,down).
```

There is no chemical potential or magnetic field. L is the number of sites,
N=L the number of electrons, and M=L/2 the number of down-spin electrons.
Thus Sz=(N-2M)/2=0. The energy is the **total unshifted Hubbard energy**,
not that of `U (n_up-1/2)(n_down-1/2)`. At half filling the latter convention
has energy `E_shifted = E - U*L/4`.

For L=2 the periodic sum contains the same hopping bond twice, with effective
two-site hopping amplitude 2. This gives the useful exact check

```text
E(L=2) = (U - sqrt(U^2+64))/2.
```

It is not the isolated dimer with hopping amplitude 1.

## Lieb-Wu equations and ground-state labels

Let `u=U/4`. For U>0, using principal real arctangents, we solve

```text
F_charge(j) = L*k_j - 2*pi*I_j
              + 2*sum_a atan((sin(k_j)-Lambda_a)/u) = 0,

F_spin(a) = 2*sum_j atan((Lambda_a-sin(k_j))/u)
            - 2*sum_b atan((Lambda_a-Lambda_b)/(2*u)) - 2*pi*J_a = 0.

E = -2*sum_j cos(k_j),
P = sum_j k_j = (2*pi/L)*(sum_j I_j + sum_a J_a)  (mod 2*pi).
```

These are the logarithmic Lieb-Wu equations; see the
[reference and equation provenance](../CITATIONS.md). The spin self-scattering
term is zero. `Lambda` in the API and root tables is the conventional Lieb-Wu
rapidity, not a rescaled quantity.

The two even-length classes have different quantum-number parity:

| Ring length | Charge labels I (L values) | Spin labels J (M values) | Momentum |
| --- | --- | --- | --- |
| L=4m+2 | `-(L-1)/2, ..., (L-1)/2`, half-odd integers | `-(M-1)/2, ..., (M-1)/2`, integers | 0 |
| L=4m | `-L/2+1, ..., L/2`, integers | `-(M-1)/2, ..., (M-1)/2`, half-odd integers | pi |

Each list advances in steps of one. Charge labels span a full Brillouin zone;
for L=4m we choose the endpoint +pi rather than -pi. The API stores labels
exactly using Uni20's `half_int`, and momentum as both an exact integer index
and `2*pi*index/L`. The parity conventions and momenta are tested against
small-ring exact diagonalization, including the fermionic sign of translation.

## How the numerical solve works

We exploit reflection symmetry, solving only for positive k and Lambda.
The partners are reflected exactly. For L=4m the roots k=0 and k=pi are fixed;
for L=4m+2 the central spin rapidity is exactly zero. In particular, the
fixed k=pi contributes `sin(k)=0` exactly, avoiding division of a floating
approximation to sin(pi) by a very small U.

The reduced system uses damped Newton updates with an analytic Jacobian.
Uni20 performs the dense linear solve in the selected scalar type; there is
no conversion through double. Spin variables are internally scaled as
`Lambda/max(1,u)` to condition the large-U equations, and converted back on
return. For U<8 we first solve at U=8 and halve U stage by stage down to the
requested value, preserving physical Lambda between stages. For U>=8 we
solve directly at the requested interaction.

The reduced dense solve has approximately 3L/4 unknowns: O(L^2) storage and
O(L^3) work per Newton update, plus O(L^2) equation/Jacobian evaluation. This
is a starting point for moderate rings, not yet a matrix-free large-L solver.

Diagnostics include separate `max|F_charge|/L` and `max|F_spin|/L`, their
maximum, accepted Newton updates, and the number of completed continuation
stages. The default tolerance is 32 times the selected type's epsilon.
The iteration budget applies to **all stages together**, not to each stage.
A zero budget evaluates the initial large-U seed without taking a Newton
step; the exact U=0 calculation requires no updates.

If the budget is exhausted or line search cannot improve the residual, the
returned state is explicitly marked unconverged (CLI exit status 2).
Even if continuation stops early, the returned roots, energy, and residuals
are evaluated at the **requested U**, not presented as a solved intermediate
interaction. The momentum field identifies the intended ground-state branch;
an unconverged energy must not be treated as a certified eigenvalue.

Weak coupling brings nearby roots close together, so a solve can stall at a
precision-dependent limit. Try higher precision, or a larger update budget
if the budget was exhausted. At extremely strong coupling the energy sums
terms of order one to obtain a result of order 1/U; cancellation limits
relative energy accuracy even if the residual is small. A small equation
residual is **not an energy-error bound**.

## Free fermions and strong coupling

At U=0, filling the lowest L/2 one-particle levels with each spin gives

```text
E = -4*cot(pi/L)  for L=4m,
E = -4*csc(pi/L)  for L=4m+2.
```

We evaluate these expressions directly. The interacting Bethe coordinates
are singular at U=0, so `charge_momenta` instead holds the occupied
free-fermion momenta, with repetitions for opposite spins.
`spin_rapidities` and both Bethe-label lists are empty, and `free_fermion`
is true. For L=4m the free ground state is degenerate; we choose a
representative with momentum pi. This is an occupation-state representative,
not a claim to reproduce the interacting singlet's limiting wavefunction.

At large U and L>=4, virtual hopping generates an antiferromagnetic XXX chain
with exchange J=4/U. In the conventions of this repository,

```text
(U/4)*E_Hubbard -> E_XXX(J=1) - L/4.
```

The offset matters: dropping it would compare different Hamiltonians.
The regression tests check this limit against the existing XXX solver, the
free-fermion limit against independently filled one-particle levels, and
finite-U energies and momentum against exact fermionic diagonalization
through six sites. Native-precision two-site tests distinguish long-double
and fp128 calculations from fp64 results padded with extra output digits.

## C++ API

```cpp
#include <bethe/hubbard.hpp>

auto state = bethe::hubbard::ground_state<long double>(32, 4.0L);
if (!state.converged) {
    // Inspect state.status, residual_norm, and iterations before using energy.
}
auto const& k = state.charge_momenta;
auto const& lambda = state.spin_rapidities;
auto labels = bethe::hubbard::ground_quantum_numbers(32); // U>0 branch
```

`State<Real>` carries the two root families, separate quantum-number lists,
energy, momentum, residual diagnostics, requested interaction, and particle
counts. `SolverOptions<Real>` supplies `residual_tolerance` and
`max_iterations`, as for the other models. Link the public `bethe::bethe`
CMake target, which includes the Uni20 linear-algebra dependency.

The nested state type intentionally differs from the spin-chain `RealState`:
future doping, spin sectors, and charge/spin excitations should extend this
module without overloading the meaning of a single rapidity vector. There
is no arbitrary-quantum-number solve or Hubbard excitation scan yet; a
complete spectrum would also require complex spin and charge-spin strings.
