# Periodic Hubbard ground states

[Back to the overview](../README.md)

For infinite-chain spinon and charge lines at or below half filling, use the
separate [thermodynamic dispersion tool](hubbard-dispersion.md), including
[doped dispersions](hubbard-doped.md), both interaction conventions and
Hamiltonian/Fermi energy references.

The Hubbard model adds mobile charge to the spin-chain problem: an electron
can hop between sites, while opposite spins pay an interaction energy when
they occupy the same site. Its nested Bethe ansatz therefore has **two root
families**, charge momenta $k_{j}$ and spin rapidities $\Lambda_{a}$.

`bethe-hubbard-pbc` is a separate front end, backed by `bethe::hubbard` in
`<bethe/hubbard.hpp>`. It calculates sector ground states on **even periodic
rings**, including repulsive half-filled spin sectors, balanced attractive
systems at every even particle number, and selected doped sectors. Not every
finite-ring shell parity is supported: see [sectors and symmetry mappings](hubbard-sectors.md).
There is an unrestricted exact free-fermion path at $U =0$. Excited-state scans
are not yet implemented. For free ends, use the separate
[`bethe-hubbard-obc` front end](hubbard-open.md), which supports every physical
particle/spin sector and odd as well as even lengths.

## First calculations

```sh
build/bethe-hubbard-pbc 6 --u 4
build/bethe-hubbard-pbc 6 --u 4 --roots
build/bethe-hubbard-pbc 32 --u 1 --precision long-double
build/bethe-hubbard-pbc 6 --u 0 --max-iterations 0
build/bethe-hubbard-pbc 16 --u 4 --particles 6
build/bethe-hubbard-pbc 16 --u -4 --particles 8 --roots
build/bethe-hubbard-pbc 16 --u 4 --sz 2
# In a binary128-enabled build:
build/bethe-hubbard-pbc 32 --u 4 --precision fp128
```

The six-site ground energy at U=4 is approximately `-3.66870617887296`.
`--particles` and `--sz` select the physical sector, defaulting to L and zero.
Unsupported sectors are rejected. `--u` is required; zero must be requested
explicitly. Symmetry mappings preserve the requested physics and report any
auxiliary root sector separately. The common [precision and presentation controls](command-line.md)
apply, including solver CPU time and native-precision parameter parsing.

## Hamiltonian and energy convention

We set the hopping amplitude t=1 and use periodic site indices:

```math
H=-\sum_{j,\sigma}\left(c^\dagger_{j,\sigma}c_{j+1,\sigma}+\mathrm{h.c.}\right)
+U\sum_j n_{j,\uparrow}n_{j,\downarrow}.
```

There is no chemical potential or magnetic field. L is the number of sites,
N the number of electrons, and M the number of down-spin electrons.
Thus Sz=(N-2M)/2. The energy is the **total unshifted Hubbard energy**,
not that of $U (n_{\mathrm{up}} -1/2)(n_{\mathrm{down}} -1/2)$. At half filling the latter convention
has energy $E_{\mathrm{shifted}} = E - U \,L /4$; in general it is $E - U \,N /2 + U \,L /4$.

For L=2 the periodic sum contains the same hopping bond twice, with effective
two-site hopping amplitude 2. This gives the useful exact check

```math
E(L=2)=\frac{U-\sqrt{U^2+64}}{2}.
```

It is not the isolated dimer with hopping amplitude 1.

## Lieb-Wu equations and ground-state labels

Let $u =U /4$. In the repulsive root sector, using principal real arctangents, we solve

```math
\begin{aligned}
F_{\mathrm{charge}}(j)&=Lk_j-2\pi I_j
+2\sum_a\arctan\!\left(\frac{\sin k_j-\Lambda_a}{u}\right)=0,\\
F_{\mathrm{spin}}(a)&=2\sum_j\arctan\!\left(\frac{\Lambda_a-\sin k_j}{u}\right)
-2\sum_b\arctan\!\left(\frac{\Lambda_a-\Lambda_b}{2u}\right)-2\pi J_a=0,\\
E_{\mathrm{roots}}&=-2\sum_j\cos k_j,\\
P_{\mathrm{roots}}&=\sum_j k_j
=\frac{2\pi}{L}\left(\sum_j I_j+\sum_a J_a\right)\pmod{2\pi}.
\end{aligned}
```

These are the logarithmic Lieb-Wu equations; see the
[reference and equation provenance](../CITATIONS.md). The spin self-scattering
term is zero. `Lambda` in the API and root tables is the conventional Lieb-Wu
rapidity, not a rescaled quantity.

At half filling and Sz=0, the two even-length classes have different quantum-number parity:

| Ring length | Charge labels I (L values) | Spin labels J (M values) | Momentum |
| --- | --- | --- | --- |
| L=4m+2 | `-(L-1)/2, ..., (L-1)/2`, half-odd integers | `-(M-1)/2, ..., (M-1)/2`, integers | 0 |
| L=4m | `-L/2+1, ..., L/2`, integers | `-(M-1)/2, ..., (M-1)/2`, half-odd integers | pi |

Each list advances in steps of one. Charge labels span a full Brillouin zone;
for L=4m we choose the endpoint +pi rather than -pi. The API stores labels
exactly using Uni20's `half_int`, and momentum as both an exact integer index
and `2*pi*index/L`. The parity conventions and momenta are tested against
small-ring exact diagonalization, including the fermionic sign of translation.
For polarized half filling, the same construction is controlled by the parity
of M, not L/2. In supported interacting doped root sectors, N is even and M
odd: charge labels are centered half-odd integers and spin labels are centered
integers. Their root momentum is zero. Exact symmetry energy and momentum
offsets are applied when the physical sector differs from the root sector.

## How the numerical solve works

We exploit reflection symmetry, solving only for positive k and Lambda.
The partners are reflected exactly. For even M at half filling, k=0 and k=pi
are fixed; for odd M, the central spin rapidity is exactly zero. In particular, the
fixed k=pi contributes $\sin (k)=0$ exactly, avoiding division of a floating
approximation to sin(pi) by a very small U.

The reduced system uses damped Newton updates with an analytic Jacobian.
Uni20 performs the dense linear solve in the selected scalar type; there is
no conversion through double. Spin variables are internally scaled as
$\Lambda /\max (1,u)$ to condition the large-U equations, and converted back on
return. For U<8 we first solve at U=8 and halve U stage by stage down to the
requested value, preserving physical Lambda between stages. For U>=8 we
solve directly at the requested interaction.

The reduced dense solve has approximately (N+M)/2 unknowns, at most 3L/4:
O((N+M)^2) storage and O((N+M)^3) work per Newton update. This
is a starting point for moderate rings, not yet a matrix-free large-L solver.

Diagnostics include separate $\max \lvert F_{\mathrm{charge}} \rvert /L$ and $\max \lvert F_{\mathrm{spin}} \rvert /L$, their
maximum, accepted Newton updates, and the number of completed continuation
stages. The default tolerance is 32 times the selected type's epsilon.
The iteration budget applies to **all stages together**, not to each stage.
A zero budget evaluates the initial large-U seed without taking a Newton
step; the exact U=0 calculation requires no updates.

If the budget is exhausted or line search cannot improve the residual, the
returned state is explicitly marked unconverged (CLI exit status 2).
Even if continuation stops early, the returned roots and residuals use the
**final root-sector interaction**, not a solved intermediate interaction.
For attractive U this is |U|. The energy includes the physical symmetry offset,
and the momentum field identifies the intended physical ground-state branch;
an unconverged energy must not be treated as a certified eigenvalue.

Weak coupling brings nearby roots close together, so a solve can stall at a
precision-dependent limit. Try higher precision, or a larger update budget
if the budget was exhausted. At extremely strong coupling the energy sums
terms of order one to obtain a result of order 1/U; cancellation limits
relative energy accuracy even if the residual is small. A small equation
residual is **not an energy-error bound**.

## Free fermions and strong coupling

At U=0, the two spin populations independently fill their lowest one-particle
levels. At half filling and Sz=0 this gives

```math
E=\begin{cases}
-4\cot(\pi/L),&L=4m,\\
-4\csc(\pi/L),&L=4m+2.
\end{cases}
```

We evaluate finite trigonometric sums directly, setting empty and filled-band
kinetic energies exactly to zero. The interacting Bethe coordinates
are singular at U=0, so `charge_momenta` instead holds the occupied
free-fermion momenta, with repetitions for opposite spins.
`spin_rapidities` and both Bethe-label lists are empty, and `free_fermion`
is true. For L=4m the free ground state is degenerate; we choose a
representative with momentum pi. This is an occupation-state representative,
not a claim to reproduce the interacting singlet's limiting wavefunction.
The same no-scattering path is used when a root sector contains only one spin
species, even at nonzero U; `free_fermion` identifies this path, not only U=0.

At large U and L>=4, virtual hopping generates an antiferromagnetic XXX chain
with exchange J=4/U. In the conventions of this repository,

```math
\frac U4 E_{\mathrm{Hubbard}}\longrightarrow E_{\mathrm{XXX}}(J=1)-\frac L4.
```

The offset matters: dropping it would compare different Hamiltonians.
The regression tests check this limit against the existing XXX solver, the
free-fermion limit against independently filled one-particle levels, and
finite-U energies and momentum against exact fermionic diagonalization
through six sites, across all supported spin populations and both signs of U.
Degenerate ground eigenspaces are Fourier-projected to check momentum.
Native-precision two-site tests distinguish long-double
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

// Doped attraction; physical N=16, Sz=0, auxiliary repulsive half-filled roots.
auto paired = bethe::hubbard::sector_ground_state<long double>(
    32, 16, uni20::half_int{0}, -4.0L);
// Always inspect paired.root_particles, root_down_spins and root_interaction
// before interpreting its roots. paired.energy and momentum are physical.
```

`State<Real>` carries the two root families, separate quantum-number lists,
energy, momentum, residual diagnostics, requested interaction, and particle
counts, together with explicit root-sector metadata and symmetry offsets.
`SolverOptions<Real>` supplies `residual_tolerance` and
`max_iterations`, as for the other models. Link the public `bethe::bethe`
CMake target, which includes the Uni20 linear-algebra dependency.

The nested state type intentionally differs from the spin-chain `RealState`:
future charge/spin excitations should extend this
module without overloading the meaning of a single rapidity vector. There
is no arbitrary-quantum-number solve or Hubbard excitation scan yet; a
complete spectrum would also require complex spin and charge-spin strings.
