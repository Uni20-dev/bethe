# Hubbard chains: finite ground states and spin-charge excitations

[All models](models.md) · [Half-filled tutorial](tutorials/hubbard-half-filled.md) ·
[Doped tutorial](tutorials/hubbard-doped.md)

Hubbard calculations have separate charge momenta and spin rapidities. A finite
ground state, a thermodynamic elementary dispersion and a two-particle continuum
edge are distinct outputs, with different restrictions.

| What to calculate | Program | Coverage and guide |
| --- | --- | --- |
| Finite periodic ground states | `bethe-hubbard-pbc` | [Periodic Hubbard](hubbard.md): repulsive half filling, balanced attraction and selected doped shells on even rings; unrestricted free limit |
| Finite free-end ground states | `bethe-hubbard-obc` | [Open Hubbard](hubbard-open.md): all physical fillings/spin projections, either interaction sign and either length parity |
| Half-filled spinon, holon and antiholon lines | `bethe-hubbard-dispersion` | [Half-filled dispersions](hubbard-dispersion.md), U>0 and zero field |
| Doped spinon and charge lines | `bethe-hubbard-dispersion` | [Doped dispersions](hubbard-doped.md), U>0 and zero field |
| Half-filled two-particle continuum edges | `bethe-hubbard-continuum` | [Continuum guide](hubbard-continuum.md): spinon pairs and charge-containing families |

## Match the Hamiltonian and energy reference

Finite-system tools use hopping $`t=1`$ and the unshifted interaction
$`U n_{\mathrm{up}}n_{\mathrm{down}}`$. The dispersion tool defaults to the
symmetric interaction convention, useful for SO(4) formulations. Particle-changing
excitations also require a choice of Hamiltonian or Fermi-level energy reference.
The tutorials show these shifts explicitly; matching the ground-state energy
alone does not fix a particle-addition/removal dispersion convention.

The [sector note](hubbard-sectors.md) explains finite periodic branch restrictions.
Open-chain symmetry mappings can solve an auxiliary sector; interpret its roots
together with the mapping flags and energy offset in the physical result.

## From elementary lines to continuum edges

The continuum solver combines constituent dispersions and optimizes at fixed
total momentum. It gives family boundaries, not spectral weights. The
[charge-continuum method note](hubbard-charge-continuum-design.md) documents
conventions, extrema searches and numerical budgets behind those calculations.
Quadrature, momentum inversion and extrema convergence are separate from the
number of arithmetic bits.
