# XXX chains: states, spinons and spectral weights

[All models](models.md) · [First calculation](getting-started.md)

The spin-$`1/2`$ Heisenberg chain is a useful starting point for both finite-size
energies and MPS excitation comparisons. The guides use $`J=1`$ and spin
operators $`\mathbf S=\boldsymbol\sigma/2`$. Match those conventions before
comparing energies or intensities.

| What to calculate | Program | Guide |
| --- | --- | --- |
| Periodic ground state, sector minima and explicit Bethe states | `bethe-xxx-pbc` | [Periodic XXX](xxx.md) |
| Free-end ground state, sector minima and real-root excitations | `bethe-xxx-obc` | [Open XXX](open-chains.md) |
| Selected periodic real-root excitations | `bethe-xxx-pbc --excitations` | [Excitation families and enumeration](excitations.md) |
| Odd-ring one-spinon branch or thermodynamic spinon kinematics | `bethe-xxx-pbc --spinons`; library formulas | [Single spinons](spinons.md) |
| Finite-ring two-spinon Lehmann spectrum | `bethe-xxx-structure-factor` | [Weights, momentum and sum rules](xxx-structure-factor.md) |
| Infinite-chain selected spinon spectral densities | `bethe-xxx-structure-factor-thermo` | [Two-spinons](xxx-structure-factor-thermo.md), [selected sectors](xxx-spinon-sectors.md) |
| Direct fixed-momentum and global spinon moments | `bethe-xxx-moments-thermo` | [Moment definitions, controls and exports](xxx-spinon-moments.md) |

## Choose a benchmark

Start with [finite chains](tutorials/xxx.md) to understand quantum numbers,
open versus periodic boundaries and the selected excitation family. For
precision-sensitive gaps, use [the precision tutorial](tutorials/xxx-precision.md).

For a dynamical comparison, distinguish the [finite-ring spectral lines](tutorials/xxx-structure-factor.md)
from the [unbroadened thermodynamic density](tutorials/xxx-structure-factor-thermo.md).
The finite-ring tool computes the supported two-spinon family. The thermodynamic
tool defaults to two spinons and can include four and six; these remain partial
responses. Broadening discrete lines is a plotting choice; it does not supply
the missing state families.

## Interpret the labels

Total spin, magnetization and a particular root family are different labels.
The excitation guide explains highest-weight states and what `all` enumerates.
An open chain has standing waves rather than a lattice momentum. For a periodic
chain, the physical momentum and any unit-cell folding must be matched to the
comparison calculation. Every model result still needs its convergence status
and numerical work limits; choosing MPFR does not enlarge the physical family.

For thermodynamic four- and six-spinon spectra, see
[selected sectors, explicit proposals and diagnostics](xxx-spinon-sectors.md).
