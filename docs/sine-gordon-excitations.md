# Sine-Gordon particle lines and asymptotic levels

[Overview](../README.md) · [Vacuum solver and conventions](sine-gordon.md) · [Exports](output.md)

Two tools cover the first excitation checkpoint:

- `bethe-sine-gordon-dispersion`: infinite-volume particle dispersions and
  two-particle lower thresholds, in both attractive and repulsive regimes.
- `bethe-sine-gordon-bethe-yang`: selected large-volume two-soliton levels.
  These are **asymptotic**, not exact finite-volume excited energies.

Both use fp64 by default, native long-double and optional MPLAPACK fp128.
All dimensional quantities use velocity=hbar=1. The mass input M is the physical
soliton mass, not the bare cosine coefficient or the lightest breather mass.

```sh
bethe-sine-gordon-dispersion --p 0.4 --points 129 --csv sg-lines.csv
bethe-sine-gordon-dispersion --p 0.5 --momentum 0.7 --precision long-double
bethe-sine-gordon-bethe-yang --p 2 --length 10 --levels 4 --json sg-pairs.json
bethe-sine-gordon-bethe-yang --p 1 --length 10 --number 1.5 --precision fp128
```

Use `--references` for literature; ordinary reports include conventions,
provenance, numerical status and CPU time. CSV/TSV/JSON share the usual streaming,
retention and overwrite-protection facilities. Half-integer Bethe labels are
exported as decimal numbers.

## Particles, charges and thresholds

We retain the vacuum solver's canonical kinetic term and coupling:

```math
p=\frac{\beta^2}{8\pi-\beta^2},\qquad
m_s=m_{\bar s}=M,\qquad m_{B_n}=2M\sin\frac{n\pi p}{2}.
```

Soliton and antisoliton have topological charges +1 and -1. Breathers are
neutral and exist for positive integers satisfying **np<1**. A pole reaching
np=1 is at the soliton–antisoliton threshold and is excluded from our stable
particle list. Thus p=1 has no breathers, while p=1/2 has only B1.
Existence is decided in the selected scalar arithmetic; use higher precision
to distinguish couplings closer to a threshold than that arithmetic resolves.
See [Fehér–Takács](../CITATIONS.md#feher-takacs-2011), Eq. (2.2).

For physical momentum k, the curves are

```math
E_a(k)=\sqrt{m_a^2+k^2},\qquad
E_{a_1\cdots a_r,\min}(k)=\sqrt{(m_{a_1}+\cdots+m_{a_r})^2+k^2}.
```

The second relation follows by minimizing at a common velocity, not at equal
momenta for unequal masses. It is an infimum for specified particle content,
with **no finite upper edge**. A lower-lying bound particle is a separate branch,
not part of the scattering continuum's onset. Charges add; energies alone do
not establish whether a chosen local operator couples to that sector.

The frontend enumerates all stable species and unordered pairs, including
identical pairs. `--branch particles|thresholds|all` selects rows. A default grid
covers k=0..5M; `--max-momentum` changes its endpoint and `--momentum` selects a
single signed value instead. There is no Brillouin zone or unit-cell folding in
this continuum theory. Very weak coupling can produce many species: enumeration
is limited to 256 breathers and one million output rows, with an explicit error
rather than silent truncation. The library can evaluate a chosen B_n directly.

## Scattering and Bethe–Yang convention

The library exposes the odd same-charge phase for real rapidity at all p>0:

```math
S_{ss}(\theta)=S_{\bar s\bar s}(\theta)=-e^{i\chi_p(\theta)},\qquad
\chi_p(0)=0,\qquad \chi'_p(\theta)=2\pi G_p(\theta).
```

The constant minus sign is kept separate. The phase integrates the existing
vacuum kernel's Fourier multiplier with sin(k theta)/k, using native compensated
quadrature and a separate analytic Fourier-tail bound. Its convention follows
[Fehér–Pálmai–Takács](../CITATIONS.md#feher-palmai-takacs-2012), Eq. (2.1).

The current Bethe–Yang solver selects **same-charge pairs at p>=1**, with total
momentum zero, rapidities ±theta, and topological winding charge +2 or -2:

```math
ML\sinh\theta+\chi_p(2\theta)=2\pi I,\qquad
I=\tfrac12,\tfrac32,\ldots,\qquad E_{\mathrm{BY}}=2M\cosh\theta.
```

These winding sectors allow the field to change by 2pi Q/beta around the circle;
they are not neutral states with strictly zero winding. `--charge` selects ±2;
`--number` is the first positive I and `--levels` adds consecutive values.
The free solution brackets the positive root in the repulsive regime. The
constant S-matrix sign accounts for the half-odd quantization. At p=1, chi=0,
giving momenta ±2pi I/L directly.

We do not solve mixed soliton–antisoliton transfer-matrix channels, breather
Bethe–Yang levels, or nonzero-total-momentum pairs here. In particular, simply
relabeling this pair as neutral would give the wrong state-selection rules.

Bethe–Yang omits wrapping/vacuum-polarization corrections. No minimum ML is
enforced: small-volume output is still a solution of the stated asymptotic
equation, **not** a certified approximation to the exact level there. Adding
the separate vacuum energy does not restore the omitted excited-state dressing.
For exact selected finite-volume levels, use the separate
[excited-state NLIE solver](sine-gordon-excited.md): p>=1, opposite rapidities,
I=0.5 or 1.5, charge +/-2. It resolves the sea and subtracts a separately
converged vacuum to obtain a gap.

## Numerical diagnostics and checks

The `levels` table separates the counting residual, Fourier quadrature estimate,
tail bound and phase status. `--tolerance` is absolute in the counting equation,
default 16384 native epsilons; each phase gets tolerance/16. Iteration and
per-phase quadrature/cutoff budgets are independent. None bounds the omitted
physical wrapping corrections. Energy and rapidity are absent on failed solves
(JSON null), with exit 2; invalid parameters exit 1 before exports are opened.
Dispersion rows likewise omit unrepresentable energies. Unrepresentable input
mass scales are rejected before exporting a spectrum.

Checks include the free Dirac spectrum, the p=1/2 closed phase
$`\chi=-\arctan(\sinh\theta)`$, the p=2 rapidity-space kernel, near-free
native-precision response, and the original exponential quantization condition.
[Independent 90-digit references](../scripts/reference_sine_gordon_excitations.py)
use a closed rapidity kernel rather than the production Fourier integral.
Library regressions cover unequal-mass threshold minimization, scaling,
breather appearance boundaries and distinct numerical-budget failures.

## The two-flavour Schwinger comparison

For equal small fermion masses at theta=0, the **leading light-sector effective
theory** has canonical beta squared 2pi, hence p=1/3. Soliton, antisoliton and B1
form a mass-degenerate triplet; B2 is a singlet with mass sqrt(3)M. B3 is at
threshold and is not counted. These are useful light-particle MPS targets.

This does not solve the full massive Schwinger model, its heavy eta mode, or a
finite-spacing gauge Hamiltonian. Corrections to the leading effective theory
and lattice artifacts are separate. The [DMRG study by Itou, Matsumoto and
Tanizaki](../CITATIONS.md#itou-matsumoto-tanizaki-2024), particularly Sec. 2,
explains the bosonization and the triplet/singlet comparison.

## Library entry points

```cpp
#include <bethe/sine_gordon_particles.hpp>
#include <bethe/sine_gordon_bethe_yang.hpp>
namespace sg = bethe::sine_gordon;
sg::ParticleSpectrum<long double> spectrum(1.0L, 0.4L);
auto energy = spectrum.energy({sg::ParticleKind::breather, 1}, 0.7L);
auto phase = sg::soliton_phase(0.7L, 2.0L);
auto pair = sg::same_charge_pair(1.0L, 10.0L, 2.0L, uni20::from_twice(1));
// Check pair.converged before using *pair.energy or *pair.rapidity.
```
