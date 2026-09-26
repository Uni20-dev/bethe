# Takhtajan–Babujian: thermodynamic spinons and continua

[Overview](../README.md) · [Finite-ring ground states](takhtajan-babujian.md) · [Output](output.md)

`bethe-tb-dispersion` supplies zero-field excitation references for the critical
spin-1 Takhtajan–Babujian (TB) chain: an elementary spinon line and two-/four-spinon
energy bounds. These are **infinite-chain kinematic curves**, not finite-ring
levels, form factors or a complete spectrum.

```sh
bethe-tb-dispersion --points 129 --csv tb.csv
bethe-tb-dispersion --branch spinon --momentum 0.3 --precision long-double
bethe-tb-dispersion --branch four-spinon --folded --json quadrupole-bounds.json
bethe-tb-dispersion --points 65 --precision fp128 --tsv tb.tsv
```

The default is fp64, with native long-double and optional MPLAPACK-enabled fp128.
Each branch's grid spans its own momentum interval. `--momentum` replaces the
grid and must be valid for every selected branch. `--references` prints the
literature; normal reports retain the Hamiltonian, conventions, provenance and
CPU time. Spin labels use Uni20's `half_int`, exported as decimal 0.5 rather
than a plotting-unfriendly fraction. CSV/TSV/JSON, streaming and retention use
the same output facilities as the [ULS tool](su3-dispersion.md).

## A factor of four worth keeping explicit

Our Hamiltonian is

```math
H=J\sum_j\left[\mathbf S_j\cdot\mathbf S_{j+1}
 -(\mathbf S_j\cdot\mathbf S_{j+1})^2\right],\qquad J\gt0.
```

The spins are 1. There is no additive constant, and the bond energies in
total-spin sectors 0, 1 and 2 are -6J, -2J and 0. The existing finite-ring
`bethe-tb-pbc` tool uses J=1.

[Vlijm–Caux](../CITATIONS.md#vlijm-caux-2014), Eq. (1.2), has a prefactor
$`J_{\mathrm{paper}}/4`$: **our J=1 means their J=4**, not J=1.
For the common Hamiltonian with coefficients cos(theta), sin(theta), use
$`J=1/\sqrt{2}`$ at $`\theta=-\pi/4`$. This is not the generic bilinear
spin-1 Heisenberg chain or the positive-biquadratic [ULS point](su3.md).

## Spinons and physical sectors

With lattice spacing and hbar equal to one,

```math
\epsilon(p)=v\sin p,\qquad v=2\pi J,\qquad 0\le p\le\pi.
```

The spinon is gapless at both endpoints and carries **spin 1/2**, although each
physical site has spin 1. A periodic chain of integer-spin sites has only
integer total spin: an isolated spinon is not one of its eigenstates.

The two-string sea's density and dressed hole energy give an independent way
to identify the scale:

```math
\rho(\lambda)=\frac{1}{2\cosh(\pi\lambda)},\qquad
p(\lambda)=2\arctan(e^{-\pi\lambda}),\qquad
\epsilon(\lambda)=4\pi J\rho(\lambda).
```

Two spinons combine into S=0 or S=1, while four allow S=0,1,2. From the singlet
ground state, a local spin operator selects S=1, a scalar bond/dimer operator
selects S=0, and a local rank-two quadrupole selects S=2. The last therefore
requires **at least four spinons**, even though its kinematic lower threshold
can coincide with the two-spinon one. These selection rules permit matrix
elements; they do not prove nonzero weight at every momentum or on an edge.
For example, the Q=0 total-spin operator annihilates the singlet.

Identical sine dispersions do not make this the XXX spectrum with a rescaled J.
TB's low-energy theory is SU(2) at level 2 (central charge 3/2), and its spinons
have additional non-Abelian state-counting structure. We do not infer finite-size
degeneracies, fusion multiplicities, or spectral weights from the curves.
See [Frahm–Stahlsmeier](../CITATIONS.md#frahm-stahlsmeier-1998).

## Two- and four-spinon bounds

For total one-site momentum $`0\le Q\le2\pi`$, the two-spinon envelope is

```math
E_{2,-}(Q)=v|\sin Q|,\qquad E_{2,+}(Q)=2v\sin(Q/2).
```

It follows by extremizing $`\epsilon(p)+\epsilon(Q-p)`$ over allowed individual
momenta. The minimum puts one spinon at a soft endpoint; the maximum shares
the total momentum equally. The shared sine/elliptic-band helper also used by
XXX/XXZ supplies this calculation; the TB layer owns normalization and sectors.

For four spinons, let $`q=\min(Q,2\pi-Q)`$. Including total-momentum aliases gives

```math
E_{4,-}(Q)=v|\sin Q|,\qquad E_{4,+}(Q)=4v\cos(q/4).
```

Soft spectator spinons preserve the minimum. Concavity of sin(p) on [0,pi]
puts the maximum at four equal momenta, choosing the alias nearest a total
of 2pi. At Q=0, the two-spinon interval collapses to zero but the four-spinon
upper edge is 4v: ignoring the modulo-2pi alias would miss it entirely.
At Q=pi their upper edges are 2v and $`2\sqrt{2}v`$ respectively.
Neither upper edge bounds the full spectrum with arbitrarily many spinons.

`--folded` takes the envelope of Q and Q+pi, corresponding to a two-site iMPS
unit cell. `p` remains in radians/site, while `cell_momentum` is
$`K=2p\pmod{2\pi}`$ in radians/cell. The lower edge is unchanged. With
$`r=\min(q,\pi-q)`$, the folded upper edges are $`2v\cos(r/2)`$ and
$`4v\cos(r/4)`$. Folding is bookkeeping for the chosen cell, not evidence of
a gapped or spontaneously dimerized ground state at the critical TB point.

## Finite chains are checks, not the thermodynamic answer

[Independent ED](../scripts/reference_tb_excitations.py) constructs the spin-1
bond matrix, projects lattice translations, and resolves total S using the
Casimir. It does not use the Bethe equations or dispersion formulas. Fixed Sz
alone is insufficient to label total S. For J=1, selected gaps are:

| L | S=1, Q=2pi/L | S=1, Q=pi | S=2, Q=0 |
|---|---:|---:|---:|
| 4 | 8.167056259933 | 3.403124237433 | 6.574697112687 |
| 6 | 6.207722752576 | 2.225739367652 | 4.686482195819 |
| 8 | 4.852117819006 | 1.661960361397 | 3.654190748086 |
| 10 | 3.948711041553 | 1.328463417081 | 2.998020406819 |

The small-momentum triplet velocity approaches 2pi. The staggered triplet has
scaling dimension 3/8, with significant logarithmic finite-size corrections;
its positive finite-L gap should not be mistaken for a thermodynamic gap.
Tests compare the existing finite-root ground solver against the independent
ED normalization while retaining its finite string deviations. Replacing those
roots by ideal two-strings is not an exact finite-size excitation method.

Native-precision regressions also check an 80-digit dispersion reference,
independent rapidity-space expressions, exchange scaling, direct momentum scans,
four-particle convolutions and folding. Unrepresentable positive energies or
bounds produce `precision_limit`, missing values (JSON null), and exit 2.
Invalid parameters exit 1; an exchange that overflows the whole elementary band
is rejected before opening exports. The implementation evaluates the soft
endpoints exactly and does not silently replace underflowing positive energies
by gapless zeros. There are no numerical quadrature or root-solving tolerances
for these analytic curves.

## Library interface

```cpp
#include <bethe/tb_dispersion.hpp>
bethe::takhtajan_babujian::SpinonDispersion<long double> band(1.0L);
auto spin = band.spin(); // Uni20 half_int, exactly 1/2
auto energy = band.energy(0.3L);
auto pair = band.continuum(3.0L);
auto four = band.four_spinon_continuum(3.0L, bethe::SpinonMomentum::folded);
```

Finite-ring excited-string enumeration, fields, open ends, spectral intensities
and general spin-s chains are not implemented by this frontend.
