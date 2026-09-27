# SU(3)/ULS: elementary excitations and continuum bounds

[Overview](../README.md) · [Finite rings](su3.md) · [Output](output.md)

`bethe-su3-dispersion` gives zero-field, infinite-chain reference curves for
the fundamental SU(3) permutation chain. Its principal application is comparing
momentum-resolved MPS excitations at the spin-1 Uimin–Lai–Sutherland point.
It reports elementary lines and **specified-particle-content continuum bounds**,
not spectral weights or finite-ring excited energies.

```sh
bethe-su3-dispersion --points 129 --csv uls.csv
bethe-su3-dispersion --branch bar3 --precision long-double
bethe-su3-dispersion --branch four-soliton --folded --json edges.json
bethe-su3-dispersion --branch 3 --momentum 0.3 --precision fp128
```

The default is fp64; long-double and optional MPLAPACK-enabled fp128 keep native
arithmetic, input parsing and output precision. `--branch all` samples each
branch over its own momentum interval; a single `--momentum` must lie in every
selected interval. `--references` prints the literature. CSV/TSV/JSON exports
include Hamiltonian, precision, momentum conventions, provenance and CPU time.

The only output table is `dispersion`, which is also the primary table.
`--csv uls.csv` and `--csv-table dispersion=uls.csv` are equivalent. The
`3`, `bar3`, `two-soliton` and `four-soliton` branches are rows within this table,
selected with `--branch`, not separate table names.

Unrepresentable continuum bounds or underflowing positive energies are omitted,
with `precision_limit` and exit 2;
invalid parameters exit 1. No optimization tolerance is needed: the bounds below
are analytic, including all stationary points and endpoint cases.

## Hamiltonian and units

```math
H_P=J\sum_j P_{j,j+1},\qquad J\gt0,
```

where each site is a fundamental SU(3) representation. In a spin-1 basis,

```math
P_{j,j+1}=\mathbf S_j\cdot\mathbf S_{j+1}
 +(\mathbf S_j\cdot\mathbf S_{j+1})^2-1.
```

Consequently the ULS Hamiltonian with equal bilinear and biquadratic coefficients
J has total energy $`E_{\rm ULS}=E_P+JL`$ on a ring; **excitation gaps are unchanged**.
For coefficients $`\cos\theta,\sin\theta`$ at $`\theta=\pi/4`$, use
$`J=1/\sqrt{2}`$. This is a nearest-neighbour model, not the SU(3) Haldane–Shastry chain.

## Elementary branches and representations

Let $`a=\pi/3`$ and $`A=4\pi J/(3\sqrt{3})`$. The two branches are

```math
\begin{aligned}
\epsilon_{3}(p)&=A[1/2-\cos(p+a)], &0\le p\le4a,\\{}
\epsilon_{\bar3}(p)&=A[\cos(p-a)-1/2], &0\le p\le2a.
\end{aligned}
```

Both vanish at their endpoints, with endpoint speed $`v=2\pi J/3`$.
Their maximum energies are $`3A/2`$ and $`A/2`$, respectively. The implementation
uses a product of sines rather than a cosine difference to preserve tiny
positive energies near a soft endpoint.

The first nested sea's hole carries highest weight $`(1,0)`$, representation 3;
the second carries $`(0,1)`$, representation bar3. The branch ranges follow
Eq. (64) of [Vörös–Penc](../CITATIONS.md#voros-penc-2021), with twice their energies.
The prose following their Eq. (65) interchanges the ranges; we use Eq. (64),
independently checked against the nested-sea densities and representation
assignment of [Doikou–Nepomechie](../CITATIONS.md#doikou-nepomechie-1998).
Our positive momentum convention reverses the negative hole momenta of the
latter reference; reflection leaves these bands and continuum envelopes invariant.

An isolated elementary line is not a local-response quasiparticle pole. The
SU(3) centre charge (triality) of a periodic L-site fundamental chain is L mod 3.
A configuration of $`n_3,n_{\bar3}`$ holes must therefore obey
$`n_3+2n_{\bar3}=L\pmod3`$. In particular, neither single branch is by itself a
state of the balanced L=0 mod 3 ring.

The smallest neutral pair has $`3\otimes\bar3=1\oplus8`$. Longer neutral
combinations include three equal representations and two such pairs. The
spin-1 spin operators and five quadrupoles together form the SU(3) adjoint 8;
from a singlet vacuum they access adjoint states, not the singlet part of the
pair. The adjoint decomposes into spin 1 and spin 2 under physical spin rotation.
Symmetry permits these channels but does not determine their matrix elements.
At exactly Q=0 the summed on-site generators annihilate the singlet, even though
the kinematic four-particle continuum has nonzero width.

## Two versus four solitons

`two-soliton` means one 3 and one bar3. `four-soliton` means two of each;
its representation product is $`(1\oplus8)\otimes(1\oplus8)`$ and includes
additional adjoint channels. Auxiliary Bethe strings select multiplets; they do
not add an independent energy to these thermodynamic hole sums.

All continuum momenta are total one-site momenta modulo $`2\pi`$. Set
$`q=\min(Q,2\pi-Q)`$, so $`0\le q\le\pi`$. The two-particle bounds are

```math
\begin{aligned}
E_{2,-}(q)&=\begin{cases}
\epsilon_{\bar3}(q),&0\le q\le2a,\\{}
\epsilon_3(q-2a),&2a\le q\le3a,
\end{cases}\\{}
E_{2,+}(q)&=\begin{cases}
\epsilon_3(q),&0\le q\le a,\\{}
2A\sin(q/2),&a\le q\le3a.
\end{cases}
\end{aligned}
```

To derive them, minimize/maximize the sum at fixed Q over
$`p_{\bar3}\in[\max(0,Q-4a),\min(2a,Q)]`$. The sum is
$`2A\sin(Q/2)\cos[p_{\bar3}-(Q-a)/2]`$. Testing its endpoints and stationary
point gives the expressions above, without a discretized search.

For two pairs, convolution of those intervals, including total-momentum aliases,
gives

```math
\begin{aligned}
E_{4,-}(q)&=\begin{cases}
\epsilon_{\bar3}(q),&0\le q\le2a,\\{}
\epsilon_{\bar3}(q-2a),&2a\le q\le3a,
\end{cases}\\{}
E_{4,+}(q)&=4A\cos(q/4).
\end{aligned}
```

For the lower bound, periodically extend the narrow branch with period 2a. It
is subadditive under momentum addition; each elementary energy bounds it from
above. A soft pair at momentum 0, 2a or 4a saturates the required shifted
two-particle minimum. For the upper bound, each pair obeys
$`E_{2,+}(Q)\le2A\sin(Q/2)`$; sharing the appropriate total-momentum alias equally
between the pairs saturates the maximum.

In particular, at Q=pi the two-particle onset is A but the four-particle onset
is A/2. Plotting only the former misses lower-energy adjoint excitations between
the soft momenta. Conversely, **neither upper edge bounds the full spectrum**:
additional particles extend it. These outputs do not calculate spectral intensity,
threshold exponents, or claim a finite-size state count.

## Three-site iMPS momentum

The table always retains one-site p (radians/site). `cell_momentum` is
$`K=3p\pmod{2\pi}`$ (radians/three-site cell). Three distinct one-site images
map to the same cell momentum. `--folded` takes their envelope for each continuum;
it does not relabel a single elementary band as a complete cell spectrum.

If r is the distance to the nearest multiple of 2a, $`0\le r\le a`$, the folded
lower bound is $`\epsilon_{\bar3}(r)`$ for both particle contents. Their upper
bounds are $`2A\cos[(a-r)/2]`$ and $`4A\cos(r/4)`$ respectively. A folded envelope
does not by itself specify internal gaps, translation eigenvectors or weights.

## Validation and finite-size caveats

Typed tests compare 80-digit momentum-space references, independently expressed
nested-sea densities, endpoint velocities, exchange scaling, direct constrained
momentum scans, pair-continuum convolutions, and three-image folding. Tests cover
fp64, long-double and fp128, including energies at momenta below machine epsilon.

[The independent ED script](../scripts/reference_su3_excitations.py) swaps colour
words, projects one-site momentum, and uses the total SU(3) Casimir to isolate
the adjoint. Fixed colour populations alone do **not** isolate a representation.
For J=1 the lowest adjoint gaps include:

| L | Q=2pi/L | Q=2pi/3 | Q=pi |
|---|---:|---:|---:|
| 6 | 2.036701665932 | 1.401149725135 | 3.605551275464 |
| 9 | 1.430375396369 | 0.932235847198 | not an allowed momentum |
| 12 | 1.090728755289 | 0.700304961905 | 2.018878869993 |

The small-momentum velocity approaches 2pi/3. The soft tower has conformal
dimension 2/3, with logarithmic finite-size corrections: its finite-ring gap is
positive even though the thermodynamic threshold vanishes. At L=12, Q=pi, the
adjoint level is already below the thermodynamic two-soliton lower edge
2.418399152312, consistent with the lower multiparticle branch. It is not equal
to the thermodynamic four-soliton onset 1.209199576156; small rings are checks of
normalization, sectors and scaling, not exact continuum-bound fixtures.

## C++ API

```cpp
#include <bethe/su3_dispersion.hpp>
bethe::su3::ExcitationDispersion<long double> band(1.0L);
auto energy = band.energy(bethe::su3::Particle::antifundamental, 0.3L);
auto edges = band.continuum(3.0L, bethe::su3::Continuum::four_soliton,
                           bethe::su3::Momentum::three_site_folded);
```

The existing finite-ring ground-state API is unchanged. Other populations,
finite-ring excited roots, fields, open ends and general SU(N) excitation tools
remain separate extensions.
