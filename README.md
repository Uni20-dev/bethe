# Bethe ansatz solvers

C++23 tools for finite-system and thermodynamic Bethe ansatz calculations, complementing
[Uni20](https://github.com/Uni20-dev/uni20) and the
[Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

High precision is a first-class feature: the solvers support
**fp64**, **long-double** (platform-dependent extended precision), and optional
**fp128** (binary128, about 34 significant decimal digits) through Uni20/MPLAPACK.
Parameters, numerical calculations, excitation gaps, and output retain the
selected precision. See [precision controls](docs/command-line.md#select-arithmetic-precision)
and [binary128 setup](docs/building.md#enable-binary128).

A calculation starts with a model, a boundary condition, and a choice of
Bethe quantum numbers. The solver finds the corresponding rapidities and
uses them to calculate energies and, for periodic chains, momenta. You can
begin with a ground state without choosing the quantum numbers yourself,
then explore magnetization sectors, excitation families, and spinons.

## Choose a model

Commands follow `bethe-<model>-<boundary>`: `pbc` means periodic boundaries,
and `obc` means open boundaries (currently free ends, with no boundary fields).
Nonspatial models such as Richardson pairing omit the boundary suffix.
Thermodynamic dispersion tools use `bethe-<model>-dispersion`.

- [Periodic XXX](docs/xxx.md): `bethe-xxx-pbc` gives ground states,
  sector minima, real-root excitations, and the odd-chain one-spinon branch.
- [Free-end XXX](docs/open-chains.md): `bethe-xxx-obc` supports
  ground states, sector minima, and real-root excitations, without lattice momentum.
- [Periodic XXZ](docs/xxz.md): `bethe-xxz-pbc` supports ground states and sector minima
  for `Delta >= 0` (also `-1 < Delta < 0` on even rings), plus restricted
  real-root excitations for `0 <= Delta <= 1`.
- [Free-end XXZ](docs/xxz-open.md): `bethe-xxz-obc` supports ground states and
  sector minima for `Delta > -1`, including massive boundary roots, and restricted
  real-root excitations for `0 <= Delta <= 1`, without lattice momentum.
- [Periodic Hubbard](docs/hubbard.md): `bethe-hubbard-pbc` gives repulsive
  half-filled spin sectors, balanced attractive ground states at any even filling,
  and selected doped sectors on even rings, plus the unrestricted U=0 limit.
- [Free-end Hubbard](docs/hubbard-open.md): `bethe-hubbard-obc` gives ground
  states at every physical filling and spin projection, for either sign of U
  and odd or even lengths, without lattice momentum.
- [Hubbard dispersions](docs/hubbard-dispersion.md): `bethe-hubbard-dispersion`
  gives half-filled spinon/holon/antiholon lines and [doped spinon/charge lines](docs/hubbard-doped.md)
  for U>0 at zero field, with symmetric (SO(4)) or unshifted interactions and
  Hamiltonian or Fermi-level energy references for iMPS comparisons. Supports
  [CSV/TSV/JSON exports](docs/output.md) alongside screen output.
- [Periodic Lieb–Liniger](docs/lieb-liniger.md): `bethe-lieb-liniger-pbc` gives
  repulsive continuum-boson ground states, specified Bethe states, and
  excitation scans within an explicit finite quantum-number window.
- [Periodic SU(3)](docs/su3.md): `bethe-su3-pbc` gives the balanced singlet
  ground state of the permutation chain for lengths divisible by three,
  also covering the spin-1 ULS point through an energy shift.
- [Periodic Gaudin–Yang](docs/gaudin-yang.md): `bethe-gaudin-yang-pbc` gives
  repulsive continuum-fermion ground states with odd populations of both spins,
  plus unrestricted free and fully polarized limits.
- [Supersymmetric t–J](docs/tj.md): `bethe-tj-pbc` gives periodic sector ground
  states at t=1, J=2 for odd spin populations, plus all no-hole and fully
  polarized sectors, with double occupancy excluded.
- [Spin-1 Takhtajan–Babujian](docs/takhtajan-babujian.md): `bethe-tb-pbc`
  gives the periodic even-ring singlet ground state of `H=sum[S.S-(S.S)^2]`,
  retaining complex-root finite-size string deviations.
- [Free-end spin-1 biquadratic](docs/biquadratic.md): `bethe-biquadratic-obc`
  gives even-chain ground energies, TL module minima and restricted real-root
  excitations of `H=-sum(S.S)^2`, with physical representation multiplicities
  kept separate from the auxiliary XXZ equations and spin labels. A separate
  [Q-system mode](docs/xxz-open-qsystem.md) adds complex-root levels and bounded
  small-chain spectrum searches. A [targeted two-string solver](docs/xxz-open-two-string.md)
  supplies the low-lying complex-root singlet on long chains (`--singlet-excitation`).
  The [ferromagnetic sign](docs/biquadratic-ferromagnetic.md) adds an exact
  one-defect band, [long-chain two-defect bound pairs](docs/biquadratic-bound-pairs.md),
  [three-defect droplets](docs/biquadratic-bound-triples.md),
  [four-defect droplets](docs/biquadratic-bound-quartets.md),
  [real-root scattering windows](docs/biquadratic-scattering.md),
  [mixed pair-plus-defect scans](docs/biquadratic-pair-defect.md), and sign-aware real/complex-root levels.
- [Richardson pairing](docs/richardson.md): `bethe-richardson` gives attractive
  reduced-BCS ground energies for distinct levels and a specified pair/blocked
  sector, using regularized variables through pair-root collisions.
- [Central spin](docs/central-spin.md): `bethe-central-spin` gives fixed-Sz
  ground energies for distinct nonzero spin-1/2 bath couplings of either
  sign, with a central field of either sign or exactly zero.
- [SU(n) fermion gas](docs/su-fermions.md): `bethe-sun-fermions-pbc` gives
  repulsive multicomponent continuum ground states when every occupied
  population is odd, plus unrestricted free and single-component limits.
- [Integrable spin ladder](docs/ladder.md): `bethe-ladder-pbc` gives periodic
  zero-field ground energies and singlet-count sector minima, with the
  required four-spin coupling and either sign of the rung exchange.
- [Haldane–Shastry ring](docs/haldane-shastry.md): `bethe-haldane-shastry-pbc`
  gives exact ground/spin-sector energies and motif spectra with Yangian
  multiplicities, for even and odd inverse-chord-square spin-1/2 rings.
- [Sutherland gas](docs/sutherland.md): `bethe-sutherland-pbc` gives exact
  periodic bosonic ground and excited energies for a specified collision
  exponent, with explicit integer labels or bounded label-window scans.

The XXX and XXZ models use spin-1/2 operators, J=1, and zero magnetic field.
Hubbard uses hopping t=1. Finite-system tools use the unshifted interaction
`U*n_up*n_down`; the dispersion tool defaults to the symmetric convention.
Spin-chain excitation scans cover explicitly supported real-root families, **not complete
spectra**: `--excitations all` means all states in that family. General complex-string
and infinite-root descendant scans are not implemented.

[Analytic thermodynamic spinon dispersions](docs/spinons.md#thermodynamic-dispersion)
are also available for XXX and gapless XXZ; these are distinct from finite-chain
calculations.
An [explicit XXZ spin-helix library API](docs/xxz-spin-helix.md) also provides
special eigenstates at commensurate couplings; these are not ground-state scans.

## Build

You need CMake 3.28+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
and Uni20's numerical dependencies. The default configuration fetches a
tested Uni20 revision:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

See [the build guide](docs/building.md) for a sibling Uni20 checkout, dependency
setup, optional binary128 support, or embedding `bethe::bethe` in another project.
For a checkout shared between hosts, keep builds on machine-local storage;
`build/` below is just the example build directory.

## First calculations

Start with the ground state of a four-site periodic XXX chain:

```sh
build/bethe-xxx-pbc 4
```

Its total energy is E=-2 in our spin-1/2 normalization. The report includes
energy per site, convergence diagnostics, and solver CPU time. To see how the
boundary condition or anisotropy changes the problem, try:

```sh
build/bethe-xxx-obc 4
build/bethe-xxz-pbc 4 --delta 0.5
build/bethe-xxz-obc 4 --delta 0.5
```

For mobile electrons rather than a spin-only chain, try
`build/bethe-hubbard-pbc 6 --u 4 --roots`; its nested Bethe ansatz has separate
charge momenta and spin rapidities. Use `bethe-hubbard-obc` for free ends;
see the [periodic](docs/hubbard.md) and [open-chain](docs/hubbard-open.md) guides.

For continuum bosons, try `build/bethe-lieb-liniger-pbc 4 --length 4 --c 1`.
Here 4 particles occupy a ring of physical length 4; it is not a four-site chain.

For three-state sites, try `build/bethe-su3-pbc 6 --roots`. This uses
`H=sum P`, where P swaps adjacent colors; its six-site energy is `-1-sqrt(13)`.

Next, distinguish the lowest state in a magnetization sector from a family
of excited states:

```sh
build/bethe-xxx-pbc 16 --sz 1
build/bethe-xxx-pbc 16 --sectors
build/bethe-xxx-pbc 16 --excitations 10 --spin 1
build/bethe-xxz-pbc 16 --delta 0.5 --excitations all --sz 1
```

XXX excitation scans select total spin with `--spin`; XXZ scans select
magnetization with `--sz`. Both include the sector minimum and report gaps
relative to the global ground state.

The default arithmetic is fp64. Select `--precision long-double` or, in an
enabled build, `--precision fp128` for higher precision. Add `--roots` to inspect
rapidities and quantum numbers.
Reports are formatted on a terminal and plain when redirected; `--format plain`
makes that choice explicit. Every frontend can also save CSV, TSV, and JSON
with provenance and native-precision values; see the [output guide](docs/output.md).
Always check convergence: the equation residual
is not an energy-error bound. See [CLI controls and diagnostics](docs/command-line.md)
for precision, tolerances, output formats, and exit statuses.

## Go deeper

For the theory behind the calculations, Jean-Sébastien Caux's
[The Bethe Ansatz](https://integrability.org/) is an excellent companion:
pedagogical derivations, integrable models, Bethe equations, and the
classification of states and excitations. Our guides explain the particular
conventions and numerical coverage implemented here; [CITATIONS.md](CITATIONS.md)
links the specific literature and sections used by the code.

The guides keep the examples, equations, numerical conventions, and limitations
together. Read them in roughly this order, or go straight to your model:

1. [Building and using Uni20](docs/building.md) — reproducible builds,
   local development, high precision, and library integration.
2. [Command-line calculations and diagnostics](docs/command-line.md) — choose
   a calculation and interpret the result.
3. [Periodic XXX states and conventions](docs/xxx.md) — connect the C++ API
   to quantum numbers, rapidities, energy, and momentum.
4. [Open XXX chains](docs/open-chains.md) — free ends, reflected scattering,
   and the changed energy normalization.
5. [XXX real-root excitations](docs/excitations.md) — what `COUNT` and `all`
   enumerate, computational cost, gaps, and failed candidates.
6. [Single-spinon dispersion](docs/spinons.md) — finite odd-chain branches
   versus the thermodynamic curve.
7. [Periodic XXZ chains](docs/xxz.md) — anisotropy, scaled rapidities,
   sector minima, and the excitation window.
8. [Free-end XXZ chains](docs/xxz-open.md) — boundary reflection phases,
   standing waves, and open-chain excitations.
9. [Periodic Hubbard rings](docs/hubbard.md) — nested charge/spin equations,
   doping, attractive interactions, symmetry mappings, and weak/strong coupling limits.
10. [Free-end Hubbard chains](docs/hubbard-open.md) — reflected scattering,
    standing waves, and unrestricted ground-state sectors.
11. [Periodic Lieb–Liniger gas](docs/lieb-liniger.md) — continuum units,
    repulsive bosons, weak coupling, and finite excitation windows.
12. [Periodic SU(3) chain](docs/su3.md) — nested real-root seas, balanced
    singlets, permutation normalization, and the spin-1 ULS mapping.
13. [Periodic Gaudin–Yang gas](docs/gaudin-yang.md) — continuum fermions,
    charge/spin roots, sector restrictions, and weak/strong-coupling limits.
14. [Supersymmetric t–J chain](docs/tj.md) — projected electrons, nested hole
    roots, the XXX limit, and fermionic momentum conventions.
15. [Spin-1 Takhtajan–Babujian chain](docs/takhtajan-babujian.md) — complex
    two-strings, finite-size deviations, and bilinear–biquadratic normalization.
16. [Richardson pairing](docs/richardson.md) — single-particle levels, blocking,
    eigenvalue variables and ground-state continuation through root collisions.
17. [Rational central spin](docs/central-spin.md) — specified bath couplings,
    magnetization sectors, and field continuation including exactly zero field.
18. [SU(n) fermion gas](docs/su-fermions.md) — component populations, arbitrary
    nesting depth, finite-ring shell restrictions and continuum limits.

For possible additions, see the [model catalogue and development proposal](docs/models.md):
known integrable families, literature, implementation restrictions, and suggested priorities.

## Source and attribution

`include/bethe/` contains the scalar-templated library, `apps/` the thin CLI
front ends, and `tests/` the regression suite. Numerical code is separate
from presentation and any future Python bindings. The public CMake target
is `bethe::bethe`; model guides include C++ examples.

This is the successor to MPToolkit's `misc/heisenberg-energy.cpp`; no legacy
copy is kept here. It retains GPL-3.0-or-later licensing and attribution.
See [COPYING](COPYING) and [References and provenance](CITATIONS.md).
The [paper archive](papers/README.md) collects redistributable reference PDFs
and links to other editions; the papers retain their own licenses.
