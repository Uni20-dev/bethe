# Bethe ansatz solvers

C++23 tools for finite-system Bethe ansatz calculations, complementing
[Uni20](https://github.com/Uni20-dev/uni20) and the
[Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

High precision is a first-class feature: all finite-chain solvers support
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

- [Periodic XXX](docs/xxx.md): `bethe-xxx-pbc` gives ground states,
  sector minima, real-root excitations, and the odd-chain one-spinon branch.
- [Free-end XXX](docs/open-chains.md): `bethe-xxx-obc` supports
  ground states, sector minima, and real-root excitations, without lattice momentum.
- [Periodic XXZ](docs/xxz.md): `bethe-xxz-pbc` supports ground states, sector minima,
  and an anisotropy-dependent real-root excitation family for `0 <= Delta <= 1`.
- [Free-end XXZ](docs/xxz-open.md): `bethe-xxz-obc` provides the corresponding
  open-chain calculations with no boundary fields or lattice momentum.
- [Periodic Hubbard](docs/hubbard.md): `bethe-hubbard-pbc` gives repulsive
  half-filled spin sectors, balanced attractive ground states at any even filling,
  and selected doped sectors on even rings, plus the unrestricted U=0 limit.
- [Free-end Hubbard](docs/hubbard-open.md): `bethe-hubbard-obc` gives ground
  states at every physical filling and spin projection, for either sign of U
  and odd or even lengths, without lattice momentum.

The XXX and XXZ models use spin-1/2 operators, J=1, and zero magnetic field.
Hubbard uses hopping t=1 and the unshifted interaction `U*n_up*n_down`.
Spin-chain excitation scans cover explicitly supported real-root families, **not complete
spectra**: `--excitations all` means all states in that family. Complex strings
and infinite-root descendants are not implemented.

[Analytic thermodynamic spinon dispersions](docs/spinons.md#thermodynamic-dispersion)
are also available for XXX and gapless XXZ; these are distinct from finite-chain
calculations.

## Build

You need CMake 3.24+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
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
makes that choice explicit. Always check convergence: the equation residual
is not an energy-error bound. See [CLI controls and diagnostics](docs/command-line.md)
for precision, tolerances, output formats, and exit statuses.

## Go deeper

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

## Source and attribution

`include/bethe/` contains the scalar-templated library, `apps/` the thin CLI
front ends, and `tests/` the regression suite. Numerical code is separate
from presentation and any future Python bindings. The public CMake target
is `bethe::bethe`; model guides include C++ examples.

This is the successor to MPToolkit's `misc/heisenberg-energy.cpp`; no legacy
copy is kept here. It retains GPL-3.0-or-later licensing and attribution.
See [COPYING](COPYING) and [References and provenance](CITATIONS.md).
