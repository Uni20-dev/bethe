# Bethe ansatz solvers

C++23 tools for finite-system Bethe ansatz calculations, complementing
[Uni20](https://github.com/Uni20-dev/uni20) and the
[Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

A calculation starts with a model, a boundary condition, and a choice of
Bethe quantum numbers. The solver finds the corresponding rapidities and
uses them to calculate energies and, for periodic chains, momenta. You can
begin with a ground state without choosing the quantum numbers yourself,
then explore magnetization sectors, excitation families, and spinons.

## Choose a model

- [Periodic XXX](docs/xxx.md): `heisenberg-energy` gives ground states,
  sector minima, real-root excitations, and the odd-chain one-spinon branch.
- [Free-end XXX](docs/open-chains.md): `heisenberg-open-energy` supports
  ground states, sector minima, and real-root excitations, without lattice momentum.
- [Periodic XXZ](docs/xxz.md): `xxz-energy` supports ground states, sector minima,
  and an anisotropy-dependent real-root excitation family for `0 <= Delta <= 1`.

The finite-chain models use spin-1/2 operators, J=1, and zero magnetic field.
Excitation scans cover explicitly supported real-root families, **not complete
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
build/heisenberg-energy 4
```

Its total energy is E=-2 in our spin-1/2 normalization. The report includes
energy per site, convergence diagnostics, and solver CPU time. To see how the
boundary condition or anisotropy changes the problem, try:

```sh
build/heisenberg-open-energy 4
build/xxz-energy 4 --delta 0.5
```

Next, distinguish the lowest state in a magnetization sector from a family
of excited states:

```sh
build/heisenberg-energy 16 --sz 1
build/heisenberg-energy 16 --sectors
build/heisenberg-energy 16 --excitations 10 --spin 1
build/xxz-energy 16 --delta 0.5 --excitations all --sz 1
```

XXX excitation scans select total spin with `--spin`; XXZ scans select
magnetization with `--sz`. Both include the sector minimum and report gaps
relative to the global ground state.

The default arithmetic is fp64. Add `--roots` to inspect rapidities and
quantum numbers, or `--precision long-double` for the platform's extended type.
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

## Source and attribution

`include/bethe/` contains the scalar-templated library, `apps/` the thin CLI
front ends, and `tests/` the regression suite. Numerical code is separate
from presentation and any future Python bindings. The public CMake target
is `bethe::bethe`; model guides include C++ examples.

This is the successor to MPToolkit's `misc/heisenberg-energy.cpp`; no legacy
copy is kept here. It retains GPL-3.0-or-later licensing and attribution.
See [COPYING](COPYING) and [References and provenance](CITATIONS.md).
