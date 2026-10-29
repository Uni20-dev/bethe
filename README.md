# Bethe ansatz solvers

C++23 tools and a scalar-templated library for finite-system and thermodynamic
Bethe ansatz calculations, complementing [Uni20](https://github.com/Uni20-dev/uni20)
and the [Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

**[Documentation and tutorials](https://uni20-dev.github.io/bethe/)** ·
[Supported models](docs/models.md) · [Build and install](docs/building.md)

Calculate ground states, selected excitation families, dispersions, continuum
edges and supported spectral weights. The [model guides](docs/models.md) state
which boundaries, sectors and root families each calculation covers.

High precision is a first-class feature: **fp64** by default, platform-native
**long-double**, **fp80** where extended80 is native, optional **fp128**
(binary128, about 34 decimal digits), and explicit **arbitrary precision**
through MPFR/MPC in every numerical frontend. Run `bethe-info` to see what your
build provides. [Precision controls](docs/command-line.md#select-arithmetic-precision)
explain how arithmetic precision differs from numerical convergence.

## First calculation

After [building or installing](docs/building.md), put the programs on `PATH`
and compare periodic and free-end four-site Heisenberg chains:

```sh
bethe-xxx-pbc 4
bethe-xxx-obc 4
```

The periodic chain has total energy $`E=-2`$ for spin-$`1/2`$ operators and
$`J=1`$. The reports include energy per site, convergence diagnostics and CPU
time. The [first-calculation guide](docs/getting-started.md) explains these
quantities, boundaries, quantum numbers, precision and exports.

Without `PATH` setup, use `./bethe-xxx-pbc` from the build directory.
Python is only needed for the optional plotting and maintenance scripts.

## Choose a route

- **Run calculations:** [first steps](docs/getting-started.md),
  [choose a supported model](docs/models.md), then
  [save and plot results](docs/output.md).
- **Compare with MPS or iMPS:** start with a
  [plotted tutorial](docs/tutorials/index.md). Each explains energy/momentum
  conventions and the physical contribution being calculated.
- **Use C++:** [embed `bethe::bethe`](docs/building.md#embed-the-library),
  use the synchronous APIs in the model guides, then compose independent
  work with [Async](docs/async-execution.md) where useful.

## Tutorials with plots

- [XXZ spinons and continua](docs/tutorials/xxz-spinons.md): gapless/gapped
  chains and two-site iMPS momentum conventions.
- [Hubbard at half filling](docs/tutorials/hubbard-half-filled.md) and
  [away from half filling](docs/tutorials/hubbard-doped.md): spin/charge lines
  and the choice of energy zero.
- [XXX spectral weights](docs/tutorials/xxx-structure-factor.md): finite-ring
  two-spinon intensities, sum rules and broadening.
- [Why high precision matters](docs/tutorials/xxx-precision.md): small energy
  scales, stable formulas and arbitrary precision.

The [full tutorial collection](docs/tutorials/index.md) includes source data,
plotting scripts and examples for every implemented model family.

## Theory, development and attribution

Jean-Sébastien Caux's [The Bethe Ansatz](https://integrability.org/) is a
pedagogical companion for the underlying theory. [CITATIONS.md](CITATIONS.md)
lists the literature used by the implementation; the
[model roadmap](docs/model-roadmap.md) surveys possible extensions.

`include/bethe/` contains the library, `apps/` the command-line frontends, and
`tests/` the regression suite. Numerical code is separate from presentation
and any future Python bindings.

Bethe succeeds MPToolkit's `misc/heisenberg-energy.cpp` and retains
GPL-3.0-or-later licensing and attribution. See [COPYING](COPYING),
[references and provenance](CITATIONS.md), and the
[redistributable paper archive](papers/README.md).
