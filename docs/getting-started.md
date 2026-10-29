# Your first calculation

[Build and install](building.md) · [Supported models](models.md) ·
[Plotted tutorials](tutorials/index.md)

This walkthrough connects a small calculation to the quantities in its report.
Use installed commands on `PATH`, or run them with `./` from your build directory.
No Python is needed. First inspect the available arithmetic and providers:

```sh
bethe-info
```

## Ground-state energy and conventions

```sh
bethe-xxx-pbc 4
```

This is the periodic spin-$`1/2`$ Heisenberg chain, with four sites and
$`H=\sum_j\mathbf S_j\cdot\mathbf S_{j+1}`$, $`J=1`$. Its ground-state
**total energy is $`-2`$**, and its energy per site is $`-1/2`$.
These use spin operators, not Pauli matrices; replacing each spin by a Pauli
matrix multiplies the exchange energy by four.

Read the energy together with the convergence status, equation residual and
iteration count. The residual measures how well the Bethe equations are
satisfied; it is not a bound on the energy error. The report also gives CPU
time. See [convergence diagnostics](command-line.md#cpu-time-and-convergence)
for timing, budgets and exit statuses.

## Change the boundaries or model

```sh
bethe-xxx-obc 4
bethe-xxz-pbc 4 --delta 0.5
```

The first command opens the chain, removing the bond between sites four and
one. Open-chain standing waves do not carry a lattice momentum. The second
keeps periodic boundaries but changes the anisotropy: the longitudinal
$`S^zS^z`$ coupling is half the transverse exchange.

Commands usually follow `bethe-<model>-<boundary>`: `pbc` is periodic and
`obc` is open. Open boundaries mean the specific boundaries in that model's
guide; for example, the quantum-group XXZ frontend has imaginary end fields,
whereas `bethe-xxz-obc` has free ends. Thermodynamic `-dispersion` tools solve
a different problem from finite-chain energies.

## Select a sector or excitation family

```sh
bethe-xxx-pbc 16 --sz 1
bethe-xxx-pbc 16 --excitations 10 --spin 1
```

The first command finds the minimum in a specified magnetization sector.
The second enumerates a supported real-root excitation family with total
spin one, including its minimum, and reports gaps relative to the global
ground state. It does not enumerate the whole spectrum. In a scan,
`--excitations all` means all states in the documented family or window.
For XXZ, magnetization `--sz` remains a good quantum number, while total spin
generally does not. [XXX states and excitations](xxx-family.md) explains the
labels and restrictions before larger scans.

## Inspect roots and save results

```sh
bethe-xxx-pbc 4 --roots
bethe-xxx-pbc 4 --csv-table states=xxx-four-site.csv
```

`--roots` shows rapidities and quantum numbers on screen. A table export writes
that table independently of screen output. Table names belong to the frontend;
its `--help` lists the valid choices. CSV/TSV carry provenance comments before
the header, and JSON carries structured metadata. The [output guide](output.md)
shows how to read them in plotting programs.

Use `--help` for options and `--references` for the literature used by a tool.
Commands run without numerical arguments show help rather than guessing a
system size.

## Increase precision deliberately

The default is binary64. In suitable configured builds:

```sh
bethe-xxx-pbc 4 --precision long-double
bethe-xxx-pbc 4 --precision fp128
bethe-xxx-pbc 4 --precision digits-100
```

The last command requests explicit MPFR arithmetic. Its availability depends
on the [build configuration](building.md#enable-arbitrary-precision).
`bits-128` means 128 significand bits; it is different from `fp128`, which has
113. More bits do not repair a missing state family, an insufficient quadrature
mesh or an iteration limit. Compare converged calculations with suitable work
budgets, and consult the [precision tutorial](tutorials/xxx-precision.md).

## Next steps

Choose a [supported model](models.md), or start a plotted tutorial on
[XXZ spinons](tutorials/xxz-spinons.md),
[Hubbard spin and charge excitations](tutorials/hubbard-half-filled.md), or
[XXX spectral weights](tutorials/xxx-structure-factor.md). For MPS comparisons,
match the Hamiltonian normalization, momentum convention, energy reference
and operator normalization before comparing curves.
