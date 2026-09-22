# Command-line calculations and diagnostics

[Back to the overview](../README.md)

The programs share precision selection, convergence controls, and
report formatting. The periodic XXX chain is a useful first example; model
and state-selection details live in the linked guides.

`bethe-xxx-pbc` is the periodic XXX front end. For free-end OBC, use the separate
`bethe-xxx-obc` program described in the [open-chain guide](open-chains.md).
For anisotropy, use [`bethe-xxz-pbc`](xxz.md) or the free-end
[`bethe-xxz-obc`](xxz-open.md); both require `--delta`.
For electrons, [`bethe-hubbard-pbc`](hubbard.md) requires `--u` and supports
either sign of the interaction on even rings. Select `--particles N` and
`--sz VALUE`; the defaults remain N=L and Sz=0. See the
[supported sector families](hubbard-sectors.md) before choosing a doped sector.
For free ends, [`bethe-hubbard-obc`](hubbard-open.md) supports every physical
N and Sz, odd or even L, and either sign of U. It defaults to N=L and the
smallest nonnegative Sz compatible with N (0 or 1/2).
For continuum bosons, [`bethe-lieb-liniger-pbc`](lieb-liniger.md) takes the
particle count as its positional argument and requires `--length ELL --c C`.
Its excitation scans also require a finite `--padding P` window; continuum
momentum is not reduced to a Brillouin zone.
For three-state sites, [`bethe-su3-pbc`](su3.md) takes L divisible by three
and selects the balanced SU(3) singlet ground state of `H=sum P`. It reports
two nested rapidity families; arbitrary sectors and excitations are not yet
available for this model.
For spin-1/2 continuum fermions, [`bethe-gaudin-yang-pbc`](gaudin-yang.md)
takes N and requires `--length ELL --c C`. Select the populations with `--sz`;
interacting mixed-spin sectors require odd populations of both spins. The
free and fully polarized limits also support other populations.

Each program's `--help` (and no-argument usage) includes relevant literature
references with links and a note on the modes they support. Normal numerical
output is unchanged. [CITATIONS.md](../CITATIONS.md) gives the full bibliography
and explains the conventions used here.

## Choose a calculation

```sh
build/bethe-xxx-pbc 16
build/bethe-xxx-pbc 64 --precision fp64 --tolerance 1e-12
build/bethe-xxx-pbc 6 --precision long-double --roots
build/bethe-xxx-pbc 15 --sz 1/2
build/bethe-xxx-pbc 16 --sectors
build/bethe-xxx-pbc 65 --spinons
build/bethe-xxx-pbc 32 --excitations 10 --spin 1
build/bethe-xxx-pbc 5 --quantum-numbers -1,1 --roots
# In a binary128-enabled build:
build/bethe-xxx-pbc 16 --precision fp128 --tolerance 1e-30
```

The default is a ground state; odd lengths return one of the degenerate
`Sz=1/2` representatives. `--sz` accepts integers, fractions such as `-3/2`,
or half-integer decimals. `--sectors` reports one lowest-energy representative
in each sector, including spin-reversed partners. `--spinons` requires odd
`N>=3` and reports the one-spinon family, not the full `Sz=1/2` spectrum.
`--quantum-numbers` takes a strictly increasing comma-separated list;
an empty string selects the fully polarized state. `--excitations COUNT|all`
scans a restricted real-root family, described in the [XXX excitation guide](excitations.md).
These five explicit XXX modes are mutually exclusive. `--roots` also prints
the exact Bethe quantum numbers.

XXZ uses the same common precision and output controls, but different mode
combinations: `--excitations` selects a family at fixed `--sz`, not `--spin`.
See [periodic XXZ excitations](xxz.md#real-root-excitations) or
[free-end XXZ excitations](xxz-open.md#real-root-excitations) before transferring
XXX commands directly to the anisotropic model.

## Select arithmetic precision

`--precision` selects `fp64` (the default), `long-double`, or optional `fp128`.
`long double` precision is platform-dependent; `fp128` uses Uni20's configured
MPLAPACK binary128 type. The fp128 CLI currently requires MPLAPACK's native
`_Float128/strfromf128` mode: other modes are rejected at compile time because
the pinned Uni20's scalar I/O narrows them to `long double`.
Calculations, tolerance parsing, and output retain the
selected precision; displayed digits are not a guarantee of energy accuracy.

For setup details and the supported MPLAPACK configuration, see
[building with binary128](building.md#enable-binary128).

## Read and save the report

Output uses Uni20's presentation layer on terminals: aligned fields and tables,
exact fractional quantum-number labels, and semantic convergence markers.
`--format auto` (the default) selects this report on a terminal and the existing
plain, whitespace-separated output when redirected. Use `--format pretty` to
save a formatted report or `--format plain` to request script-oriented output
explicitly. Both retain the selected type's full round-trip precision, including
binary128; no displayed values are narrowed to `double`.

The pretty report honors `UNI20_COLOR`, `NO_COLOR`, `UNI20_GLYPHS`, and
`UNI20_CHARSET`. For example, `UNI20_GLYPHS=ascii UNI20_COLOR=never` selects
ASCII table rules and status markers without ANSI color. Terminal width (or
`COLUMNS`) selects aligned tables or labeled records; numeric strings are never
split or truncated, so a single long value can still exceed a very narrow
terminal. Scans separate energies from convergence diagnostics, and root tables
identify their sector, hole, or excitation level.

## CPU time and convergence

Output includes energy, momentum (periodic systems only), normalized equation
residual, convergence status, update count, and solver CPU time in seconds.
CPU time measures process CPU consumption during state construction and solving,
not elapsed wall time;
for scans it covers the whole scan, including the ground reference and energy
ordering for excitations. Report formatting and output
are excluded. Both output formats include it, in metadata before scan tables
(a `# CPU time:` comment for spin-chain plain scans). Very short runs may report zero at the clock's
resolution; an unavailable or wrapped CPU clock is reported as `unavailable`.

`--max-iterations` defaults to 10000 **per state**;
zero evaluates only the initial guess (zero roots for the XXX/XXZ solvers).
Hubbard counts accepted Newton updates across all continuation stages and
uses a large-U seed; its exact U=0 path needs no updates.
SU(3) counts accepted Newton updates from a filled-sea density seed.
Gaudin–Yang counts accepted updates across its strong-to-weak continuation;
its exact free path needs no updates.
Exit status is 0 if all states converged, 2 if any exhausted their budget
or a Newton line search stalled, and 1 for invalid input or another
error. A nonconverged result is explicitly marked as an unconverged
estimate. There is no silent precision fallback. Run `--help` for the options.

The residual measures how closely the rapidities satisfy the Bethe equations;
it is not a bound on the error in the energy. The default tolerance is 32
times the selected type's epsilon. Periodic spin-chain solvers normalize by N,
while the open spin-chain solvers normalize by 2N. Hubbard normalizes both
equation families by L for periodic rings or 2(L+1) for free ends, and reports
their maximum. Even a stopped continuation reports residuals at the requested
root-sector U. Lieb–Liniger uses a component-scaled dimensionless residual,
with a weak-coupling scale that resolves roots of order sqrt(c*ell);
see its [numerical-method guide](lieb-liniger.md#numerical-method-and-precision).
SU(3) normalizes both nested equation families by L and reports the maximum
over its independent reflection-reduced equations, as explained in the
[SU(3) guide](su3.md#numerical-method-and-failure-reporting).
Gaudin–Yang scales weak-coupling charge residuals by their root or sqrt(c*ell)
and otherwise uses N; see its [diagnostics guide](gaudin-yang.md#numerical-method-and-diagnostics).
See the model guides for the equations.
Increasing precision can help resolve closely spaced levels, but always
inspect the convergence status before treating an energy as an eigenvalue.

Continue with [periodic XXX conventions](xxx.md), [free-end chains](open-chains.md),
XXZ with [periodic](xxz.md) or [free-end](xxz-open.md) boundaries, or Hubbard
with [periodic](hubbard.md) or [free-end](hubbard-open.md) boundaries.
