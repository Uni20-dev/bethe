# Command-line calculations and diagnostics

[Back to the overview](../README.md)

The programs share precision selection, convergence controls, and
report formatting. The periodic XXX chain is a useful first example; model
and state-selection details live in the linked guides.

`bethe-xxx-pbc` is the periodic XXX front end. For free-end OBC, use the separate
`bethe-xxx-obc` program described in the [open-chain guide](open-chains.md).
For anisotropy, use [`bethe-xxz-pbc`](xxz.md) or the free-end
[`bethe-xxz-obc`](xxz-open.md); both require `--delta`.
Both support ground states and sector minima for finite `Delta>=0`;
their excitation and specified-real-root modes still require `0<=Delta<=1`.
Massive free-end ground states report an explicit boundary-root coordinate
where needed; see the [boundary-root guide](xxz-open-massive.md).
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
For projected lattice electrons, [`bethe-tj-pbc`](tj.md) fixes t=1, J=2 and
takes L with optional `--particles N --sz VALUE`. It supports doped sectors
with odd populations of both spins, plus every no-hole or fully polarized
sector. Defaults are N=L and the smallest nonnegative Sz compatible with N.
For the spin-1 bilinear–biquadratic TB point,
[`bethe-tb-pbc`](takhtajan-babujian.md) takes even L>=4 and selects the
zero-field singlet ground state of `H=sum[S.S-(S.S)^2]`. Its `--roots` report
retains finite deviations and the real and imaginary parts of each rapidity.
For the pure spin-1 biquadratic chain,
[`bethe-biquadratic-obc`](biquadratic.md) takes even N>=2 and selects the
free-end singlet ground state of `H=-sum(S.S)^2` by default. `--through-lines`
selects a TL module, `--sectors` lists module minima, and `--excitations COUNT|all`
scans the restricted real-root family (default through-lines=2). The output
includes physical multiplicities, not a decomposition into SU(2) multiplets.
Its `--roots` output belongs to the auxiliary XXZ model with opposite end
fields, not the physical spin-1 chain. Complex-root levels are not yet included.
For reduced BCS pairing, [`bethe-richardson`](richardson.md) requires
`--levels E0,E1,... --pairs M --g G`, with optional zero-based `--blocked`
indices. Levels are single-particle energies; g>=0 is attractive. This is
not a chain, so no boundary suffix or lattice momentum applies. Use
`--variables` for regularized eigenvalue variables, not pair rapidities.
For the rational central-spin model, [`bethe-central-spin`](central-spin.md)
requires `--couplings A1,A2,... --field B --sz SZ`. All spins are 1/2 and Sz
includes the central spin. Distinct nonzero couplings and the field can
have either sign; exactly zero field is supported. This nonspatial model
also has no boundary suffix; `--variables` prints regularized variables.
For more than two continuum fermion components,
[`bethe-sun-fermions-pbc`](su-fermions.md) takes
`--populations N1,N2,... --length ELL --c C`. Every occupied interacting
population must be odd; free and single-component limits allow any counts.
`--roots` displays all nested seas, with the physical-to-nesting component
mapping retained in the report.
For Wang's integrable ladder, [`bethe-ladder-pbc`](ladder.md) takes the
number of rungs L>=2 and requires `--rung JR`. The leg coefficient is 1
and the four-spin coefficient is fixed at 4; this is not an ordinary
Heisenberg ladder. Use `--singlets NS` for one singlet-count sector or
`--sectors` for all sector minima; triplet populations are minimized too.
`--roots` prints the selected state's highest-weight representative and
retains the distinction between it and a physical SU(4) descendant.

Each program's `--help` (and no-argument usage) includes relevant literature
references with links and a note on the modes they support. Normal numerical
output is unchanged. [CITATIONS.md](../CITATIONS.md) gives the full bibliography
and explains the conventions used here.

### Uni20 help and option parsing

`bethe-hubbard-dispersion` is the first frontend using Uni20's shared CLI and
presentation module. `--help` (or `-h`) now shows a styled program banner,
grouped options with their actual defaults and constraints, examples,
conventions, and the same application-owned literature references. Color and
width follow Uni20's terminal policy; redirected help is plain by default.
`--version` gives the Bethe version/revision and `--build-info` adds Uni20's
configured compiler and dependency information. Neither includes a bibliography
or starts a calculation. No arguments still print help to stderr and exit 1.

Both `--u 4` and `--u=4` are accepted. Scalar options reject duplicate
occurrences; `--csv`, `--tsv`, and `--json` remain repeatable, with one path per
occurrence. A value such as `--csv=--help` is a filename, not a help request.
Real-valued parameters remain strings until precision selection, so putting
`--precision fp128` last cannot lose digits. Help/version/build information is
handled before required-option and value validation, without opening exports;
a missing value detected while gathering arguments can still be an error.
Invalid arguments produce a concise stderr diagnostic rather than full help.

The other frontends retain their existing parsers for now. The
[Hubbard dispersion guide](hubbard-dispersion.md) and [output guide](output.md)
describe the pilot's numerical and file-output conventions, which are unchanged.

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
zero evaluates only the initial guess. XXX and gapless nonnegative XXZ start
from zero roots. Supported negative-Delta XXZ uses a free-fermion hyperbolic seed,
and massive open XXZ uses coupling continuation; these count accepted Newton
updates, with one shared budget across any continuation stages.
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
times the selected type's epsilon. The logarithmic spin-chain solvers normalize
by N (periodic) or 2N (open). Supported negative-Delta XXZ instead uses
[rank-subtracted, scaled equations](xxz-negative.md); massive open XXZ includes
a [regularized boundary residual](xxz-open-massive.md) when needed. These
conventions are explicitly reported and are not interchangeable. Hubbard normalizes both
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
