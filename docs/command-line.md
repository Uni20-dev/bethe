# Command-line calculations and diagnostics

[Back to the overview](../README.md)

The programs share precision selection and report formatting; convergence
controls depend on the solver. Start with the periodic XXX example below,
then follow a model guide for its Hamiltonian, supported states and limitations.

- [Choose a calculation](#choose-a-calculation)
- [Help, references and option syntax](#help-references-and-option-syntax)
- [Select arithmetic precision](#select-arithmetic-precision)
- [Read and save the report](#read-and-save-the-report)
- [CPU time and convergence](#cpu-time-and-convergence)
- [Model-specific entry points](#model-specific-entry-points)

## Choose a calculation

Examples assume the executables are on `PATH`. From a build directory, you can
instead use `./bethe-xxx-pbc` etc.; see [running the programs](building.md#run-the-programs).
No particular source or build directory is required for a calculation.

```sh
bethe-xxx-pbc 16
bethe-xxx-pbc 64 --precision fp64 --tolerance 1e-12
bethe-xxx-pbc 6 --precision long-double --roots
bethe-xxx-pbc 15 --sz 1/2
bethe-xxx-pbc 16 --sectors
bethe-xxx-pbc 65 --spinons
bethe-xxx-pbc 32 --excitations 10 --spin 1
bethe-xxx-pbc 5 --quantum-numbers -1,1 --roots
# In a binary128-enabled build:
bethe-xxx-pbc 16 --precision fp128 --tolerance 1e-30
```

### XXX state-selection modes

The default is a ground state; odd lengths return one of the degenerate
$`S^z =1/2`$ representatives. The five explicit modes are mutually exclusive:

| Option | Calculation or selection |
| --- | --- |
| `--sz VALUE` | A sector minimum; accepts integers, fractions such as `-3/2`, or half-integer decimals |
| `--sectors` | One lowest-energy representative in each sector, including spin-reversed partners |
| `--spinons` | The one-spinon family for odd $`N \ge 3`$, not the full $`S^z =1/2`$ spectrum |
| `--quantum-numbers LIST` | A strictly increasing comma-separated list; an empty string selects the fully polarized state |
| `--excitations COUNT\|all` | A restricted real-root family; see the [XXX excitation guide](excitations.md) |

Add `--roots` to print the roots and exact Bethe quantum numbers.

For **spectral weights**, use the separate
[`bethe-xxx-structure-factor`](xxx-structure-factor.md) frontend on an even ring.
It always scans the complete supported two-spinon family, and reports its
contribution to the zz or raising-channel sum rules.
For an infinite-chain spectral density, use
[`bethe-xxx-structure-factor-thermo`](xxx-structure-factor-thermo.md), with
`--momentum Q` for a frequency cut or `--momentum-points N` for a heat-map grid.

XXZ uses the same common precision and output controls, but different mode
combinations: `--excitations` selects a family at fixed `--sz`, not `--spin`.
See [periodic XXZ excitations](xxz.md#real-root-excitations) or
[free-end XXZ excitations](xxz-open.md#real-root-excitations) before transferring
XXX commands directly to the anisotropic model.

<a id="uni20-help-and-option-parsing"></a>

## Help, references and option syntax

### Information without calculation

All frontends use Uni20's shared CLI and presentation module. These options
need no model parameters and do not start a calculation:

| Option | Output |
| --- | --- |
| `--help`, `-h` | Program banner, grouped options with defaults and constraints, examples and conventions |
| `--references` | Literature links and notes on the modes they support |
| `--version` | Bethe version and revision |
| `--build-info` | Version plus Uni20's configured compiler and dependency information |

Only `--references` includes a bibliography. Color and width follow Uni20's
terminal policy; redirected help is plain by default. No arguments print help
to stderr and exit 1.

References include links and applicability notes. Normal numerical output
contains no bibliography; [CITATIONS.md](../CITATIONS.md) collects the full
bibliography and explains the conventions.

### Parsing rules

Both `--u 4` and `--u=4` are accepted. Scalar options reject duplicate
occurrences; `--csv`, `--tsv`, and `--json` remain repeatable, with one path per
occurrence. A value such as `--csv=--help` or `--csv=--references` is a filename,
not an information request.

Real-valued parameters remain strings until precision selection, so putting
`--precision fp128` last cannot lose digits. Help/version/build information and
references are handled before required-option and value validation, without
opening exports; a missing value detected while gathering arguments can still
be an error.
If both `--help` and `--references` are given, ordinary help takes precedence.
Invalid arguments produce a concise stderr diagnostic rather than full help.

Existing positional arguments and scientific option names remain available;
options may also precede positional arguments. Integer counts and half-integers
are parsed exactly, with negative counts, trailing junk and overflow rejected.
An explicitly empty list remains distinct from an omitted option.
All frontends support typed tables and file exports; see the
[output guide](output.md) for table names, metadata, and numeric encoding.

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
decimal half-integer labels, and semantic convergence markers.
`--format auto` (the default) selects this report on a terminal and a plain
human-readable report when redirected. Use `--format pretty` to save a formatted
report, or `--format csv|tsv|json` for scripts. Add repeatable `--csv FILE`,
`--tsv FILE` or `--json FILE` exports independently of screen output.
All retain the selected type's full round-trip precision, including
binary128; no displayed values are narrowed to `double`.

The pretty report honors `UNI20_COLOR`, `NO_COLOR`, `UNI20_GLYPHS`, and
`UNI20_CHARSET`. For example, `UNI20_GLYPHS=ascii UNI20_COLOR=never` selects
ASCII table rules and status markers without ANSI color. Terminal width (or
`COLUMNS`) selects aligned tables or labeled records; numeric strings are never
split or truncated, so a single long value can still exceed a very narrow
terminal. Root and auxiliary rows link to state records through `state_id`;
excitation scans keep the ground reference and failed estimates separate
from ranked levels. See the [table schemas](output.md#named-tables).

## CPU time and convergence

### Parallel state scans

`bethe-xxx-structure-factor` and the `--excitations` modes of
`bethe-xxx-pbc`, `bethe-xxx-obc`, `bethe-xxz-pbc`, and `bethe-xxz-obc` accept
`--threads N` (positive integer, default 1). The limit belongs to a Uni20
oneTBB scheduler, not one operating-system thread per state. Other modes of
the four energy tools reject an explicit `--threads` option.

Each state is solved independently; ordered result collection, failure
diagnostics and moment sums are deterministic. Worker threads do not write
tables. Library users can also select a scheduler for the shared real-root
scanner used by other models; their frontends do not yet expose this option.
See the [XXX structure-factor guide](xxx-structure-factor.md#parallel-calculations)
for examples. Set BLAS/OpenMP thread limits to one to avoid nested threading.

### Timing

Output includes energy, momentum (periodic systems only), normalized equation
residual, convergence status, update count, and solver CPU time in seconds.
CPU time measures process CPU consumption during state construction and solving,
not elapsed wall time (parallel workers' CPU times add together);
for scans it covers the whole scan, including the ground reference and energy
ordering for excitations. Report formatting and output
are excluded. Human reports include it once; exports carry it in the table
summary (trailing comments for CSV/TSV). Very short runs may report zero at the clock's
resolution; an unavailable or wrapped CPU clock is reported as `unavailable`.

### Iteration budgets and exit status

For the finite-root solvers below, `--max-iterations` defaults to 10000
**per state**; zero evaluates only the initial guess. Other solvers have their
own controls: for example, the [exact sine-Gordon excited solver](sine-gordon-excited.md)
defaults to 32768 nonlinear updates per state. Check the frontend's `--help`.

The iteration count depends on the model:

- **XXX/XXZ:** accepted Newton updates, sharing one budget across any
  continuation stages. XXX and gapless nonnegative XXZ start from zero roots;
  supported negative-Delta XXZ uses a free-fermion hyperbolic seed, and massive
  open XXZ uses coupling continuation.
- **Hubbard:** accepted Newton updates across all continuation stages, from
  a large-U seed. The exact U=0 path needs no updates.
- **SU(3):** accepted Newton updates from a filled-sea density seed.
- **Gaudin–Yang:** accepted updates across strong-to-weak continuation.
  The exact free path needs no updates.

Exit codes distinguish numerical failure from invalid input:

| Code | Meaning |
| --- | --- |
| `0` | Successful calculation; all requested states converged |
| `2` | Incomplete or nonconverged calculation, for example an exhausted budget or stalled line search |
| `1` | Invalid input or another error |

A nonconverged estimate is explicitly marked; some solvers withhold failed
quantities entirely. There is no silent precision fallback.

### Residuals and tolerances

The residual measures how closely the rapidities satisfy the Bethe equations;
it is not a bound on the error in the energy. The finite-root solvers described
here default to 32 times the selected type's epsilon; integral-equation and
NLIE solvers have their own error budgets and refinement checks.

The logarithmic spin-chain solvers normalize
by N (periodic) or 2N (open). Supported negative-Delta XXZ instead uses
[rank-subtracted, scaled equations](xxz-negative.md); massive open XXZ includes
a [regularized boundary residual](xxz-open-massive.md) when needed. These
conventions are explicitly reported and are not interchangeable.

Hubbard normalizes both
equation families by L for periodic rings or 2(L+1) for free ends, and reports
their maximum. Even a stopped continuation reports residuals at the requested
root-sector U.

Lieb–Liniger uses a component-scaled dimensionless residual,
with a weak-coupling scale that resolves roots of order sqrt(c*ell);
see its [numerical-method guide](lieb-liniger.md#numerical-method-and-precision).

SU(3) normalizes both nested equation families by L and reports the maximum
over its independent reflection-reduced equations, as explained in the
[SU(3) guide](su3.md#numerical-method-and-failure-reporting).

Gaudin–Yang scales weak-coupling charge residuals by their root or sqrt(c*ell)
and otherwise uses N; see its [diagnostics guide](gaudin-yang.md#numerical-method-and-diagnostics).

See the model guides for the equations. Increasing precision can help resolve
closely spaced levels, but always inspect the convergence status before
treating an energy as an eigenvalue.

## Model-specific entry points

These notes highlight differences from the XXX example, not every available
program. See the [model catalogue](models.md#candidate-index) and
[overview](../README.md) for additional models and thermodynamic tools.

### Spin-1/2 chains

`bethe-xxx-pbc` is the periodic XXX front end. For free-end OBC, use the separate
`bethe-xxx-obc` program described in the [open-chain guide](open-chains.md).

For anisotropy, use [`bethe-xxz-pbc`](xxz.md) or the free-end
[`bethe-xxz-obc`](xxz-open.md); both require `--delta`.
Both support ground states and sector minima for finite $`\Delta \ge 0`$;
their excitation and specified-real-root modes still require $`0\le \Delta \le 1`$.
Massive free-end ground states report an explicit boundary-root coordinate
where needed; see the [boundary-root guide](xxz-open-massive.md).

### Hubbard electrons

For electrons, [`bethe-hubbard-pbc`](hubbard.md) requires `--u` and supports
either sign of the interaction on even rings. Select `--particles N` and
`--sz VALUE`; the defaults remain N=L and Sz=0. See the
[supported sector families](hubbard-sectors.md) before choosing a doped sector.

For free ends, [`bethe-hubbard-obc`](hubbard-open.md) supports every physical
N and Sz, odd or even L, and either sign of U. It defaults to N=L and the
smallest nonnegative Sz compatible with N (0 or 1/2).

### Continuum gases

For continuum bosons, [`bethe-lieb-liniger-pbc`](lieb-liniger.md) takes the
particle count as its positional argument and requires `--length ELL --c C`.
Its excitation scans also require a finite `--padding P` window; continuum
momentum is not reduced to a Brillouin zone.

For spin-1/2 continuum fermions, [`bethe-gaudin-yang-pbc`](gaudin-yang.md)
takes N and requires `--length ELL --c C`. Select the populations with `--sz`;
interacting mixed-spin sectors require odd populations of both spins. The
free and fully polarized limits also support other populations.

For more than two continuum fermion components,
[`bethe-sun-fermions-pbc`](su-fermions.md) takes
`--populations N1,N2,... --length ELL --c C`. Every occupied interacting
population must be odd; free and single-component limits allow any counts.
`--roots` displays all nested seas, with the physical-to-nesting component
mapping retained in the report.

### Supersymmetric t–J chain

For projected lattice electrons, [`bethe-tj-pbc`](tj.md) fixes t=1, J=2 and
takes L with optional `--particles N --sz VALUE`. It supports doped sectors
with odd populations of both spins, plus every no-hole or fully polarized
sector. Defaults are N=L and the smallest nonnegative Sz compatible with N.

### SU(3)/ULS and spin-1 TB chains

For three-state sites, [`bethe-su3-pbc`](su3.md) takes L divisible by three
and selects the balanced SU(3) singlet ground state of `H=sum P`. It reports
two nested rapidity families; arbitrary sectors and finite-ring excitations are
not yet available in that frontend. [`bethe-su3-dispersion`](su3-dispersion.md)
provides thermodynamic elementary lines and two-/four-soliton continuum bounds,
with optional three-site-cell folding.

For the spin-1 bilinear–biquadratic TB point,
[`bethe-tb-pbc`](takhtajan-babujian.md) takes even L>=4 and selects the
zero-field singlet ground state of `H=sum[S.S-(S.S)^2]`. Its `--roots` report
retains finite deviations and the real and imaginary parts of each rapidity.
[`bethe-tb-dispersion`](tb-dispersion.md) instead supplies infinite-chain
spinon lines and two-/four-spinon bounds, with explicit spin labels and
optional two-site momentum folding.

### Biquadratic chains and TL modules

For the pure spin-1 biquadratic chain,
[`bethe-biquadratic-obc`](biquadratic.md) takes even N>=2 and selects the
free-end singlet ground state of `H=-sum(S.S)^2` by default. `--through-lines`
selects a TL module, `--sectors` lists module minima, and `--excitations COUNT|all`
scans the restricted real-root family (default through-lines=2). The output
includes physical multiplicities, not a decomposition into SU(2) multiplets.

Its `--roots` output belongs to the auxiliary XXZ model with opposite end
fields, not the physical spin-1 chain. The [Q-system](xxz-open-qsystem.md)
and [two-string singlet](xxz-open-two-string.md) modes include selected complex roots.

With `--ferromagnetic`, the [excitation tools](biquadratic-ferromagnetic.md)
also offer an analytic one-defect band and targeted bound pairs/triples.
Ferro real-root scans accept odd/even N; `--real-window WIDTH` restricts
`--excitations` to a [high-label scattering window](biquadratic-scattering.md).

For [mixed bound-pair-plus-defect states](biquadratic-pair-defect.md), use
`--pair-defect I1,...,Ir,J` or `--pair-defects COUNT|all --real-defects R`
(default R=1). The latter can restrict all label ranges with `--mixed-window WIDTH`.
Both fix ell=N-2(R+2); these are not full module spectra, and partial scans return exit 2.
Explicit `--quantum-numbers` also accepts odd N; AF ground/sector and Q-system
modes retain their even-N restriction.

### Pairing and central spin

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

### Integrable ladder

For Wang's integrable ladder, [`bethe-ladder-pbc`](ladder.md) takes the
number of rungs L>=2 and requires `--rung JR`. The leg coefficient is 1
and the four-spin coefficient is fixed at 4; this is not an ordinary
Heisenberg ladder. Use `--singlets NS` for one singlet-count sector or
`--sectors` for all sector minima; triplet populations are minimized too.
`--roots` prints the selected state's highest-weight representative and
retains the distinction between it and a physical SU(4) descendant.

### Sine-Gordon field theory

For [sine-Gordon excitations](sine-gordon-excitations.md),
`bethe-sine-gordon-dispersion --p P` supplies particle lines and pair thresholds;
`bethe-sine-gordon-bethe-yang --p P --length L` supplies asymptotic same-charge
pair levels for P>=1. The latter omits wrapping corrections and is distinct
from the [finite-volume vacuum](sine-gordon.md) calculation.

`bethe-sine-gordon-excited --p P --length L --number 0.5` instead solves the
[exact two-hole NLIE](sine-gordon-excited.md) for P>=1, I=0.5 or 1.5. Its
`gap` table subtracts a separately converged vacuum; use
`--csv-table gap=gap.csv` to export that table directly.
