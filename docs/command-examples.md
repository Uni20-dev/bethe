# Command examples by model

[Supported models and scope](models.md) · [Shared syntax and precision](command-line.md) ·
[Plotted tutorials](tutorials/index.md)

These are entry commands, not a complete option or state-family catalogue.
Use each model guide and the executable's `--help` for supported regimes and
valid table names. Commands assume installed executables on `PATH`.

## Model-specific entry points

These notes highlight differences from the XXX example, not every available
program. See the [model catalogue](models.md#current-coverage) and
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
Explicit `--precision bits-N` or `digits-N` is supported with Uni20 MPFR scalars, including all three existing branches and root exports.

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
The frontend supports explicit `bits-N`/`digits-N` precision with Uni20 MPFR scalars, including fixed singlet and magnetization sectors and partial scans.

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

`bethe-kondo-response` also accepts `bits-N` and `digits-N` with Uni20 MPFR
scalars enabled. Its series and quadrature budgets remain separate from the
arithmetic precision; see [Kondo controls](kondo.md#library-use-and-numerical-contract).

Hubbard thermodynamic MPFR uses explicit precision for quadrature and extrema.
All paths need MPFR scalars; doped dispersion backgrounds use normal Uni20
CPU/provider LU dispatch. See [Hubbard dispersions](hubbard-dispersion.md#explicit-precision-in-c).

`bethe-tasep-pbc` accepts `bits-N` and `digits-N` with Uni20 MPFR/MPC
scalars, using normal CPU/provider LU dispatch. Seed and Newton work budgets
remain independent of arithmetic precision; see [TASEP](tasep.md#explicit-precision).

`bethe-asep-pbc` also accepts explicit MPFR precision for its TASEP seed and
continuation; see [ASEP controls and shared seeds](asep.md#explicit-precision).
