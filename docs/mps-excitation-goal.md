# Goal: exact excitation benchmarks for MPS and iMPS

Develop native-precision Bethe-ansatz reference data for momentum-resolved MPS
excitation calculations, prioritizing short-range lattice Hamiltonians and
thermodynamic dispersions. Next sequence: **SU(3)/ULS excitations → spin-1
Takhtajan–Babujian excitations → sine-Gordon excitations**. ULS is the first
priority because it is a particularly demanding non-Abelian benchmark for MPS.

Goal statement: complete the required unchecked milestones below, following the
checkpoint completion rules. The earlier XXZ, XYZ and Potts sequence is complete;
the three new milestones below define the next goal.

Unchecked milestones define proposed work, not existing capabilities. Mark
milestones complete only when their documented public API/frontend and validation land.

## 1. SU(3)/ULS excitations

- [x] Extend the [SU(3) permutation-chain calculation](su3.md) from its balanced
  singlet ground state to zero-field thermodynamic elementary excitation branches,
  with physical momentum ranges and SU(3) representation labels.
- [x] Provide the relevant multiparticle continuum boundaries, auditing which
  combinations belong to physical periodic-chain sectors and which are accessible
  to local spin or quadrupolar operators. Do not equate an elementary branch with
  an isolated pole in a local response.
- [x] Document the conversion between permutation and spin-1 ULS Hamiltonians,
  including additive constants, energy scales and one-site versus three-site
  momentum folding for iMPS comparisons.
- [x] Validate dispersions, velocities and continuum edges against independent
  exact results; use selected finite-size BA/ED spectra as checks, with their
  finite-size corrections and representation content identified.

First deliverable: a public native-precision dispersion/continuum API and frontend,
not general finite-size excited-state enumeration. Reuse the existing nested-model
conventions and common continuum machinery where applicable. Starting references:
[Sutherland](../CITATIONS.md#sutherland-1975),
[Doikou–Nepomechie](../CITATIONS.md#doikou-nepomechie-1998), and
[SU(3) dynamical-spectrum comparisons](https://arxiv.org/abs/2107.09588).

Implemented: [ULS excitation guide](su3-dispersion.md), native-precision
`ExcitationDispersion` and `bethe-su3-dispersion`, with elementary 3/bar3 lines,
two-/four-soliton envelopes, three-site folding and CSV/TSV/JSON exports.
Validation: 1114 full-suite tests, 20 targeted checks in the fp128-enabled build,
and a pinned-Uni20 frontend build/export check. Independent momentum/Casimir ED
through L=12 verifies adjoint sectors and scaling, including a level below the
two-soliton onset; density identities and constrained scans check the analytic
curves. GitHub preserves the guide's 31 math expressions.

## 2. Spin-1 Takhtajan–Babujian excitations

- [x] Extend the [TB calculation](takhtajan-babujian.md) to zero-field
  thermodynamic spinon dispersions and the lowest relevant multiparticle
  continuum boundaries, with spin labels and physical momentum conventions.
- [x] Preserve the existing bilinear-minus-biquadratic Hamiltonian normalization;
  distinguish elementary fractional excitations from allowed periodic-chain
  multiplets and local-response sectors.
- [x] Validate exact dispersions and velocities independently, and compare
  continuum predictions with selected finite-size BA/ED levels. Account for
  finite-size effects and string deviations rather than treating ideal strings
  as exact finite-chain solutions.
- [x] Provide a public native-precision API/frontend and an MPS comparison guide,
  sharing kinematic and output helpers with ULS and the existing spinon tools.

General finite-size excited-string enumeration and dynamical spectral weights are
follow-ups, not prerequisites. Starting reference:
[Vlijm–Caux](../CITATIONS.md#vlijm-caux-2014).

Implemented: [TB excitation guide](tb-dispersion.md), native-precision
`SpinonDispersion` and `bethe-tb-dispersion`, with spin-1/2 lines, two-/four-spinon
envelopes, two-site folding and CSV/TSV/JSON exports. The guide distinguishes
local spin and quadrupole sectors, SU(2)_2 state counting, and the factor-four
Hamiltonian conversion from Vlijm–Caux. Shared helpers handle spinon kinematics,
streaming output completion and independent ED translation projection.
Validation: 1126 full-suite tests, 71 targeted checks in the fp128-enabled build,
and pinned-Uni20 TB/ULS frontend export checks. Independent spin/momentum-resolved
ED through L=10 checks normalization and finite-size scaling; native rapidity
identities and momentum scans check the analytic bounds. GitHub preserves the
guide's 16 math expressions.

## 3. Sine-Gordon excitations

- [ ] Extend the [sine-Gordon vacuum calculation](sine-gordon.md) with soliton,
  antisoliton and stable breather masses/dispersion branches, existence conditions,
  topological charges and multiparticle thresholds.
- [ ] Expose the scattering data needed for selected excited sectors and add
  large-volume Bethe–Yang levels as an explicitly asymptotic first checkpoint.
- [ ] Implement selected finite-volume excited levels using an appropriate
  excited-state NLIE/TBA, with an audited coupling/sector domain and state-selection
  rules. Do not present Bethe–Yang levels as exact finite-volume results.
- [ ] Validate against the free-fermion point, exact breather mass ratios,
  large-volume scattering quantization and appropriate ultraviolet conformal
  limits; retain separate quadrature, cutoff and nonlinear-solve diagnostics.
- [ ] Document the two-flavour Schwinger connection at equal small fermion masses:
  at theta=0 the leading light-sector theory has p=1/3, a triplet and a singlet
  with mass ratio sqrt(3). Label this as a scaling-limit benchmark, not an exact
  solution of the full massive or finite-spacing lattice Schwinger model.

Reuse the native-precision vacuum kernel and mass/coupling conventions. Publish
the thermodynamic/scattering and finite-volume stages as separate checkpoints;
both are required, but a complete finite-volume spectrum is not. Starting points:
[existing sine-Gordon references](sine-gordon.md#physical-convention) and the
[two-flavour Schwinger DMRG study](https://arxiv.org/abs/2407.11391).

## Completed sequence

The following milestones record the earlier completed goal. Their validation
counts are historical checkpoint results, not a claim about the current suite.

### Massive XXZ dispersions

- [x] Add `bethe-xxz-dispersion`, reusing the existing gapless zero-field
  formulas and extending to the zero-field antiferromagnet at Δ>1.
- [x] Provide bulk ground-state energy density, single-spinon ε(k), the
  single-spinon gap, and two-spinon lower/upper continuum edges.
- [x] Specify spin/sector labels, energy normalization, physical momentum and
  two-site-unit-cell folding. Distinguish a topological single spinon connecting
  different asymptotic vacua from a two-spinon excitation in a fixed vacuum sector.
- [x] Validate against independent exact expressions and isotropic/Ising limits,
  including precision-sensitive cases near Δ=1.

Implemented: [XXZ dispersion guide](xxz-dispersion.md), native-precision
`SpinonDispersion` and bulk-energy APIs, shared elliptic-band continuum
kinematics, and CSV/TSV/JSON frontend. Checkpoint validation: 1070 tests in the
full fp64/long-double suite and 21 targeted tests in the fp128-enabled build;
GitHub's Markdown renderer preserves the guide's 25 math expressions.

Finite-field/dressed dispersions are a follow-up, not required for this milestone.
Starting references: [Caux–Mossel–Pérez Castillo](https://arxiv.org/abs/0806.3069)
and the [symmetry-resolved MPS excitation framework](https://arxiv.org/abs/1802.07197).

### XYZ thermodynamic excitations

- [x] Audit the exact excitation formulas and select an explicitly documented
  parameter region; map its elliptic conventions to our existing XYZ Hamiltonian.
- [x] Add elementary thermodynamic dispersion branches, gaps and sector labels,
  reusing the native-precision elliptic-function infrastructure.
- [x] Add the bound-state branches present in the supported region, with
  explicit existence conditions and continuum thresholds.
- [x] Validate normalization and branches against XXZ/XY limits and independent
  references; use finite-size spectra as checks with finite-size effects identified.

Do not require a complete finite-ring XYZ spectrum or support for every coupling
region. Starting reference: [Johnson–Krinsky–McCoy](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.8.2526).

Implemented: [XYZ excitation guide](xyz-dispersion.md), native-precision spinon
and bound-branch API, and `bethe-xyz-dispersion` with parity/momentum labels and
CSV/TSV/JSON exports. Validation: 1084 full-suite tests, 21 targeted checks in
the fp128-enabled build, independent complex-rapidity references, XY and both
massive XXZ limits, and finite-ring checks with the splitting retained.

### Critical ferromagnetic three-state Potts chain

- [x] Fix the critical Hamiltonian, normalization and boundary conditions;
  implement the ground state and selected low-lying momentum/Z₃-resolved levels.
- [x] Audit physical-root selection and state counting before defining any
  excited-state enumeration or completeness claim.
- [x] Validate small chains against independent exact diagonalization and
  finite-size scaling against the known conformal spectrum.

This adds a discrete-symmetry, local-dimension-three benchmark. Generic off-critical
or chiral Potts models are outside this first scope. Starting reference:
[Dasmahapatra et al.](https://arxiv.org/abs/hep-th/9304150).

Implemented: [Potts guide and selection audit](potts.md), native-precision
vacuum/charged-one-hole API and `bethe-potts-pbc`, with physical momentum/Z₃
labels, finite-ring gaps, scaled dimensions and CSV/TSV/JSON exports. This is
2L selected charged levels, not a complete spectrum; neutral excitations remain
a follow-up. Independent charge/momentum-resolved clock diagonalization,
75-digit references, an excluded-XXZ-state counterexample and c=4/5 scaling
validate the selection. Closing regression: 1100 full-suite tests, 83 targeted
checks in the fp128-enabled build across the three checkpoints and shared XXZ
solver, and a pinned-Uni20/no-tests frontend build. Separate L=8,9 oracle audits
agree for every selected level. All three guides pass GitHub math rendering.

## Completion rules for every checkpoint

- Preserve fp64, native long-double and optional fp128 arithmetic throughout;
  include tests that detect precision narrowing and explicit numerical failures.
- Reuse Uni20 solves, CLI, tables and run metadata, and share numerical helpers
  between models where their mathematical contracts agree.
- Provide plotting-ready CSV/TSV/JSON, literature citations via `--references`,
  and concise documentation of units, sectors, supported ranges and limitations.
- Separate elementary branches, multiparticle thresholds and spectral weights;
  energies alone do not establish an observable's spectral intensity.
- Run relevant regressions, update the [model catalogue](models.md), and commit
  and push after each model or major user-facing checkpoint.

Wang-ladder excitations are deferred outside this goal. Spectral weights, general
finite-temperature dynamics, complete finite-size spectra and a general massive
Schwinger solver are not requirements of this goal.
