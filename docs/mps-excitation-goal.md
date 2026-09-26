# Goal: exact excitation benchmarks for MPS and iMPS

Develop native-precision Bethe-ansatz reference data for momentum-resolved MPS
excitation calculations, prioritizing short-range lattice Hamiltonians and
thermodynamic dispersions. Main sequence: **massive XXZ → XYZ excitations →
critical three-state Potts**, with an optional integrable-ladder checkpoint.

Goal statement: complete the required unchecked milestones below, following the
checkpoint completion rules. The optional ladder work does not block completion.

Unchecked milestones define proposed work, not existing capabilities. Mark
milestones complete only when their documented public API/frontend and validation land.

## 1. Massive XXZ dispersions

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

## 2. XYZ thermodynamic excitations

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

## Optional checkpoint: integrable-ladder excitations

- [ ] For Wang's ladder in its rung-singlet regime, expose the one-triplon
  dispersion, then selected two-triplon scattering/bound branches.
- [ ] Keep the required four-spin interaction and translation-by-one-rung
  momentum convention explicit; validate against small-system Hamiltonians.

This provides a conventional-particle benchmark alongside the topological
branches. It is not the generic Heisenberg ladder. Insert it when useful without
blocking the main sequence. Reference: [Wang](https://arxiv.org/abs/cond-mat/9901168).

## 3. Critical ferromagnetic three-state Potts chain

- [ ] Fix the critical Hamiltonian, normalization and boundary conditions;
  implement the ground state and selected low-lying momentum/Z₃-resolved levels.
- [ ] Audit physical-root selection and state counting before defining any
  excited-state enumeration or completeness claim.
- [ ] Validate small chains against independent exact diagonalization and
  finite-size scaling against the known conformal spectrum.

This adds a discrete-symmetry, local-dimension-three benchmark. Generic off-critical
or chiral Potts models are outside this first scope. Starting reference:
[Dasmahapatra et al.](https://arxiv.org/abs/hep-th/9304150).

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

SU(3)/ULS and spin-1 TB excitations remain later candidates. Spectral weights,
general finite-temperature dynamics and complete finite-size spectra are not
requirements of this goal.
