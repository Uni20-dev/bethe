# Bethe-ansatz model catalogue and development proposal

[Overview](../README.md) · [Bibliography](../CITATIONS.md)

Initial literature survey: 2026-09-22. Original implementation baseline: `6a414cc`;
the coverage below is updated as each implementation lands.
This is a living, **non-exhaustive catalogue**, with an emphasis on models
whose energies would be useful alongside Uni20 and MPToolkit. A model being
exactly solvable in the literature does not mean this repository solves it,
or that every boundary condition, coupling, or excited state is covered.

The first new implementations are **repulsive periodic Lieb–Liniger**,
the **periodic SU(3) balanced-singlet ground state**, and selected
**repulsive Gaudin–Yang** ground-state sectors. Next are wider XXZ ground-state
coverage and the supersymmetric t–J model; extending existing models'
boundary/state coverage remains valuable as well.
These priorities and difficulty assessments are our engineering judgments,
not conclusions of the cited papers or a committed implementation schedule.

## Reading and maintaining the catalogue

The entry IDs below are stable handles for issues and future model guides.
Status has a deliberately narrow meaning:

- **Implemented (limited):** executable/API coverage exists; follow its guide.
- **Proposed:** a useful first implementation is outlined, but no solver exists.
- **Watch:** established literature, with scope or numerical machinery still
  to investigate before implementation.
- **Related:** an exact benchmark worth considering, but not necessarily a
  conventional finite-size Bethe-root calculation.

For a proposed slice, **small** means mostly equations, conventions, and
validation using existing numerical tools; **medium** includes new sector or
nested-root handling; **large** includes complex roots, singular solutions,
functional equations, or a different spectral problem. These are relative
technical risks, not estimates in days. A small ground-state implementation
can have a much harder excitation extension.

When updating an entry, record the precise Hamiltonian and energy shift,
statistics/local representation, integrable coupling restrictions, boundary
conditions, supported sectors/root families, finite-size versus thermodynamic
scope, reference IDs, tests, and the next unresolved step. Once code exists,
link its guide and implementation commit here; do not turn an entire family
green because one slice works. New entries should have at least one original
solution or explicit Bethe-equation source. Before coding, do a separate
equation-level audit of the chosen finite-size branch and normalization.

Bibliographic metadata belongs in [data/citations.json](../data/citations.json).
Catalogue-only references appear in CITATIONS.md, **not in any executable's
help**, until a supported method actually uses them. See
[citation maintenance](citations.md). The links in this document identify
literature starting points; they are not a claim that we implement a cited
paper's correlation functions, thermodynamics, or full spectrum.

## Current coverage

| ID | Model | Status and actual scope | Main missing pieces |
| --- | --- | --- | --- |
| `xxx` | Spin-1/2 nearest-neighbor XXX | Implemented (limited): [PBC](xxx.md), [free ends](open-chains.md), sector minima and restricted real-root excitations; periodic one-spinon family | Complex strings, full spectrum, twists/boundary fields |
| `xxz` | Spin-1/2 nearest-neighbor XXZ | Implemented (limited): [PBC](xxz.md) and [free-end](xxz-open.md) ground states/sectors at `Delta>=0`; restricted excitations at `0<=Delta<=1`; [massive boundary roots](xxz-open-massive.md) | Negative anisotropy, massive excitations, additional root families, twists/boundary fields |
| `hubbard` | One-band Hubbard, hopping t=1 | Implemented (limited): [PBC](hubbard.md) on even rings with sector restrictions; [free ends](hubbard-open.md) at every physical filling/Sz and either sign of U | Hubbard excitations; remaining PBC shell branches and odd rings |
| `lieb-liniger` | Continuum contact-interacting bosons | Implemented (limited): [repulsive PBC](lieb-liniger.md), ground state, explicit labels, and bounded excitation scans | Hard walls, attraction, thermodynamics |
| `su-n` | Fundamental SU(n) permutation chain | Implemented (limited): [SU(3) PBC](su3.md), balanced singlet ground state for L>=3 divisible by three, J=1 | Other populations/lengths, excitations, general n, open boundaries |
| `gaudin-yang` | Equal-mass spin-1/2 continuum delta-interacting fermions | Implemented (limited): [repulsive PBC](gaudin-yang.md), odd populations of both spins; unrestricted free and fully polarized limits | Other periodic shell branches, excitations, attraction, hard walls, thermodynamics |

Analytic thermodynamic XXX/XXZ spinon dispersions are separate existing
facilities; they do not constitute a general thermodynamic Bethe ansatz
(TBA) solver. Similarly, `--excitations all` covers the documented finite
family, not the entire Hilbert space.

## Candidate index

Only entries marked **Implemented** are available; the other first scopes
are proposals, not CLI capabilities.
Boundary conditions shown are starting targets, not a classification of all
integrable boundaries in the literature.

| ID | Candidate and first scope | Status | Relative effort |
| --- | --- | --- | --- |
| `lieb-liniger` | Repulsive one-component Bose gas on a ring; ground state and bounded real-root excitation scans | [Implemented (limited)](lieb-liniger.md) | First slice complete |
| `su-n` | Fundamental SU(3) antiferromagnetic permutation chain, PBC, balanced ground state | [Implemented (limited)](su3.md) | First slice complete |
| `gaudin-yang` | Equal-mass repulsive spin-1/2 delta-interacting Fermi gas, PBC, selected ground-state sectors | [Implemented (limited)](gaudin-yang.md) | First slice complete |
| `tj-susy` | Projected t–J chain at J=2t, PBC, selected ground-state sectors | Proposed | Medium–large |
| `spin-s-tb` | Spin-1 Takhtajan–Babujian chain, PBC ground state with finite-size string deviations | Proposed | Large |
| `richardson` | Reduced BCS pairing Hamiltonian, specified levels and pair number | Proposed | Medium–large |
| `gaudin-magnet` | Rational Gaudin/central-spin spectra at specified couplings and magnetization | Watch | Medium–large |
| `multicomponent-gas` | SU(kappa) fermions or the equal-coupling Bose–Fermi mixture, PBC | Watch | Medium–large |
| `integrable-ladder` | A specified SU(4)-type ladder with the required four-spin interaction, PBC | Watch | Medium after SU(n) |
| `q-boson` | Integrable q-boson hopping/phase model, PBC at fixed particle number | Watch | Medium |
| `xyz` | Zero-field spin-1/2 XYZ chain, PBC finite-size spectrum | Watch | Large |
| `haldane-shastry` | Inverse-chord-square spin chain, exact finite-ring levels | Related | Small for energies; larger for state counting |
| `sutherland` | Trigonometric inverse-square gas on a circle, exact pseudomomentum energies | Related | Small–medium |
| `kondo` | The integrable continuum single-impurity Kondo problem | Watch | Large; different physical scope |
| `sine-gordon` | Integrable quantum field theory, finite-volume ground energy via nonlinear integral equations | Watch | Large |
| `asep` | Periodic asymmetric exclusion process, relaxation spectrum | Watch | Large; non-Hermitian, not an energy spectrum |

## Leading proposals

### `lieb-liniger`: the simplest new interacting family

The continuum Bose gas has contact interactions. In units `hbar^2/(2m)=1`,
one common convention is `H=-sum_j d_j^2 + 2c sum_(i<j) delta(x_i-x_j)`.
For c>0 its finite-ring Bethe momenta are real, with `E=sum_j k_j^2`.
The original solution and the two excitation branches are in
[Lieb–Liniger I](../CITATIONS.md#lieb-liniger-1963) and
[Lieb II](../CITATIONS.md#lieb-1963-excitations).

Implemented in [lieb_liniger.hpp](../include/bethe/lieb_liniger.hpp) and
`bethe-lieb-liniger-pbc`: physical length ell, N bosons, c>0, consecutive ground
labels, explicit labels, and scans over N occupied slots in an N+2P window.
Implementation commit:
[`4ba617c`](https://github.com/Uni20-dev/bethe/commit/4ba617c9eaebaea9f25af5fe1f8e2d2ec14e1baf).
See the [guide](lieb-liniger.md) for units, parity, residual scaling and limits.
Tests cover two-/three-body analytic values, weak coupling, the
impenetrable-boson (Tonks–Girardeau) limit, boosts, all three precisions, finite
windows, and convergence to a separately discretized bulk integral equation.
The implementation checkpoint passed 285 tests with GCC 13 and fp128 enabled,
198 with Clang 20 Release without MPLAPACK, and the new CLI suite in an
app-only build using the published Uni20 pin (no sibling-source override).

Hard walls are a natural second slice with reflected scattering, but require
their own equations: [Gaudin (1971)](../CITATIONS.md#gaudin-1971).
Attractive c is a separate bound-state problem, not a sign toggle on a
real-root solver. Nor does `all` make sense without an energy or quantum-number
cutoff: even a fixed-N continuum system has infinitely many levels.

**Next slice:** hard-wall reflection equations, or a thermodynamic ground-state
and type-I/type-II dispersion API. The finite-ring first slice is complete;
attraction remains a separate bound-state project.

### `su-n`: permutation chains and the spin-1 ULS point

Take one fundamental n-state degree of freedom per site and the
antiferromagnetic Hamiltonian `H=J sum_j P_(j,j+1)`, J>0, where P exchanges
neighboring colors. This is the SU(n) permutation chain associated with
[Sutherland's multicomponent solution](../CITATIONS.md#sutherland-1975).
It is a direct lattice benchmark for non-Abelian tensor-network calculations.

SU(3) is especially attractive: in a spin-1 basis,
`P = S_i.S_j + (S_i.S_j)^2 - 1`. Thus it also supplies the
Uimin–Lai–Sutherland (ULS) bilinear–biquadratic point, after an explicit
per-bond energy shift. The identity follows from the two-site total-spin
eigenvalues; this is not the generic spin-1 Heisenberg chain.

Implemented in [su3.hpp](../include/bethe/su3.hpp) and `bethe-su3-pbc`:
periodic, J=1, balanced singlet for L>=3 divisible by three, with finite real
roots. The state keeps explicit populations and two separate root/label
arrays, with M1=2L/3 and M2=L/3. Implementation commit:
[`3970675`](https://github.com/Uni20-dev/bethe/commit/39706753e542123802714dee014ed1b828d98a8d).
See the [guide](su3.md) for the equations,
normalization, ULS conversion, and numerical limits. The equation audit uses
[Doikou–Nepomechie](../CITATIONS.md#doikou-nepomechie-1998), including the
filled-sea selection in Sec. 2.3; our energy is twice the paper's plus L.

Tests cover three-/six-site analytic roots or energies in every precision,
independent color-space ED and translation through L=9, original multiplicative
and logarithmic equations, Jacobians, the SU(2)/XXX reduction, and the spin-1
permutation identity. Finite chains through L=192 approach the known bulk
energy density. The implementation checkpoint passed 315 tests with GCC 13
and fp128, 219 with Clang 20 Release without MPLAPACK, and the new frontend
and citation checks against the published Uni20 pin in an app-only build.

**Next slice:** audit other color sectors and hole/excitation families,
including lengths not divisible by three. General n adds n-1 nesting levels;
complex strings, descendants, twists, and open ends need their own state and
equation treatment. The two-level SU(3) state does not change the Hubbard API
or claim an arbitrary-rank solver.

### `gaudin-yang`: the continuum counterpart of nested Hubbard

This is the equal-mass spin-1/2 Fermi gas with delta interactions, solved by
[Gaudin](../CITATIONS.md#gaudin-1967) and
[Yang](../CITATIONS.md#yang-1967); Yang's
[scattering construction](../CITATIONS.md#yang-1968) also treats attraction.
It combines continuum charge momenta with an auxiliary spin-rapidity family.

The first slice should be repulsive periodic ground states in explicitly
validated particle/spin sectors, followed by polarization and charge/spin
excitations. Reuse nested Newton and continuation concepts, not Hubbard's
`sin(k)` kernels, energy formula, or finite-ring shell selection. Validate
the noninteracting and fully polarized gases, the two-body problem, and
strong-coupling behavior; the dilute Hubbard limit is an additional check
only after matching units and taking a controlled continuum limit.

Attraction introduces paired complex charge roots. The bipartite-lattice
Shiba shortcut used for Hubbard is not a corresponding continuum sign-change
mapping; finite-size paired states need their own treatment. Multicomponent
fermions add more nesting levels and, for attraction, larger bound complexes;
see [Lee–Guan–Batchelor](../CITATIONS.md#lee-2011).

**Implemented first slice:** [commit 6b7fc3c](https://github.com/Uni20-dev/bethe/commit/6b7fc3c)
adds `bethe-gaudin-yang-pbc` and the scalar-templated
`bethe::gaudin_yang::ground_state(N_up,N_down,ell,c)` API. Interacting
mixed-spin states require odd populations of both spins; c=0 and full
polarization use exact free occupations for any particle count. The
[model guide](gaudin-yang.md) specifies physical units, nested labels,
free-current degeneracy, spin reversal, and incomplete-continuation diagnostics.
The explicit equations and sector choice follow
[Oelkers et al.](../CITATIONS.md#oelkers-2006), especially Eq. (25) and Sec. 5.

Validation includes the original rational equations, an independent
two-body jump-condition oracle, exact free energies, weak/strong-coupling
coefficients, the XXX spin-root limit, controlled dilute-Hubbard convergence,
analytic-Jacobian checks, and fp64/long-double/fp128 precision tests.
The complete suite at this checkpoint passes 344 GCC Debug tests with fp128
and 239 Clang Release tests without MPLAPACK; the front end also passes its
CLI checks against the published Uni20 pin without a sibling-source override.

**Next slice:** other periodic shell branches need a separate finite-ring
state-selection audit, not just shifted centered labels. Then consider
charge/spin excitations or hard-wall reflection equations. Attraction still
requires paired complex roots; no lattice particle-hole shortcut applies.

### `tj-susy`: a useful strongly correlated lattice benchmark

The local states are empty, up, and down: double occupancy is projected out.
For the convention
`H=-t sum(projected hopping+h.c.) + J sum(S_i.S_(i+1)-n_i*n_(i+1)/4)`,
target the supersymmetric **J=2t** point, not arbitrary J/t. Use the explicit
alternative nested Bethe ansätze in
[Essler–Korepin](../CITATIONS.md#essler-korepin-1992).

Start with periodic ground states in a specified filling/magnetization
family. The choice of grading/reference state changes the root description;
we should choose a numerically suitable formulation and derive its particle
counts and energy offsets explicitly. Hubbard's quantum-number sea cannot
simply be copied. Validate the projected Fock-space Hamiltonian, the
fully polarized free-fermion sector, and the no-hole reduction to XXX with
the `-J/4` bond shift. Integrable open boundaries exist, but should follow
as a separate boundary-equation project, not be inferred from the PBC code.

**Next decision:** compare the available gradings on a few small systems
before choosing the first supported real/complex root family.

### `spin-s-tb`: spin-1 first, with genuine complex-root support

The integrable higher-spin chains of
[Babujian](../CITATIONS.md#babujian-1982) have specially chosen polynomial
exchange interactions. At spin 1 the Takhtajan–Babujian (TB) point is
proportional to `sum[S_i.S_(i+1) - (S_i.S_(i+1))^2]`.
It is neither the ULS point above nor the generic bilinear spin-1 chain.

Although there is only one nesting level, the antiferromagnetic ground
state already involves complex two-string patterns. Finite-size deviations
matter: [Vlijm–Caux](../CITATIONS.md#vlijm-caux-2014) provides a numerical
starting point. Solving ideal string centers alone must be labelled an
approximation, not an exact finite-chain energy calculation.

Start with even periodic spin-1 chains, small enough for independent exact
diagonalization, and check the original complex equations after reconstructing
the roots. Broader spin, excitations, and open boundaries follow only after
singular/colliding roots and missing-state diagnostics are under control.

**Next dependency:** a validated complex-root or regularized string-deviation
solver. This machinery also benefits missing XXX/XXZ and Hubbard excitations.

### `richardson`: finite pairing spectra without a spatial chain

The reduced BCS Hamiltonian has specified single-particle levels and a
uniform pair-scattering coupling. Richardson's solution and its Gaudin
generalizations are surveyed by the model developers in
[Dukelsky–Pittel–Sierra](../CITATIONS.md#dukelsky-2004).
This would broaden the library to pairing benchmarks at fixed particle number.

Start with distinct levels, a fixed pair number, and the lowest state in a
specified blocked-level sector (singly occupied levels do not participate in
pair scattering). Pair rapidities can become complex and collide with poles
as coupling varies. Consider regularized or eigenvalue-based equations as
well as direct root continuation; a naive complex Newton replacement is not
enough. Check zero coupling, one pair, and exact diagonalization of small
pairing spaces across root collisions.

**Next decision:** choose the regularized variables and state-tracking rule.
Inputs are levels, degeneracies, blocked occupations, and coupling; PBC/OBC
and lattice momentum are not appropriate interface concepts here. Keep a
separate front end rather than forcing this into a chain-shaped CLI.

## Extensions of current models, rather than new solvers

- **Wider XXZ anisotropy:** the [periodic ground-state/sector path](xxz.md#easy-axis-ground-states-delta1)
  now includes `Delta>1`, retaining a continuous atan2 scattering phase in
  scaled real coordinates. The source is [Dugave et al.](../CITATIONS.md#dugave-2015),
  Eqs. (1.1)-(1.2), with energy normalized from Pauli matrices to spin operators.
  Tests cover odd/even sector energies and momentum against ED, original
  trigonometric equations, native precision, continuity to XXX, and Ising limits.
  The [massive free-end module](xxz-open-massive.md) now solves the
  ground-state sectors through the same ground-state API and `bethe-xxz-obc`,
  retaining the finite-size boundary-root deviation explicitly.
  An [internal negative-anisotropy engine](xxz-negative.md) now covers
  `-1<Delta<0` sector ground roots for even rings and open chains, with
  rank-subtracted equations scaled to remain discriminating near Delta=-1.
  It is not yet exposed by the public ground-state API or frontends.
  The [polynomial equations](xxz-polynomial.md) and an internal
  [momentum-constrained adaptive driver](xxz-odd-continuation.md) now handle
  conjugate-pair formation without extracting roots, with independent
  spin-basis energy checks through 21 sites at selected couplings.
  **Next:** physical-state admissibility and global sector selection,
  then public integration; only subsequently extend excitation classification.
  The broader family is established by [Yang–Yang](../CITATIONS.md#yang-yang-1966), but the
  old real-coordinate domain and excitation window do not extend unchanged.
  Treat `Delta<=-1` and polarized/degenerate limits separately. These remain
  XXZ-module extensions, not differently named models.
- **Spinless t–V fermions / nearest-neighbor hard-core bosons:** useful physical
  interfaces to XXZ, not independent Bethe engines. Under the usual convention
  the anisotropy is V/(2t), but the hopping-sign gauge, fermionic parity twist
  on a ring, density-dependent constants, and OBC endpoint terms must be
  derived for the declared Hamiltonian. Never reuse a spin-chain energy or
  momentum without this translation.
- **Uniform fields and chemical potential:** a commuting term gives
  `E(h,mu)=E(0,0)-h*Sz-mu*N` within a solved sector. The global minimum still
  requires comparison across all relevant sectors; incomplete PBC Hubbard
  coverage cannot give an unrestricted grand-canonical minimum. This needs
  sector scanning/reporting, not new Bethe equations.
- **Flux/twists:** useful for persistent currents and stiffness; the boundary
  phases enter the equations and can change the lowest branch. See
  [Shastry–Sutherland](../CITATIONS.md#shastry-sutherland-1990).
  Track crossings and symmetry mappings at the requested twist instead of
  retaining the zero-flux quantum numbers unconditionally.
- **Integrable boundary fields:** start with diagonal boundaries preserving
  the conserved projection; consult the reflection-matrix construction of
  [de Vega–González-Ruiz](../CITATIONS.md#de-vega-gonzalez-ruiz-1994) and the
  existing [open Hubbard reference](../CITATIONS.md#deguchi-yue-1997).
  Boundary bound states can change the root content. General boundary fields
  are not obtained by merely adding a post-solve energy term.
- **Hubbard excitations and missing ring sectors:** start with explicitly
  admissible nested label changes and small-ring ED; add complex spin and
  charge–spin bound states only with appropriate equations. Keep energies
  produced by symmetry mappings distinct from physical attractive roots.
  The [Hubbard guide](hubbard.md) and its literature are the starting point.

## Wider catalogue: useful, but not the next default targets

The following entries deliberately remain short until selected for an
equation-level feasibility study. Each has a concrete restriction and an
initial deliverable, rather than an unqualified claim of model support.

| ID | Literature and integrable restriction | Useful first deliverable / main obstacle |
| --- | --- | --- |
| `gaudin-magnet` | Rational Gaudin magnets; a concrete central-spin realization is treated by [Faribault–Schuricht](../CITATIONS.md#faribault-schuricht-2013) | Small-system energies in fixed magnetization sectors; share regularization ideas with Richardson, but specify the exact integrable coupling/field family, not arbitrary spin-bath interactions |
| `multicomponent-gas` | Equal-mass SU(kappa) delta fermions, [Lee et al.](../CITATIONS.md#lee-2011); equal-mass Bose–Fermi mixture with equal repulsive Bose–Bose/Bose–Fermi couplings, [Imambekov–Demler](../CITATIONS.md#imambekov-demler-2006) | PBC repulsive ground energies at fixed component counts; more nesting/statistics bookkeeping. A trapped local-density calculation would be an approximation, not an exact trapped BA solution |
| `integrable-ladder` | [Wang's ladder](../CITATIONS.md#wang-1999), with its required exchange and four-spin terms | Ground energies and rung-sector competition; an application of higher-rank nesting. Not the ordinary two-leg Heisenberg ladder at generic couplings |
| `q-boson` | Deformed boson hopping and its phase-model limit, [Bogoliubov–Izergin–Kitanine](../CITATIONS.md#bogoliubov-1997) | PBC finite-N energies, checked in a bosonic occupation basis; define the deformed local algebra and coupling limits. Not the standard Bose–Hubbard chain |
| `xyz` | Zero-field spin-1/2 XYZ / eight-vertex family, [Baxter](../CITATIONS.md#baxter-1973) | Small periodic-chain energies via elliptic/functional equations; Sz is not generally conserved. Native-precision elliptic functions and branch selection are substantial new requirements |
| `haldane-shastry` | Spin-1/2 inverse-chord-square exchange on a ring, [Haldane](../CITATIONS.md#haldane-1988) and [Shastry](../CITATIONS.md#shastry-1988) | Exact energy and spinon-state enumeration benchmarks; implement the known spectral rules, not artificial Newton roots. Multiplet/motif counting is distinct from the nearest-neighbor XXX problem |
| `sutherland` | The trigonometric inverse-square gas, [Sutherland](../CITATIONS.md#sutherland-1971), in the Calogero–Sutherland family | Exact pseudomomentum energies on a circle; specify statistics, coupling branch, and normalization. Keep algebraic/asymptotic-BA descriptions distinct from generic finite-range Bethe equations |
| `kondo` | Integrable continuum single-impurity Kondo Hamiltonian, [Andrei](../CITATIONS.md#andrei-1980) | Impurity ground-energy/thermodynamic benchmarks after fixing band regularization and bulk subtraction; not a generic finite-band Anderson impurity or Kondo lattice |
| `sine-gordon` | Continuum sine-Gordon theory, e.g. [Destri–de Vega](../CITATIONS.md#destri-de-vega-1992) | Finite-volume ground-state scaling function; requires integral/functional equations, mass/coupling conventions, and vacuum-energy subtraction. Large-volume Bethe–Yang quantization alone omits finite-volume corrections |
| `asep` | Periodic asymmetric simple exclusion process, [Gwa–Spohn](../CITATIONS.md#gwa-spohn-1992) | Relaxation gap at fixed particle number; eigenvalues can be complex and are decay rates, not Hermitian energies. Requires a separate spectral-state/reporting contract |

Other recognized families can be added as research needs arise: anisotropic
multicomponent/Perk–Schultz chains, integrable higher-spin anisotropic chains,
supersymmetric extended Hubbard models, and quantum Toda systems. They are
**not yet individual proposals here**: choose an explicit Hamiltonian,
representation, boundary, and source before creating a tracked entry.

## Thermodynamics is a separate capability

There are two different questions: "what is the energy of this finite-chain
state?" and "what is the bulk energy density or finite-temperature free
energy?" Root-density equations and TBA answer the second question without
enumerating finite-system eigenstates. The natural first target is the
repulsive Lieb–Liniger Yang–Yang equation,
[Yang–Yang (1969)](../CITATIONS.md#yang-yang-1969); it has a much simpler
hierarchy than spinful/nested thermodynamics.

Propose a separate `thermodynamics` layer with native-precision quadrature,
cutoff/grid convergence, and a clearly specified canonical or grand-canonical
ensemble. Check the zero-temperature equation of state and free/strong-coupling
limits. Compare thermodynamic and large finite-N results, but never substitute
one for the other silently. Quadrature and domain-cutoff errors must be tracked
separately from nonlinear residuals. This is a better first integral-equation
project than attempting a universal TBA engine for every catalogue entry.

## Shared machinery and module boundaries

Continue with **one model module and thin, boundary-appropriate front ends**.
Share numerical mechanisms when their contracts really agree:

| Reusable facility | Models that motivate it | Boundary of responsibility |
| --- | --- | --- |
| Real-root Newton, Jacobian checks, continuation | Hubbard, Lieb–Liniger, Gaudin–Yang, SU(3) | Each model owns its residuals, allowed branches, seeds, energy, and labels |
| Nested-family state bookkeeping | SU(n), Gaudin–Yang, t–J, mixtures | Do not assume exactly two families or identify all auxiliary roots with physical momenta |
| Complex-root / regularized equations | TB, attractive gases, spin-chain/Hubbard strings, Richardson | Preserve conjugacy, handle singularities, and verify original equations; ideal strings are not automatically exact |
| Quantum-number enumeration and state tracking | All finite-size excitation solvers | A scan window, multiplicity policy, and failed/missing-state counts are model-specific |
| Quadrature and nonlinear integral equations | Bulk LL, nested TBA, sine-Gordon | Separate numerical truncation errors and thermodynamic outputs from finite-root results |

Keep scalar precision, presentation, timing, and citation rendering shared.
Continue using Uni20 for numerical linear algebra; if complex solves or
special functions expose a genuine missing facility, specify its precision,
failure semantics, and tests in cooperation with Uni20 development. Do not
introduce a silent double-precision dependency into a high-precision model.
No language change or giant runtime-switched "all models" solver is needed.

## Suggested sequence and first acceptance evidence

1. **Repulsive Lieb–Liniger PBC:** ground state and explicit labels, followed
   by finite-window excitations; two-body and limiting-case checks in every
   precision. Hard walls are the next boundary slice, not part of the first
   implementation by implication.
2. **SU(3) PBC balanced singlet:** first slice complete, with small-chain ED,
   native-precision analytic checks, and documented nested labels.
   **Gaudin–Yang PBC:** the first odd-population sector slice is also complete,
   with independent two-body/limiting-case checks, dilute-Hubbard convergence,
   and documented physical and auxiliary labels.
3. **Wider XXZ ground-state coverage** and **t–J at J=2t**: extend lattice
   benchmarks while keeping branch/parameter restrictions explicit.
4. **A focused complex-root project:** start with a small known XXX string
   and the spin-1 TB ground state, or Richardson pairing if that application
   matters more. Success means finite-size deviations and collisions are
   tested, not just that an ideal-string energy looks plausible.
5. **Lieb–Liniger thermodynamics**, then whichever wider entry a concrete
   Uni20/MPToolkit calculation needs. The catalogue can grow independently
   of this suggested order.

For each implementation, require an independent Hamiltonian or limiting-case
check, original-equation residual checks, native long-double/fp128 tests that
would detect narrowing, and explicit failure reporting. Ground-state claims
need state-selection evidence, not just a converged root set. Completeness
claims additionally need counting, degeneracy, and singular-state tests.

## What this catalogue must not imply

Generic spin-1 bilinear exchange, arbitrary bilinear–biquadratic angle,
generic Bose–Hubbard chains, arbitrary t–J ratios, and arbitrary ladder or
boundary couplings are **not supplied with a Bethe solver merely by their
relation to an integrable example**. An exactly known ground state (such as
an MPS parent-Hamiltonian state) is also not an exact Bethe spectrum. The
tracked objects are explicit integrable Hamiltonians and well-defined
solution sectors, not broad model names stripped of their restrictions.
