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
**repulsive Gaudin–Yang** ground-state sectors. Wider XXZ ground-state
coverage is now public for massive anisotropy and for negative anisotropy
on even rings and free-end chains. Negative-anisotropy **odd periodic rings
are deferred until needed**, not a prerequisite for other models.
The supersymmetric t–J chain now has a first periodic ground-state slice.
The spin-1 Takhtajan–Babujian chain also has an even-ring singlet solver
retaining finite-size complex-root deviations.
The pure spin-1 biquadratic chain now has even, free-end ground states,
TL module minima and restricted real-root excitations with physical multiplicities,
plus selected complex-root levels, budgeted Q-system searches and a targeted
long-chain two-string singlet branch.
The ferromagnetic sign adds an exact odd/even one-defect band and positive
gap, targeted odd/even long-chain two-/three-defect bound droplets and real-root
scattering windows, plus sign-aware even-chain Q-system discoveries.
Richardson pairing now supplies attractive ground energies in specified
blocked-level sectors, using regular variables through pair-root collisions.
The rational central-spin model now has fixed-magnetization ground states
for distinct nonzero couplings, including exactly zero central field.
The multicomponent continuum entry also has a repulsive SU(κ) fermion
ground-state slice for odd occupied populations, with arbitrary nesting depth.
Extending existing models' boundary/state
coverage remains valuable as well.
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
| `xxz` | Spin-1/2 nearest-neighbor XXZ | Implemented (limited): [PBC](xxz.md) ground states/sectors at `Delta>=0`, also `-1<Delta<0` on even rings; [free ends](xxz-open.md) at `Delta>-1`; restricted excitations at `0<=Delta<=1`; [massive boundary roots](xxz-open-massive.md) | Negative odd rings, massive/negative excitations, additional root families, twists/boundary fields |
| `hubbard` | One-band Hubbard, hopping t=1 | Implemented (limited): [PBC](hubbard.md) and [free-end](hubbard-open.md) ground states; [half-filled](hubbard-dispersion.md) and [doped](hubbard-doped.md) thermodynamic spinon/charge lines for U>0, zero field | Continuum thresholds, finite-field/string dispersions; finite-size excitations; remaining PBC shell branches and odd rings |
| `lieb-liniger` | Continuum contact-interacting bosons | Implemented (limited): [repulsive PBC](lieb-liniger.md) and [hard walls](lieb-liniger-open.md), ground states, explicit labels, bounded excitation scans; [bulk ground state and type-I/type-II curves](lieb-liniger-thermo.md); [grand-canonical and fixed-density finite-T equilibrium and temperature scans](lieb-liniger-thermal.md) | Attraction, form factors |
| `q-boson` | Deformed boson hopping on a lattice | Implemented (limited): [fixed-N PBC ground states and excitation scans](q-boson.md), eta>=0, free and phase limits | Open boundaries, thermodynamics |
| `su-n` | Fundamental SU(n) permutation chain | Implemented (limited): [SU(3) PBC](su3.md), balanced singlet ground state for L>=3 divisible by three, J=1 | Other populations/lengths, excitations, general n, open boundaries |
| `gaudin-yang` | Equal-mass spin-1/2 continuum delta-interacting fermions | Implemented (limited): [repulsive PBC](gaudin-yang.md), odd populations of both spins; unrestricted free and fully polarized limits | Other periodic shell branches, excitations, attraction, hard walls, thermodynamics |
| `tj-susy` | Projected t–J electrons, J=2t | Implemented (limited): [PBC](tj.md), t=1, doped odd N_up and N_down on either length parity; every no-hole and fully polarized sector | Other doped shell branches, excitations, open boundaries |
| `spin-s-tb` | Integrable spin-1 bilinear–biquadratic chain | Implemented (limited): [even PBC](takhtajan-babujian.md), singlet ground state of H=sum[S.S-(S.S)^2], with finite two-string deviations | Odd lengths, sectors, excitations, higher spins, open boundaries |
| `temperley-lieb` | TL singlet-projector chains; spin-1 pure biquadratic model | Implemented (limited): [even free ends](biquadratic.md), AF ground state, module minima, real-root excitations, Q-system searches and targeted two-string singlet; [ferro one-defect band](biquadratic-ferromagnetic.md), [bound pairs](biquadratic-bound-pairs.md), [three-defect droplets](biquadratic-bound-triples.md), [real-root scattering windows](biquadratic-scattering.md) and [mixed pair-plus-defect scans](biquadratic-pair-defect.md) (odd/even), sign-aware Q-system levels; representation multiplicities, [physical-spin content API and CLI](biquadratic-spin-content.md), and generic lambda>2 TL API | More general mixed strings and larger droplets, other complex-root families, singlet ranks, odd-chain AF spinon branch, PBC twists, other representations |
| `richardson` | Reduced BCS pairing | Implemented (limited): [attractive pairing](richardson.md), distinct doublet levels, fixed pair count and blocked levels, ground energy through root collisions | Repeated levels/higher degeneracies, excitations, pair-root output, repulsive coupling |
| `gaudin-magnet` | Rational spin-1/2 central spin | Implemented (limited): [sector minima](central-spin.md), distinct nonzero bath couplings of either sign, central field of either sign or zero | Repeated/zero couplings, higher local spins, excitations, general Gaudin charges |
| `multicomponent-gas` | Equal-mass SU(κ) delta fermions | Implemented (limited): [repulsive PBC](su-fermions.md), odd occupied populations, any number of components, unrestricted free/single-component limits | Other periodic shells, attraction, excitations, hard walls, TBA |
| `bose-fermi` | Equal-mass, equal-repulsion scalar Bose–Fermi gas | [Ground-state library/frontend](bose-fermi.md): PBC, odd fermion population; unrestricted pure/free limits | Other shells, excitations, thermodynamics |
| `integrable-ladder` | Wang's spin-1/2 ladder with four-spin exchange | Implemented (limited): [zero-field PBC](ladder.md), global and singlet-count sector ground energies, either sign of J_r, leg coefficient 1 and four-spin coefficient 4 | Fields, excitations, open ends, other integrable ladder families |

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
| `tj-susy` | Projected t–J chain at J=2t, PBC, selected ground-state sectors | [Implemented (limited)](tj.md) | First slice complete |
| `spin-s-tb` | Spin-1 Takhtajan–Babujian chain, PBC ground state with finite-size string deviations | [Implemented (limited)](takhtajan-babujian.md) | First slice complete |
| `temperley-lieb` | Spin-1 pure biquadratic even free-end ground and real-root excited levels, with TL multiplicities | [Implemented (limited)](biquadratic.md) | Ground/excitation slice complete |
| `richardson` | Reduced BCS pairing Hamiltonian, specified levels and pair number | [Implemented (limited)](richardson.md) | First slice complete |
| `gaudin-magnet` | Rational spin-1/2 central-spin sector ground energies at specified couplings and field | [Implemented (limited)](central-spin.md) | First slice complete |
| `multicomponent-gas` | Repulsive SU(kappa) fermions on a ring, fixed odd occupied populations | [Implemented (limited)](su-fermions.md) | Fermion first slice complete |
| `bose-fermi` | Equal-mass scalar Bose–Fermi mixture with equal repulsive BB/BF couplings, PBC | [Ground-state library/frontend implemented](bose-fermi.md), odd fermion shells | First slice complete |
| `integrable-ladder` | Wang's SU(4)-type ladder with its required four-spin interaction, zero-field PBC sector minima | [Implemented (limited)](ladder.md) | First slice complete |
| `q-boson` | Integrable q-boson hopping/phase model, PBC at fixed particle number | [Ground-state and excitation library/frontend implemented](q-boson.md) | Canonical finite-size scans complete |
| `xyz` | Zero-field spin-1/2 XYZ chain, PBC finite-size spectrum | [Even-chain regular ground branch and frontend implemented](xyz.md) | Excited spectrum remains large |
| `haldane-shastry` | Inverse-chord-square spin-1/2 ring, ground/sector minima and bounded complete motif spectra | [Implemented](haldane-shastry.md) | First slice complete |
| `sutherland` | Periodic bosonic inverse-square gas, exact pseudomomentum energies and bounded label scans | [Implemented](sutherland.md) | First slice complete |
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

The [hard-wall library and CLI](lieb-liniger-open.md) add ground states, selected
positive-integer labels and finite-window scans using reflected scattering.
It shares Newton iteration and rational kernels with the ring, but owns its
boundary equations: [Gaudin (1971)](../CITATIONS.md#gaudin-1971).
Attractive c is a separate bound-state problem, not a sign toggle on a
real-root solver. Nor does `all` make sense without an energy or quantum-number
cutoff: even a fixed-N continuum system has infinitely many levels.

The [thermodynamic library and frontend](lieb-liniger-thermo.md) now supply
Q, energy per length, chemical potential, and type-I/type-II dispersions using
mesh-verified integral equations and physical-momentum inversion.
The repulsive finite-size and zero-temperature bulk slices are complete;
The [finite-temperature Yang–Yang library](lieb-liniger-thermal.md) also supplies
grand-canonical pressure, density, energy, and entropy with separate mesh and
cutoff checks, plus fixed-density inversion and a temperature-scan frontend;
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

Implemented in [tj.hpp](../include/bethe/tj.hpp) and `bethe-tj-pbc`, using
Sutherland's BFF grading: `M1=N_h+min(N_up,N_down)`, `M2=N_h`. The doped
real-root family requires odd populations of both spins, on odd or even L.
No-hole states use `H_tJ=2*H_XXX-L/2`; fully polarized states are exact free
fermions for every particle count. The [guide](tj.md) derives the label
parities, energy shift, and fermionic translation phase and states the
unsupported sectors explicitly.

Tests compare energies and momenta with an independent projected Fock-space
Hamiltonian through L=8, check the original rational equations and analytic
Jacobian, and retain native precision in analytic limits and solves through
L=96. Invalid shell branches are rejected, not replaced by a nearby filling
or spin. There is no mu-mu self-scattering term; SU(3)'s nested equations
cannot be reused unchanged.

**Next slice:** audit other doped population parities and their real/complex
root families before extending the public state selector. Excitations and
integrable open boundaries need separate labels/equations. Arbitrary J/t
is outside this model's implemented integrable point.

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

Implemented in [takhtajan_babujian.hpp](../include/bethe/takhtajan_babujian.hpp)
and `bethe-tb-pbc`: even L>=4, zero-field singlet, bilinear coefficient 1
(J=4 in Vlijm–Caux). The filled sea's centres and positive deviations are
solved together with an analytic Jacobian; ideal strings are not substituted
for the finite roots. [The guide](takhtajan-babujian.md) records the branch,
equations, normalization, precision and failure contract.

Tests check original complex equations through L=128 in all scalar types,
independent spin-basis energies and translation at L=4,6,8, native-precision
E_4=-11-sqrt(41), the Jacobian and incomplete-solve diagnostics.

**Next slice:** broken-string excitations require real and three-string
roots, singular-solution handling and new label branches. Other spin,
magnetization, odd lengths and open boundaries are not covered by this
ground-state implementation.

### `temperley-lieb`: shared energies, different representation multiplicities

The spin-1 pure biquadratic chain `H=-sum(S.S)^2` realizes the open TL
algebra at loop weight 3, with `e_i=(S_i.S_(i+1))^2-1`. The corresponding
spin-1/2 XXZ representation has Delta=3/2 **and opposite end fields**:
ordinary zero-field `bethe-xxz-obc` is not the reference Hamiltonian.
The physical energy is `2*E_ref-7*(N-1)/4` in our spin-half exchange-1
normalization. Within a TL module the energies agree after this mapping,
but representation multiplicities and physical spin labels do not.

Implemented in [temperley_lieb.hpp](../include/bethe/temperley_lieb.hpp),
[biquadratic.hpp](../include/bethe/biquadratic.hpp), and `bethe-biquadratic-obc`:
even N>=2, free-end ground state, TL module minima and restricted real-root
excitations. `--through-lines` selects a module, `--sectors` lists minima,
and `--excitations COUNT|all` enumerates the supported label family (default
ell=2). Each level carries the spin-1 representation multiplicity, while
gaps refer to the global singlet ground state. Here `all` never denotes the
full module spectrum. A separate [Q-system solver](xxz-open-qsystem.md) adds
selected complex-root levels (`--q-seed`) and bounded small-chain module
searches (`--q-spectrum`, validated through N=8; larger N experimental); missing/unverified levels produce an
explicitly incomplete result, not a lowest-level guarantee. The lower TL layer supports
the same real-root slice at general lambda>2 and a checked integer multiplicity
helper for singlet-projector spin-chain representations.
See the [guide](biquadratic.md) for the Bethe equations, energy shifts,
precision contract, and full small-chain ED spectral checks. The separate
[two-string solver](xxz-open-two-string.md) targets a low-lying complex-root
singlet on long chains (`--singlet-excitation`), retaining finite string
deviations in logarithmic coordinates; it does not enumerate singlets.

For `H=+sum(S.S)^2`, [ferromagnetic excitations](biquadratic-ferromagnetic.md)
reuse the same root equations with reversed physical energy ordering and
the exact degenerate ground reference. The full one-defect module and global
positive gap are analytic, for odd and even lengths. General module minima
are not obtained by reversing a truncated AF list: real-root scans can miss
complex-root minima, while Q-system coverage remains budget-dependent.
The [bound-pair solver](biquadratic-bound-pairs.md) targets the ell=N-4
two-string family directly on odd/even long chains, retaining both signs of
the finite deviation. The [three-string solver](biquadratic-bound-triples.md)
adds ell=N-6 droplets with a complex finite-size deviation, sharing the Newton
driver and reporting code. Neither enumerates the scattering spectrum.
[Real-root scattering windows](biquadratic-scattering.md) reuse the existing
equations with a bounded high-label range, avoiding a chain-length-sized
enumeration. Direct gaps also preserve ordering when total energies round
to the same value. Selected states and bounded [pair-plus-real-root scans](biquadratic-pair-defect.md)
reuse the two-string solver; combined with the other two three-defect
families, it reproduces the entire ell=N-6 module spectrum for N=6,...,10.
The mixed CLI uses the same cluster tables and native-precision exports,
with candidate budgets and explicit partial-convergence reporting.
The library and CLI extend to one pair plus several selected real roots,
with shared budgeted combination scans; four-/five-defect examples are
validated as ED subsets, not full spectra.
A [two-bound-pair solver and CLI scans](biquadratic-two-pairs.md) add another
four-defect subset, retaining both signed logarithmic deviations and sharing
the single-pair equations, energy mapping and Newton driver.
A [three-string plus real-root solver and CLI scans](biquadratic-triple-defect.md)
reuse the isolated triple and a general external-scattering kernel;
the [four-string droplet solver and CLI](biquadratic-bound-quartets.md)
completes the five four-defect topologies, with joint module-spectrum
validation through N=10 (not a general completeness proof).

Literature: [Barber–Batchelor](../CITATIONS.md#barber-batchelor-1989),
[Albertini](../CITATIONS.md#albertini-2000), and
[Aufgebauer–Klümper](../CITATIONS.md#aufgebauer-klumper-2010).

The [physical SU(2) decomposition](biquadratic-spin-content.md) of open spin-1
TL multiplicity spaces is available in the shared library and optional
`--spin-content` table, separately from auxiliary XXZ spin. **Next slice:**
higher-defect complex-root families and state-dependent spectral weights.
Odd free-end chains require a one-spinon/domain-wall branch, not a doubled
dimer ground state. PBC requires twisted XXZ references and different
representation bookkeeping; general TL, RSOS and boundary-algebra models
are not automatically covered by this first spin-chain implementation.

### `richardson`: finite pairing spectra without a spatial chain

The reduced BCS Hamiltonian has specified single-particle levels and a
uniform pair-scattering coupling. Richardson's solution and its Gaudin
generalizations are surveyed by the model developers in
[Dukelsky–Pittel–Sierra](../CITATIONS.md#dukelsky-2004).
This would broaden the library to pairing benchmarks at fixed particle number.

Implemented in [richardson.hpp](../include/bethe/richardson.hpp) and
`bethe-richardson` for g>=0, distinct ascending single-particle energies,
one time-reversed doublet per level, and specified pair/blocked sectors.
The diagonal pair-scattering term is included. The
[guide](richardson.md) derives the energy convention and regularized
quadratic equations from [Faribault et al.](../CITATIONS.md#faribault-2011).

Adaptive ground-state continuation uses a number-constrained rectangular
QR correction in native precision, with branch-distance and energy-bound
checks. Tests include all small pair sectors of irregular/clustered levels,
blocked levels, analytic limits and the exact four-level g=2/3 root collision.
Incomplete results explicitly retain their reached coupling; they are not
reported as energies at the requested target.

**Next slice:** optional rapidity recovery, excited occupation seeds, and
the derivative equations needed for repeated levels/higher degeneracies.
PBC/OBC and lattice momentum are not appropriate interface concepts here.

### `gaudin-magnet`: a central-spin realization

Implemented in [central_spin.hpp](../include/bethe/central_spin.hpp) and
`bethe-central-spin`: `H=B*S0^z+sum A_j*S0.Sj`, all spins 1/2, distinct
nonzero A_j of either sign, fixed total Sz, and any finite central field
including B=0. There are no bath fields/interactions. See the
[guide](central-spin.md) for the rational Bethe equations, normalization
and compactified eigenvalue-variable continuation.
Implementation checkpoint:
[`8ff2f52`](https://github.com/Uni20-dev/bethe/commit/8ff2f52).

This reuses the regularized-equation and native QR machinery developed for
Richardson, but uses a central-spin-specific high-field ground seed.
The connected star's sign-gauged exchange matrix gives a noncrossing
sector-ground branch. Tests compare all small sectors with independent
spin-basis diagonalization, check the two-spin formula and original
one- and two-root rational equations, and cover zero-field multiplets, clustered
couplings, larger baths and all three scalar types. Incomplete solves
report only a reached-field energy; the infinite-field seed has none.
The implementation checkpoint passed 747 tests with GCC 13 and fp128,
511 with Clang 20 Release without MPLAPACK, and the new CLI/citation checks
in an apps-only build using the published Uni20 pin. An additional 180
fixed-seed random small-sector checks agreed with independent diagonalization.

**Next slice:** repeated couplings/grouped spins, decoupled zero-coupling
bath spins, excited high-field seeds or observable matrix elements.
This does not implement arbitrary Gaudin Hamiltonians or central-spin
dynamics from the reference paper.

### `multicomponent-gas`: arbitrary-rank repulsive fermions

Implemented in [su_fermions.hpp](../include/bethe/su_fermions.hpp) and
`bethe-sun-fermions-pbc` for equal masses, equal positive contact couplings
`2c delta`, periodic circumference ell, and fixed component populations.
Every occupied interacting population must be odd; c=0 and single-component
limits allow arbitrary counts. Empty components are retained as metadata
but removed from the nested problem. The [guide](su-fermions.md) derives the
label parity and equations from [Sutherland](../CITATIONS.md#sutherland-1968)
and [Lee et al.](../CITATIONS.md#lee-2011).

Strong-to-weak continuation solves all real-root seas in native precision,
with analytic Jacobians and relative weak-cluster residuals at every level.
Incomplete solves retain only a previously converged coupling, never a
seed energy or an intermediate energy mislabeled as the requested c.
Tests include rational equations, the Gaudin–Yang reduction, the one-particle-
per-component Lieb–Liniger equivalence, free/weak/strong limits, a separate
plane-wave Hamiltonian for populations 3,1,1, and six-component 30-particle
sweeps. The spin-chain SU(3) frontend remains a different model.
The implementation checkpoint passed 774 tests with GCC 13 and fp128,
530 with Clang 20 Release without MPLAPACK, and the new CLI/citation checks
in an apps-only build using the published Uni20 pin.

**Next slice:** other finite-ring shell branches and excited seas/strings.
The distinct equal-coupling [Bose–Fermi mixture library](bose-fermi.md) now
implements odd-fermion ground shells using its own graded equations, not
by treating bosons as another fermionic color.

### `integrable-ladder`: singlet/triplet sector competition

Implemented in [ladder.hpp](../include/bethe/ladder.hpp) and
`bethe-ladder-pbc`: periodic L>=2 rungs, spin-1/2 leg coefficient 1,
four-spin coefficient 4, any finite rung exchange J_r, no field.
The [guide](ladder.md) gives the exact SU(4) permutation identity,
energy constants, shifted finite-ring labels and state-selection rules.
The rung coupling is a chemical potential for singlets; all triplet
populations are minimized within each specified singlet-count sector.

Three nested real seas are solved in native precision for the packed
ground branches of compatible Young diagrams, including displaced seas
and SU(4) descendants. This extra bookkeeping is essential: already at
L=6, the populations 4,1,1,0 have a lower descendant from the 4,2,0,0
multiplet than their own highest-weight sea. Reports preserve this
distinction rather than claiming descendant rapidities are finite.
Incomplete scans retain candidate upper bounds only. The all-singlet
product for J_r>=4 has a direct analytic ground-state path.

Tests compare all singlet-count sectors through seven rungs against
independent permutation matrices and small literal spin-basis ladder
Hamiltonians, as well as the original rational equations, momentum,
SU(2)/SU(3) reductions and algebraic native-precision regressions.
This is not a complete excited-state solver or a general nested-Bethe
completeness proof; budgets cover the entire multiplet scan.
The implementation checkpoint passed 795 tests with GCC 13 and fp128,
545 with Clang 20 Release without MPLAPACK, and the new CLI/citation checks
in an apps-only build using the published Uni20 pin. A 32-rung fp128
global scan also converged within the default budgets.

**Next slice:** field-dependent triplet populations or explicit excited
branches, with new finite-ring state-selection checks. The ordinary
two-leg Heisenberg ladder at generic couplings is not integrable here.

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
  The [negative-anisotropy engine](xxz-negative.md) now covers
  `-1<Delta<0` sector ground roots for even rings and open chains, with
  rank-subtracted equations scaled to remain discriminating near Delta=-1.
  Both paths are exposed through the ground-state APIs and existing frontends,
  with explicit hyperbolic coordinates and scaled residuals. The periodic
  ground result is distinct from the real-excitation type; negative odd rings
  remain internal and are deferred until a concrete calculation needs them.
  The [polynomial equations](xxz-polynomial.md) and an internal
  [momentum-constrained adaptive driver](xxz-odd-continuation.md) now handle
  conjugate-pair formation without extracting roots, with independent
  spin-basis energy checks through 21 sites at selected couplings.
  A separate [quantum-Wronskian diagnostic](xxz-wronskian.md) tests generic-q
  consistency; nonzero projected spin helices cover the phantom collision
  without incorrectly applying the finite-root theorem there.
  An [all-sector candidate scan](xxz-odd-continuation.md#comparing-sectors-without-assuming-the-answer)
  now compares every folded sector and withholds the minimum if any solve
  fails, with all-sector ED checks through 17 sites at selected couplings.
  A [regular-state criterion](xxz-regularity.md) now tests the Gaudin
  nonzero-vector hypotheses directly in coefficient space.
  [Projected spin helices](xxz-spin-helix.md) are available as explicit
  eigenstates at commensurate couplings; a separate diagnostic recognizes
  the all-phantom odd-ring collision.
  The [mixed-phantom reduction](xxz-phantom.md) checks the twisted finite
  equations after endpoint deflation; general nonzero lifting remains open.
  A [coordinate dressing map](xxz-phantom-wave.md) now constructs lifted
  amplitudes and exposes the nonzero-lift issue with explicit kernel tests.
  A [finite-root amplitude evaluator](xxz-coordinate-wave.md) supplies its
  input using subset sums, with native-precision direct Hamiltonian checks.
  [Native root recovery](polynomial-roots.md) now connects the continued and
  reduced polynomials to explicit vectors, including tested mixed-phantom
  lifts with three and four finite roots. It reports unresolved clusters
  instead of treating them as distinct roots.
  A bounded [nonzero-amplitude diagnostic](xxz-phantom-check.md) now adds
  numerical mixed-phantom witnesses to the sector scan, with propagated
  root-uncertainty estimates and visible work limits.
  Independent free-sea and projected-helix variational upper bounds now reject
  continuation energies incompatible with a sector minimum; passing these
  necessary checks is not a proof of minimality.
  **Deferred:** do not extend this internal investigation or expose it as a
  public ground-state solver without a concrete need. Preserve the code and
  tests; the [restart checklist](xxz-negative.md#deferred-odd-ring-work-restart-checklist)
  records the branch-selection, singular-state, numerical-failure, and API
  issues. This work does not block other models. Negative-Delta excitation
  classification remains a separate later extension.
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
  Half-filled repulsive [thermodynamic elementary lines](hubbard-dispersion.md)
  and [doped zero-field lines](hubbard-doped.md) are now separate implemented
  functionality, with both interaction conventions and Fermi energy references.
  Multiparticle threshold minimization remains a follow-up; the elementary
  holon is not automatically a sector minimum.

## Wider catalogue: useful, but not the next default targets

The following entries deliberately remain short until selected for an
equation-level feasibility study. Each has a concrete restriction and an
initial deliverable, rather than an unqualified claim of model support.

| ID | Literature and integrable restriction | Useful first deliverable / main obstacle |
| --- | --- | --- |
| `gaudin-magnet` | Rational Gaudin magnets; a concrete central-spin realization is treated by [Faribault–Schuricht](../CITATIONS.md#faribault-schuricht-2013) | [Central-spin sector ground states implemented](central-spin.md); general Gaudin charges, repeated couplings and higher spins remain open, not arbitrary spin-bath interactions |
| `multicomponent-gas` | Equal-mass SU(kappa) delta fermions, [Lee et al.](../CITATIONS.md#lee-2011); equal-mass Bose–Fermi mixture with equal repulsive Bose–Bose/Bose–Fermi couplings, [Imambekov–Demler](../CITATIONS.md#imambekov-demler-2006) | [Odd-population fermion sectors](su-fermions.md) and the [odd-fermion Bose–Fermi ground-state library](bose-fermi.md) implemented; other shells remain open. A trapped local-density calculation would be an approximation, not an exact trapped BA solution |
| `integrable-ladder` | [Wang's ladder](../CITATIONS.md#wang-1999), with its required exchange and four-spin terms | [Periodic zero-field sector minima implemented](ladder.md); fields and excitations remain open. Not the ordinary two-leg Heisenberg ladder at generic couplings |
| `q-boson` | Deformed boson hopping and its phase-model limit, [Bogoliubov–Izergin–Kitanine](../CITATIONS.md#bogoliubov-1997) and [Pozsgay](../CITATIONS.md#pozsgay-2014-q-boson) | [Fixed-N PBC ground-state and excitation library/frontend implemented](q-boson.md), with small-sector full-spectrum checks, free/phase limits and continuum scaling. Not the standard Bose–Hubbard chain |
| `xyz` | Zero-field spin-1/2 XYZ / eight-vertex family, [Baxter](../CITATIONS.md#baxter-1973) and [Zhang–Klümper–Popkov](../CITATIONS.md#zhang-klumper-popkov-2024) | [Native-precision even-chain regular ground solver and frontend implemented](xyz.md), checked against small-chain exact energies and XXZ limits; singular solutions and excited spectrum remain. Sz is not generally conserved |
| `haldane-shastry` | Spin-1/2 inverse-chord-square exchange on a ring, [Haldane](../CITATIONS.md#haldane-1988) and [Shastry](../CITATIONS.md#shastry-1988) | [Implemented](haldane-shastry.md): exact motif energies, momenta, Yangian dimensions and spin-sector minima. Full motif scans have an explicit budget; SU(2) decomposition, wavefunctions and correlations remain future work |
| `sutherland` | The trigonometric inverse-square gas, [Sutherland](../CITATIONS.md#sutherland-1971), in the Calogero–Sutherland family | [Implemented](sutherland.md): periodic scalar bosons on the specified Jastrow collision branch, exact ground and explicit-label energies, and bounded label-window spectra. Spinful/statistics variants, wavefunctions and other collision domains remain future work |
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

The [grand-canonical Lieb–Liniger library](lieb-liniger-thermal.md) now implements
this first TBA slice, with native-precision quadrature and independent nonlinear,
mesh, and cutoff checks. Pressure derivatives, the zero-temperature limit, and
the Tonks fugacity series provide validation. Canonical fixed-density inversion
and a temperature-scan frontend are also available. This model-specific layer
is not a universal TBA engine, and thermodynamic answers must never silently replace
finite-system eigenstates.

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
