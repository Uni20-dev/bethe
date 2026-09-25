# Periodic q-boson states

The C++ library `bethe/q_boson.hpp` implements the fixed-particle-number ground
state and real-root excited states of the repulsive q-boson hopping model. The
frontend `bethe-q-boson-pbc` exposes ground states and excitation scans. It is **not** the
ordinary Bose–Hubbard model.

## Command-line use

```sh
build/bethe-q-boson-pbc 16 --particles 8 --eta 0.5 --roots
build/bethe-q-boson-pbc 16 --particles 8 --phase --json phase.json
build/bethe-q-boson-pbc 2 --particles 2 --eta 1e-60 --precision fp128
build/bethe-q-boson-pbc 5 --particles 3 --eta 0.5 --excitations all --roots --json spectrum.json
build/bethe-q-boson-pbc --references
```

The positional argument is the number of sites L. `--particles N` is required
and may exceed L; vacuum N=0 is allowed. Choose either finite `--eta` (including
zero) or `--phase`. An infinite `--eta` is rejected in favor of the explicit
phase option. Energies are in units of the specified unit hopping.

`--roots` adds a `roots` table to the `states` table. The labels $I$ are decimal
half-integers in exports. Free bosons have coincident zero momenta; interacting
ground roots are strictly ordered. On failure, energy and momenta are missing
(JSON null, CSV/TSV empty), while the exact requested labels remain available.
Exit codes are 0 for convergence, 2 for a failed numerical solve, and 1 for
invalid input or output errors. Invalid physical input is checked before output
files are opened, including with `--force`.

`--tolerance` and `--max-iterations` set the solver controls described below.
The shared [output options](output.md) provide the presentation-layer overview,
CPU timing, metadata and JSON/CSV/TSV exports. To export multiple tables, use
`--json result.json` or explicit destinations such as
`--csv-table states=states.csv --csv-table roots=roots.csv` with `--roots`.
References are displayed only with `--references`.

`--excitations COUNT|all` scans all `binomial(L+N-1,N)` canonical states,
retaining the lowest COUNT converged levels (including the ground state), or
all of them. COUNT does **not** limit the number of solves. The default
`--max-candidates 10000` rejects larger families before solving or opening
output files. This option requires `--excitations`; there is no padding window.

Scans add a `reference` table for the ground state and, on failure, a `failed`
table containing the first failed candidate. `states` contains only converged
levels, with energy gaps relative to the ground reference. If that reference
fails, gaps are null. Any failed candidate makes the scan incomplete (exit 2),
even if all requested retained levels were found. Metadata records candidate,
converged, and retained counts. State IDs are unique across these tables;
`roots` refers to those IDs and includes integer mode $m$ and native-precision
`deviation` as well as $I$ and lifted $k$. On failure both $k$ and `deviation`
are null. Default ground-state mode still emits one state, with gap zero on
success.

## Hamiltonian and conventions

For L>=2 periodic sites and N>=0 particles, use lattice spacing and hopping
amplitude one:

```math
\begin{aligned}
H&=-\sum_j\left(B_j^\dagger B_{j+1}+B_{j+1}^\dagger B_j-2N_j\right),\\
B\lvert n\rangle&=\sqrt{[n]_q}\,\lvert n-1\rangle,\\
[n]_q&=\frac{1-e^{-2\eta n}}{1-e^{-2\eta}},\qquad q=e^\eta.
\end{aligned}
```

Here $N_{j}$ counts particles; it is not `B_j^dagger B_j`. The $+2N$ diagonal
shift is included. To compare with a Hamiltonian containing only the hopping
terms, subtract `2N` from the returned energy. For L=2 the periodic sum contains
both bonds, doubling the hopping between its two sites.

This is the normalization in [Pozsgay (2014)](../CITATIONS.md#pozsgay-2014-q-boson),
equations (2.1), (2.13), and the energy immediately below (2.13). It is twice
the Hamiltonian used in [Bogoliubov–Izergin–Kitanine](../CITATIONS.md#bogoliubov-1997),
equation (1.1). Neither chemical potential nor external potential is included.

At eta=0, $[n]_q =n$: all ground-state momenta and the shifted energy are zero.
At eta=+infinity, $[n]_q =1$ for n>0: this is the *phase model*, not hard-core
bosons (multiple occupancy remains allowed). Its ground momenta are exactly
$k_{j} =2\,\pi \,I_{j} /(L +N)$ with consecutive $I_{j} =j -(N -1)/2$.

## Bethe equations and state selection

The finite-eta logarithmic equations used here are

```math
\begin{aligned}
Lk_j+\sum_{l\ne j}\theta(k_j-k_l)&=2\pi I_j,\\
\theta(d)&=2\operatorname{atan2}\!\left(\sin(d/2),\tanh\eta\cos(d/2)\right),\\
E&=\sum_j4\sin^2(k_j/2).
\end{aligned}
```

Ordered ground roots lie inside $(-\pi,\pi)$. Their differences may exceed pi;
the `atan2` branch must not be replaced by a principal `atan(tan(...))`.
Integers label odd N, half-odd integers even N. The centered consecutive labels
give total momentum zero. The phase-model formula follows from the same branch
and agrees with (3.4) in Bogoliubov–Izergin–Kitanine.

Our state-selection reasoning is separate from checking a root residual:
the connected occupation-basis hopping matrix has nonpositive off-diagonals
and a unique fixed-N ground state. The centered branch connects continuously
to its free-boson ground state. The logarithmic Jacobian is L times the identity
plus a positive weighted graph Laplacian, so this ordered branch has no local
singularity for eta>0. Small-sector independent diagonalization checks that
the selected energy is the lowest, not merely an eigenvalue.

## C++ API and numerical behavior

Excited states are selected by sorted free-boson mode labels
`0 <= m_0 <= ... <= m_(N-1) < L`, with
$I_{j} =m_{j} +j -(N -1)/2$. At zero deformation these are precisely bosonic occupation
patterns, with momenta $2\,\pi \,m_{j} /L$. The canonical family has
`binomial(L+N-1,N)` candidates, equal to the fixed-N occupation-space dimension.
The roots are lifted continuously, not individually folded into a Brillouin
zone: they remain ordered with total spread less than $2\,\pi$ for eta>0.
In the phase limit,
$k_{j} =(2\,\pi \,I_{j} +P_{\mathrm{lifted}})/(L +N)$, where `P_lifted=2*pi*sum(m_j)/L`.
The returned physical `momentum` is reduced to $(-\pi,\pi]$, using integer mode
sums before floating-point conversion.

```cpp
#include <bethe/q_boson.hpp>

auto state = bethe::q_boson::ground_state(16, 8, 0.5); // sites, particles, eta
if (state.converged) {
    auto energy = *state.energy;
    auto const& roots = state.momenta;
}

std::vector<std::size_t> modes{0, 0, 2};
auto excited = bethe::q_boson::solve_modes<double>(5, modes, 0.5);
auto scan = bethe::q_boson::real_excitations(5, 3, 0.5,
    {.count = 10, .max_candidates = 1000});
```

`solve_real` accepts the corresponding `QuantumNumbers` instead of modes.
The shared excitation scanner solves the entire canonical family and retains
up to `count` lowest converged levels, including the ground state and distinct
degenerate states. Set `count` to the candidate count to retain everything.
`excitation_count(L,N,max_candidates)` checks the combinatorial budget before
any solves or root allocation; this scan is intended for small sectors.
Failed candidates are counted and one diagnostic example is retained. Gaps
are absent if the ground reference failed. A converged candidate family is a
numerical result, not by itself a proof of spectral completeness.

Use native `long double` or enabled Uni20 fp128 types in place of double.
The optional fourth argument is `SolverOptions<Real>`, with residual tolerance
$32\,\epsilon$ and 10000 accepted Newton updates by default. A zero budget only
checks the seed; exact free/vacuum/one-particle/phase results can still succeed.
Positive infinity explicitly requests the phase model; negative eta, NaN,
L<2, and nonpositive/nonfinite tolerances are invalid.

The existing physical-domain Newton driver, native dense solve, and compensated
sums are shared with the continuum solvers. The phase kernel uses `tanh(eta)`
and scaled trigonometric factors rather than exponentially large hyperbolic
functions. Roots are solved as $k_{j} =2\,\pi \,m_{j} /L +u_{j}$; `modes` and `deviations`
retain weakly split clusters even when rounded `momenta` coincide. For small
eta, the known rank phase is canceled against the labels before arithmetic.
The residual is scaled by $L \,\max (\sqrt{\tanh (\eta)/L },\lvert u_{j} \rvert)$, retaining a
relative test as deviations approach zero. The positive energy uses sine
squares instead of subtracting $2-2\,\cos (k)$.

`iteration_limit`, `stalled`, and `precision_limit` leave `energy` empty. Roots
on failure are diagnostics, not a converged eigenstate. Nonzero but insufficiently
resolved subnormal energies are rejected. Wavefunctions, correlation functions,
OBC, and q<1 are not supported.

Tests compare full energy multisets and joint energy/cos(momentum) spectra
against occupation-space diagonalization for L=2..5 and N=0..4, including free
and phase limits and several finite deformations. They also cover a native
two-site/two-particle formula, the original complex multiplicative equations,
the exact phase limit, weak eta=epsilon², larger rings, continuum scaling to
Lieb–Liniger, and explicit validation/iteration/precision failures.
