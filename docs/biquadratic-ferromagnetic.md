# Excitations of the ferromagnetic biquadratic chain

[Overview](../README.md) · [Biquadratic/TL conventions](biquadratic.md) · [Output and exports](output.md)

The ferromagnetic option uses **H=+sum (S_i.S_(i+1))^2**, the opposite sign
to this tool's default. The spins are 1 and the physical chain has free ends.
There is no bilinear term or added constant. With $e_{i} =3\,P_{\mathrm{singlet}}$,

```math
H_F=(N-1)+\sum_i e_i,\qquad E_0=N-1,\qquad E-E_0=\sum_i e_i\ge0.
```

The ground space is known exactly. The focus here is the positive-energy
structure above it: the one-defect band, levels in other TL modules, and the
complex-root states that real-root scans miss. All modes support fp64,
native long double, and enabled fp128, with the shared screen/JSON/CSV/TSV output.

## A first exploration

```sh
# Complete one-defect band, analytically (odd N also supported):
build/bethe-biquadratic-obc 64 --ferromagnetic --one-defect
build/bethe-biquadratic-obc 65 --ferromagnetic --one-defect --json band.json
# Target bound-pair modes on a long chain (not all two-defect states):
build/bethe-biquadratic-obc 128 --ferromagnetic --bound-pairs 8
# Three-defect droplets, not all three-defect states:
build/bethe-biquadratic-obc 128 --ferromagnetic --bound-triples 8
# Mixed bound pair plus a separate defect, a 6-by-6 label window:
build/bethe-biquadratic-obc 128 --ferromagnetic --pair-defects 8 --mixed-window 6
# Real-root scattering near the two-defect low-energy edge:
build/bethe-biquadratic-obc 129 --ferromagnetic --through-lines 125 --excitations all --real-window 8
# Small two-defect module, including complex-root states:
build/bethe-biquadratic-obc 6 --ferromagnetic --q-spectrum --through-lines 2 --roots
# Compare it with the restricted real-root family in the same module:
build/bethe-biquadratic-obc 6 --ferromagnetic --excitations all --through-lines 2 --roots
# Selected real-root level; labels and roots retain the auxiliary XXZ convention:
build/bethe-biquadratic-obc 8 --ferromagnetic --quantum-numbers 7 --roots
```

Both scans default to $\ell =N -2$ in ferromagnetic mode. Set $\ell =N -4$ to
explore two defects, $\ell =N -6$ for three, etc. The number $M =(N -\ell)/2$
counts **TL singlet defects, not physical spin lowerings**. Through-lines
are not a physical SU(2) total-spin quantum number.
The diagrammatic through-lines are the unpaired strands; the remaining sites
form M non-crossing singlet arcs, so N=ell+2M. Some TL literature calls the
through-lines themselves "defects"; our M-defect terminology instead
counts singlet insertions M. These are algebraic sectors, not counts of
fixed local singlet bonds in an eigenstate.

Without state-selection options, `--ferromagnetic` reports only the exact
ground-space reference and the positive spectral gap. `--through-lines N-2`
alone selects the exact first positive level. For ell=N-4, use the targeted
`--bound-pairs` family, `--bound-triples` for ell=N-6, or `--bound-quartets`
for ell=N-8. These select bound branches, not a general module-minimum search; other nontrivial
module minima still require a Q-system investigation. `--sectors` and the AF-specific
`--singlet-excitation` are rejected with this sign, not silently reinterpreted.

## What the one-defect band means

Take a singlet on bond i and fully polarized spins elsewhere, and call the
state `|i>`, for `i=1,...,N-1`. These states span one copy of the $\ell =N -2$
module. They are independent but not orthogonal: adjacent states overlap by
1/3. Direct action of the singlet projectors gives

```math
(H_F-E_0)|i\rangle=3|i\rangle+|i-1\rangle+|i+1\rangle,\qquad
|0\rangle=|N\rangle=0.
```

Thus coefficients proportional to $\sin (\pi \,j \,i /N)$ give

```math
\begin{aligned}
k_j&=\frac{\pi j}{N},\quad j=1,\ldots,N-1,\\
\mathrm{gap}_j&=3+2\cos k_j,\qquad E_j=N-1+\mathrm{gap}_j.
\end{aligned}
```

The program emits decreasing j, hence increasing energy. Its `wave_number`
column is this **OBC standing-wave coordinate**, not a translation eigenvalue.
The gap is evaluated separately from the extensive energy; subtracting two
large printed total energies is not needed. These are exact expressions
evaluated at the chosen floating-point precision, not a Newton calculation.

Every TL eigenvector carries the spin-1 representation multiplicity
$m _{N -2}$, where $m_{0} =1$, $m_{1} =3$, and $m _{\ell +1}=3\,m_{\ell} -m _{\ell -1}$.
This includes other physical-spin states beyond the polarized-background
representative. It is not one SU(2) multiplet. If levels from other modules
coincide, their multiplicities must be added separately.

The **global positive gap**, not just a one-defect variational estimate, is

```math
\Delta_N=3-2\cos(\pi/N)\longrightarrow1.
```

This follows by TL equivalence from
[Koma-Nachtergaele, Proposition 2](../CITATIONS.md#koma-nachtergaele-1997),
rescaling their XXZ Hamiltonian by $2\,\Delta =3$. It is the gap above the
*entire* ground manifold. Physical one-spin-lowering states on the fully
polarized background remain inside that manifold and cost zero energy.
They are not the positive-energy band tabulated here.

The full band has N-1 rows. `--max-candidates` bounds its allocation (default
10000); raise it explicitly for larger bands. Analytic modes accept odd and
even N>=2. `--roots` is not available in analytic modes; use
`--excitations all --roots` for the equivalent one-root Bethe family.

## Why one band is not the full low-energy spectrum

For N=4, the ground energy is 3 with multiplicity 55. All the distinct
positive-energy levels, in increasing physical energy, are:

| Energy | TL defects M | Through-lines ell | Multiplicity | Root family |
| --- | --- | --- | --- | --- |
| $6-\sqrt{2}$ | 1 | 2 | 8 | One real root |
| $(15-\sqrt{17})/2$ | 2 | 0 | 1 | Complex roots |
| `6` | 1 | 2 | 8 | One real root |
| $6+\sqrt{2}$ | 1 | 2 | 8 | One real root |
| $(15+\sqrt{17})/2$ | 2 | 0 | 1 | Two real roots |

The two-defect singlet already interleaves the one-defect band. A real-root
scan of ell=0 returns only the **higher** singlet, even with `all`.
To see the missing level:

```sh
build/bethe-biquadratic-obc 4 --ferromagnetic --q-spectrum --through-lines 0 --roots
build/bethe-biquadratic-obc 4 --ferromagnetic --excitations all --through-lines 0
```

The Q-system's first state now has complex roots. `q_coefficients` and `roots`
use the reordered `state_id`; an exported polynomial can be supplied to
`--q-seed` to follow that branch. This is a practical entry point to studying
complex-root families, not an automatic classification of bound states.

For the selected two-string family, [the bound-pair solver](biquadratic-bound-pairs.md)
now follows modes directly on long chains, including both signs of the finite
string deviation. `--bound-pairs COUNT|all` fixes ell=N-4, supports odd and
even N>=4, and excludes scattering states. The lowest branch approaches
E-E0=5/3, below the two-separated-defect threshold 2.

The [three-string solver](biquadratic-bound-triples.md) adds three-defect
droplets in ell=N-6, with a complex finite-size deviation.
`--bound-triples COUNT|all` supports odd/even N>=6. Its lowest branch
approaches E-E0=2, below 3 for three separated defects and 8/3 for a
separated pair plus one defect. Neither targeted family enumerates scattering.

[Real-root scattering windows](biquadratic-scattering.md) select the high-I
edge of the existing real family with `--real-window WIDTH`. They make
few-defect long-chain scans practical without enumerating every combination.
`all` then means all combinations within the chosen window, not all
excited states; mixed string-plus-real-root branches remain excluded.

## Ordering, references, and coverage

Sign reversal leaves the auxiliary XXZ equations and roots unchanged.
It reverses physical energies and TL energies:

```math
E_F=-E_{\mathrm{AF}}=-2E_{\mathrm{ref}}+\frac{7(N-1)}4,\qquad
\mathrm{tl\_energy}=\sum_i e_i=E_F-(N-1).
```

The `reference_energy` column still denotes the unchanged auxiliary XXZ
energy. The separate `reference` table denotes the global **physical**
ground reference. Real-family scans select in ferromagnetic energy order
*before* retaining COUNT levels. Their exact vacuum reference cannot fail
because of a Newton iteration budget, although excited-state solves can.
Failed states have no published gap and scans exclude them from ranking.

Q-system searches include complex roots, retaining the existing budget,
admissibility and numerical completeness checks. Validation covers even
N<=8; larger systems are experimental. An incomplete search is a set of
discoveries, **not necessarily the lowest levels**, even after sorting.
Selected real roots and ferro real-family/window scans accept odd/even N,
as do the targeted pair/triple solvers, analytic band and ground space.
Q-system modes and AF ground/sector helpers still require even N.

Multiplicities use checked uint64 arithmetic, becoming unavailable on overflow
rather than approximate. This does not invalidate an energy or gap.
See [Zhou et al.](../CITATIONS.md#zhou-2025-biquadratic) for the ground-space
structure and Fibonacci degeneracies.

## Library and next steps

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;
auto first = ferro::one_defect_level<long double>(65, 64);
auto gap = ferro::spectral_gap<long double>(65);
auto real = ferro::real_excitations<long double>(6, 2, {.count=10});
auto complex = ferro::qsystem::spectrum<long double>(6, 2);
auto pair = ferro::bound_pair<long double>(128, 1);
auto triple = ferro::bound_triple<long double>(128, 1);
auto scattering = ferro::real_excitations_window<long double>(129, 125, 8);
auto mixed = ferro::pair_defect<long double>(128, 125, 123);
auto pair_two_real = ferro::pair_with_real_roots<long double>(128, std::vector<std::size_t>{123, 124}, 121);
auto two_pairs = ferro::two_bound_pairs<long double>(128, {121, 122});
// Inspect real.converged() and complex.complete() before interpreting coverage.
```

Targeted two-, three-, and four-defect bound families are available on long chains.
Real-root scattering windows and [pair-plus-defect states and CLI scans](biquadratic-pair-defect.md)
are also available. One pair with several real roots is supported through
selected labels or `--pair-defects ... --real-defects R` scans.
A [two-bound-pair solver and CLI scans](biquadratic-two-pairs.md)
are also available (`--two-pairs J1,J2` or `--two-pair-states COUNT|all`,
with optional `--pair-window WIDTH`). A selected-state
[triple-plus-real-root solver and CLI scans](biquadratic-triple-defect.md) add
another four-defect family (`--triple-defect I,J` or `--triple-defects COUNT|all`
with optional `--mixed-window WIDTH`). A selected-state
[bound-quartet solver and CLI](biquadratic-bound-quartets.md) completes the
five four-defect topologies, jointly checked against complete module ED
through N=10. Use `--bound-quartets COUNT|all` for the selected four-string
family. Higher-defect families remain open work.
The [physical-spin decomposition](biquadratic-spin-content.md) is now available
through the shared library and optional `--spin-content` table for multiplicity
spaces whose total dimension fits uint64.
Other open work is state-dependent spectral weights/form factors,
and periodic twists with genuine momentum labels. Energies and multiplicities
alone do not determine which branches an operator or a chosen iMPS ground
state couples to. None of those additional capabilities is implied here.
