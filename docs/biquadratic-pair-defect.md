# A bound pair scattering with real defects

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Real-root scattering](biquadratic-scattering.md) · [Bound triples](biquadratic-bound-triples.md)

Three TL singlet insertions need not all bind, or all remain separate. A
two-string plus one real root describes the intermediate **pair-plus-defect**
branch in ell=N-6. Selected states and CLI family scans support this branch
on odd/even N>=6, in fp64, native long double and enabled fp128.
`--excitations` still means the purely real-root family; it does not contain
these mixed-string levels.
The same CLI also supports one pair plus several real roots, as described
[below](#a-pair-with-several-real-defects).

## Command line and bounded scans

```sh
# One explicitly selected pair of Bethe labels:
build/bethe-biquadratic-obc 128 --ferromagnetic --pair-defect 125,123 --roots
# Complete small mixed family, not the entire three-defect module:
build/bethe-biquadratic-obc 8 --ferromagnetic --pair-defects all
# Lowest eight converged levels among the 36 selected label combinations:
build/bethe-biquadratic-obc 129 --ferromagnetic --pair-defects 8 --mixed-window 6 \
  --precision long-double --json mixed.json --csv mixed.csv \
  --tsv-table labels=mixed-labels.tsv
```

`--pair-defects COUNT|all` scans the full rectangle
$`1\le I \le N -3, 1\le J \le N -5`$ by default. `--mixed-window WIDTH` instead selects
the highest WIDTH labels on **each** axis: `I=N-2-WIDTH,...,N-3` and
`J=N-4-WIDTH,...,N-5`. This gives WIDTH squared candidates and requires
`1<=WIDTH<=N-5`. A high-label window makes long-chain low-energy searches
practical, but does not certify that excluded labels have higher energy.

`COUNT` retains at most that many converged levels **after** solving every
candidate in the selected rectangle. `all` retains every converged candidate
in that rectangle, not the full excited spectrum. The `--max-candidates`
budget (default 10000) limits the work before allocation and checks the
product without overflowing. Requesting one output level does not bypass
the scan budget. A selected `--pair-defect I,J` requires no enumeration.
These options require `--ferromagnetic` and exclude other state selectors.
The default single-real-root family fixes ell=N-6. `--mixed-window` applies
only to `--pair-defects`.

Sorting uses the direct gap, not the extensive total energy: the latter
can round to the same value for distinct states on long chains. Exact
numerical gap ties use ascending `(I,J)`. These are rankings within the
computed set, not general module minima or a global excitation ranking.

The shared cluster reporter supplies `states`, `reference`, `string`,
`labels`, and optional `roots` tables. JSON exports all tables; `--csv`
and `--tsv` export `states`, while `--csv-table NAME=FILE` and
`--tsv-table NAME=FILE` select auxiliary tables. Metadata accompanies the
exports, and `--no-retain` uses the same streaming table interface. The
solver still stores the budgeted candidate set for sorting.

`labels` records I, J and alpha in one row per real root, grouped by
`state_id` and ordered by I;
`string` records the pair center, sign, L and optional representable
deviation. `roots` lists the real roots first, followed by the conjugate
pair. The `states.row` column is a display row number, **not a physical
mode label**. L and the sign retain the deviation when exported rounded
roots no longer resolve it; an underflowed deviation is null, not exact zero.

If any candidate fails, the command returns exit status 2 and records
the total converged/failed counts. Retained converged levels come first,
followed by **one unconverged diagnostic sample**, whose gap is null.
Its energy and TL energy are estimates only. That extra diagnostic row
does not count toward COUNT. Missing candidates could lie below the retained
ones: a partial result is not a verified lowest-level list. The exact ground
reference remains valid even if every candidate fails.

## Library API

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;

std::size_t const n = 128;
auto state = ferro::pair_defect<long double>(n, n - 3, n - 5);
if (state.reference.converged) {
    auto gap = state.tl_energy; // E-E0, E0=N-1; evaluated directly
    auto energy = state.energy;
    auto alpha = state.reference.rapidities[0];
    auto a = state.reference.center;
    auto L = state.reference.log_deviation;
}
```

The two arguments after N are integer Bethe labels **I,J**, not energy
ranks or physical lattice momenta:

```math
I=1,\ldots,N-3\quad\text{(real root)},\qquad
J=1,\ldots,N-5\quad\text{(bound pair)}.
```

The high-label corner $`(N -3,N -5)`$ approaches the separated-cluster threshold
`5/3 + 1 = 8/3`. This is above the three-defect droplet edge 2 and below the
three-unbound-defect edge 3. It is not the global first excitation, and the
labels do not imply a universal ordering of finite-size levels. All of these
counts refer to TL insertions, **not physical spin flips**. Operator spectral
weights and physical-spin decomposition remain separate questions.

## Same equations, different labels

The auxiliary model is quantum-group XXZ with opposite end fields at
Delta=3/2, not the usual zero-field open XXZ chain. We reuse the
[two-string equations](xxz-open-two-string.md) with one real root:

```math
\begin{aligned}
v&=\frac{i\alpha}{2},\qquad u_\pm=\frac{\eta+d}{2}\pm\frac{ia}{2},\\{}
\eta&=\mathrm{arcosh}\Delta,\qquad d=\mathrm{sign}\,e^{-L},\\{}
0&\lt \alpha,a\lt \pi,\qquad |d|\lt \eta.
\end{aligned}
```

Replace the real-root right-hand side by $`2\,\pi \,I`$, and the pair-phase
right-hand side by $`2\,\pi \,J`$. The modulus equation is unchanged. Retaining
the unsquared pair equation fixes `sign(d)=(-1)^(N-J-1)`; the product-phase
equation alone loses this sign information. The real root may cross the
pair center: no additional ordering condition $`\alpha \lt a`$ is imposed.
In the continuous logarithmic convention, the alpha=pi endpoint has
I=N-2 and the a=pi endpoint has J=N-4; those excluded endpoint roots give
the label bounds above.

These are our regularized equations derived from
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020), not an ideal-string
energy formula. L and the separate deviation sign remain authoritative
when exp(-L) underflows. The existing analytic Jacobian, damped Newton
driver, iteration budget, native precision and direct energy-shift mapping
are shared with the AF singlet and isolated bound-pair solvers. There are
three unknowns, so storage and work per Newton step do not grow with N.

`reference.real_labels` stores I and `reference.string_label` stores J.
The common cluster result's `mode` is zero for this mixed branch: there is
no one-dimensional mode ordering. An unconverged result retains estimates
and residuals, **not verified energies**. No tolerance relaxation, hidden
extra iterations or change of precision is performed.

The lower API
`bethe::xxz::quantum_group::two_string::pair_defect(N,Delta,I,J,options)`
accepts finite Delta>1. Branch collapse or convergence failure can occur;
the biquadratic Delta=3/2 coverage below is not a guarantee for other Delta.

## Checks and limits of coverage

For N=6,...,10, independent exact diagonalization verifies every level of
the ell=N-6 module when the three families are combined:

| Family | Candidate count |
| --- | ---: |
| Three positive real roots | binomial(N-3,3) |
| One pair plus one real root | (N-3)(N-5) |
| One three-string | N-5 |

Their sum is the module dimension `binomial(N,3)-binomial(N,2)`. Matching
the dimension alone would not establish completeness: the regression
compares the **entire sorted spectrum including multiplicities**. This is
small-chain numerical evidence, not a completeness proof for arbitrary N.

Further tests check the original complex equations for both deviation signs
and odd/even chains, the analytic Jacobian, native-precision six-site
characteristic polynomial, zero/one-iteration failure diagnostics, and the
8/3 threshold up to N=100000. The generic finite-size string description
still permits numerical failure and does not provide physical SU(2)
labels, scattering amplitudes or form factors.

## A pair with several real defects

The library and CLI also accept **one pair and a selected real-root set**,
with M=r+2 insertions and ell=N-2M.

```sh
# Two real roots I1=123, I2=124, and pair label J=121:
build/bethe-biquadratic-obc 128 --ferromagnetic --pair-defect 123,124,121 --roots
# All six one-pair-plus-two-real-root candidates on eight sites:
build/bethe-biquadratic-obc 8 --ferromagnetic --pair-defects all --real-defects 2
# Two real roots in a six-label window, times six pair labels: 90 candidates:
build/bethe-biquadratic-obc 129 --ferromagnetic --pair-defects 8 \
  --real-defects 2 --mixed-window 6 --precision fp128 --json mixed-four.json
```

For selected states, `--pair-defect I1,...,Ir,J` infers r from the labels;
the **last** label belongs to the pair. For scans, `--real-defects R`
sets r (default 1). Both require r>=1 and M=r+2<=N/2. Use `--bound-pairs`
for an isolated pair; the library additionally accepts an empty real-root set.

Without a window, the candidate count is
`binomial(N-M,r)*(N-2M+1)`. With `--mixed-window WIDTH`, it is
`binomial(WIDTH,r)*WIDTH`: choose r distinct real labels from
`N-M-WIDTH+1,...,N-M`, independently of the pair label in
`N-2M-WIDTH+2,...,N-2M+1`. The window requires
`r<=WIDTH<=N-2M+1`. Near the largest allowed M, no such shared-width window
may exist; omit the window to scan the full family within the candidate
budget, or select explicit labels. No unbudgeted fallback scan occurs.

All candidates are solved before selecting COUNT converged levels. Sorting
uses direct gaps, with lexicographic `(I1,...,Ir,J)` ties, and the same
failure-sample convention as the r=1 scan. The shared combination iterator
and overflow-checked binomial count are also used by the real-root-only
excitation scans. Table schemas are unchanged: several `labels` rows now
share each state ID, and the root table has r+2 rows per state.

```cpp
std::vector<std::size_t> labels{123, 124};
auto four_defects = ferro::pair_with_real_roots<long double>(128, labels, 121);
// M=4, ell=120: one pair plus two separate defects, not a four-string.
if (four_defects.reference.converged) {
    auto gap = four_defects.tl_energy;
}
```

The domain is $`2\le M \le N /2`$, ordered distinct real labels in `1,...,N-M`,
and pair label `J=1,...,N-2M+1`. The endpoint limits of the same logarithmic
equations give these ranges: a real root at alpha=pi would have label
N-M+1, while the pair center at a=pi would have J=N-2M+2. Both endpoints
are excluded. The deviation sign remains $`(-1)^{N -J -1}`$.

No new Bethe equations or Newton implementation are introduced. Explicit
sea labels, a label-dependent pair seed and edge-scaled coordinates extend
the existing two-string system. For fixed r the work is independent of N;
in r it uses a dense (r+2)-by-(r+2) Jacobian, quadratic storage and cubic
linear-solve work per update. Matrix-size checks precede label allocation.
Conditioning, branch topology and precision can still prevent convergence.

An empty label set reproduces an isolated pair. One label reproduces
`pair_defect`, which delegates to this API. Consecutive labels
`1,...,N/2-2` and J=1 reproduce the even-chain AF singlet's auxiliary
solution, though the general API uses a different Newton scaling. Its
biquadratic adapter always maps to the **ferromagnetic** convention.
The lower XXZ API is
`two_string::pair_with_real_roots(N,Delta,labels,J,options)`.

For Delta=3/2, all pair-plus-two-real-root candidates match distinct module
ED levels for N=8,...,10 (6, 20 and 45 candidates respectively). The ten
pair-plus-three-real-root candidates at N=10 also match ED. These are
**subsets**, not complete four- or five-defect spectra. A separate
[two-pair solver and CLI scans](biquadratic-two-pairs.md) cover another four-defect
family. A selected [triple-plus-one-real-root solver](biquadratic-triple-defect.md)
covers another; triples with more real roots and larger droplets remain missing.
Independent original-equation checks cover r=2,3,4, both deviation signs,
odd/even chains and all supported precisions. High-label branches approach
the separated-cluster threshold $`5/3+r`$, checked up to N=100000. This
threshold is not a claim about the minimum of the whole TL module.
