# Haldane–Shastry spin ring

`bethe-haldane-shastry-pbc` evaluates exact finite-ring spectral rules for the
spin-1/2 antiferromagnetic inverse-chord-square chain. This is not the
nearest-neighbor XXX model: there are no nonlinear root solves or convergence
tolerances. Arithmetic supports fp64, long-double and optional fp128.

```sh
build/bethe-haldane-shastry-pbc 16
build/bethe-haldane-shastry-pbc 15 --sz 3/2
build/bethe-haldane-shastry-pbc 8 --motif 1,3,6
build/bethe-haldane-shastry-pbc 12 --levels 10 --precision fp128
build/bethe-haldane-shastry-pbc 8 --levels all --format csv
build/bethe-haldane-shastry-pbc 6 --motif 3 --spin-content --format json
```

## Normalization

We use spin-1/2 operators, lattice spacing one, `J=1`, zero field, and

```math
H=\left(\frac\pi N\right)^2\sum_{i\lt j}
\frac{\mathbf S_i\cdot\mathbf S_j}{\sin^2[\pi(i-j)/N]}.
```

Each unordered pair occurs once, including at $`N =2`$. The fully polarized
energy is $`E_{F} =\pi ^{2}\,(N ^{2}-1)/(24\,N)`$. Ground energies are

```math
E_0=\begin{cases}
-\dfrac{\pi^2(N^2+5)}{24N},&N\text{ even},\\{}
-\dfrac{\pi^2(N^2-1)}{24N},&N\text{ odd}.
\end{cases}
```

Even rings have one singlet ground state. Odd rings have two opposite-momentum
spin-1/2 multiplets, four states in total; the default output includes both.
These are finite-ring formulas, not thermodynamic approximations.

## Motifs, energies and counting

An admissible motif is an increasing list of positions $`m`$ in $`[1,N -1]`$,
with no adjacent occupied positions. For example, `1,3,6` is allowed at
$`N =8`$, while `1,2` is not. An empty motif is the ferromagnetic multiplet.
The spectral rules and Yangian counting are reviewed in
[Jiang–Lamers–Miao, Secs. 2.2–2.3](https://arxiv.org/html/2606.20168v2).
Their Hamiltonian uses $`(1-P_{\mathrm{ij}})/(4\,\sin ^{2})`$; convert by
$`H_{\mathrm{ours}} =E_{F} -2\,(\pi /N)^{2}\,H_{\mathrm{theirs}}`$.

For a motif with $`M`$ entries:

```math
\begin{aligned}
E&=E_F-\left(\frac\pi N\right)^2\sum_m m(N-m),\\{}
\mathrm{momentum\_index}&=\sum_m m\pmod N,\\{}
P&=\frac{2\pi}{N}\,\mathrm{momentum\_index},\\{}
\mathrm{spinons}&=N-2M,\qquad S_{\max}=\frac{N-2M}{2}.
\end{aligned}
```

One row denotes a **Yangian multiplet**, generally containing several ordinary
SU(2) multiplets. `S_max` is its maximum total spin, not a claim that every
state in that row has that spin. The full dimension is
$`(N +1)`$ for the empty motif, otherwise
`m_first*(N-m_last)*product(neighboring_gap-1)`.
For example, motif `{2}` at $`N =4`$ contains a singlet and a triplet: its
dimension is four and $`S_{\mathrm{max}} =1`$, not `3/2`.

### Resolving ordinary total-spin multiplets

The library supplies `spin_decomposition(N, motif)`. The CLI's optional
`--spin-content` adds an SU(2)-resolved table alongside the Yangian levels.
Remove the sites in `motif` and `motif+1` from the ordered list `1,...,N`.
Each remaining consecutive run of length l contributes one spin-l/2 factor.
The Yangian multiplet's SU(2) content is the tensor product of these factors,
decomposed by the usual Clebsch–Gordan rule. This is the Drinfeld-polynomial
string construction in [Jiang–Lamers–Miao, Eqs. (2.16)–(2.18)](https://arxiv.org/html/2606.20168v2#S2.SS3).
These strings label representation factors, not bound-state Bethe strings.

For example, at N=6 motif `{3}` leaves two runs of length two, so its content
is spin-1 tensor spin-1: one singlet, one triplet and one quintet. All nine
states have the motif's energy and momentum. An empty motif instead gives
one spin-N/2 irrep, even for a very large ring.

```cpp
auto content = bethe::haldane_shastry::spin_decomposition(6, {3});
if (content.complete) {
  for (auto const& term : content.multiplets) {
    // term.spin is a native half_int; term.multiplicity counts SU(2) irreps.
    // Each such irrep contributes one state to every allowed Sz in [-S,S].
  }
}
```

The reusable `bethe/spin_multiplets.hpp` helper owns the integer spin coupling;
the Haldane–Shastry module only constructs its motif-dependent factors.
`bethe::spin::tensor_product(factors, options)` accepts nonnegative `half_int`
spins and returns increasing-spin multiplicities. Its empty product is a
singlet. There is no floating-point rounding in this calculation.

`DecompositionOptions::max_updates` defaults to 1000000 and counts individual
Clebsch–Gordan multiplicity additions, independently of motif-enumeration
budgets. On `work_limit` or `count_overflow`, `complete` is false and no
partial spin decomposition is published. Counts use uint64; a total dimension
can overflow even when each multiplicity fits, in which case `dimension` is
absent but the completed decomposition remains valid. Negative factors and
unrepresentable sums of twice-spin are invalid arguments.

The CLI `spin_content` table links to `levels` through `state_id`. Each row
contains `spin`, `multiplicity`, `updates`, `complete` and `status`. The work
count is per motif and is repeated across its spin rows; do not sum these
repeated values. Even when `--sz` selects a ground-sector motif, the content
describes its **whole** Yangian multiplet, not only the selected projection.
`--levels COUNT` still counts motifs, not SU(2) irreps.

`--max-spin-updates` sets the addition budget per motif and requires
`--spin-content`. On failure, one diagnostic row with missing spin and
multiplicity replaces that motif's decomposition. Other motifs are still
attempted, all valid energy rows remain available, and the overall outcome
is partial (exit 2). On motif-enumeration refusal, both tables are empty.
JSON preserves both named tables; for separate rectangular CSV/TSV files use:

```sh
build/bethe-haldane-shastry-pbc 8 --levels all --spin-content \
  --csv-table levels=levels.csv --csv-table spin_content=spins.csv \
  --json spectrum.json
```

## Ground sectors and excitation scans

`--sz` selects an allowed half-integer spin projection. The minimizing motif
has $`M =N /2-\lvert S^z \rvert`$ entries, packed with spacing two as symmetrically about
$`N /2`$ as the integer lattice allows. Odd rings have two placements unless
$`M =0`$. Returned dimensions still describe the **whole** Yangian multiplet,
not just the selected projection. Every reported gap references the global
ground energy, including when selecting a nonminimal spin sector.

`--levels COUNT` returns the lowest COUNT motifs, including ground motifs.
`--levels all` enumerates every motif and thus describes the entire finite
Hilbert-space spectrum with multiplicities. Equal energies remain separate
rows; COUNT can cut through a degeneracy. Ordering uses exact integer energy
coefficients, then momentum index, then the motif list.

The number of motifs grows as $`F _{N +1}`$. Before enumeration we check it
against `--max-motifs` (default 100,000). If the budget is insufficient,
exit status 2 reports the refusal and publishes no alleged lowest levels.
Ground-sector and specified-motif calculations do not require a spectrum scan.
The explicit size range $`2\le N \le 1,000,000`$ keeps energy numerators in signed
64-bit arithmetic and bounds single-motif storage. A multiplet dimension that
exceeds unsigned 64-bit range is unavailable, not wrapped; its energy remains
valid. This is separate from scalar rounding of energies and momenta.

`auto`, `pretty`, `plain`, `csv`, `tsv`, and `json` are available. Delimited output
has `#` metadata, CPU time and a rectangular table; the motif cell is a
space-separated list (or `empty`). The named `levels` table includes zero-based
`state_id`; half-integer spins use decimal values such as `0.5` in exports.
Independent exports, e.g. `--csv levels.csv --json levels.json`, and streaming
use the shared [output options](output.md). Use `--references` for literature,
also collected in [CITATIONS.md](../CITATIONS.md).

## Library and checks

```cpp
#include <bethe/haldane_shastry.hpp>
namespace hs = bethe::haldane_shastry;
auto ground = hs::ground_levels<double>(15);
auto sector = hs::sector_ground_levels<double>(16, uni20::half_int{2});
auto selected = hs::evaluate<long double>(8, {1, 3, 6});
auto all = hs::spectrum<double>(12); // Check all.complete.
```

Tests compare the complete motif spectrum and momentum cosines against an
independently constructed spin-basis Hamiltonian through $`N =8`$, check every
magnetization-sector minimum through $`N =9`$, and verify that motif dimensions
sum to $`2^N`$. Native-precision tests cover even/odd formulas, exact ordering,
invalid labels, count overflow and budget refusal. Spin-content tests compare
the complete magnetization-resolved energy and momentum-cosine spectra with
independent spin-basis diagonalization through N=8. Across all motifs through
N=16, the multiplicities also reproduce the SU(2) decomposition obtained by
counting spin words and subtracting adjacent magnetization-sector dimensions.
Gaps are formed from exact
integer differences before floating conversion.

Wavefunctions, correlations,
open-chain variants and a thermodynamic spinon frontend are not implemented.
