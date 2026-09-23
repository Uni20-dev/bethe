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
```

## Normalization

We use spin-1/2 operators, lattice spacing one, `J=1`, zero field, and

```text
H = (pi/N)^2 sum_{i<j} S_i.S_j / sin^2(pi*(i-j)/N).
```

Each unordered pair occurs once, including at `N=2`. The fully polarized
energy is `E_F=pi^2*(N^2-1)/(24*N)`. Ground energies are

```text
even N: E0 = -pi^2*(N^2+5)/(24*N)
odd  N: E0 = -pi^2*(N^2-1)/(24*N).
```

Even rings have one singlet ground state. Odd rings have two opposite-momentum
spin-1/2 multiplets, four states in total; the default output includes both.
These are finite-ring formulas, not thermodynamic approximations.

## Motifs, energies and counting

An admissible motif is an increasing list of positions `m` in `[1,N-1]`,
with no adjacent occupied positions. For example, `1,3,6` is allowed at
`N=8`, while `1,2` is not. An empty motif is the ferromagnetic multiplet.
The spectral rules and Yangian counting are reviewed in
[Jiang–Lamers–Miao, Secs. 2.2–2.3](https://arxiv.org/html/2606.20168v2).
Their Hamiltonian uses `(1-P_ij)/(4*sin^2)`; convert by
`H_ours=E_F-2*(pi/N)^2*H_theirs`.

For a motif with `M` entries:

```text
E = E_F - (pi/N)^2 sum_m m*(N-m)
momentum_index = sum_m m modulo N
P = 2*pi*momentum_index/N
spinons = N-2*M; S_max = (N-2*M)/2.
```

One row denotes a **Yangian multiplet**, generally containing several ordinary
SU(2) multiplets. `S_max` is its maximum total spin, not a claim that every
state in that row has that spin. The full dimension is
`(N+1)` for the empty motif, otherwise
`m_first*(N-m_last)*product(neighboring_gap-1)`.
For example, motif `{2}` at `N=4` contains a singlet and a triplet: its
dimension is four and `S_max=1`, not `3/2`.

## Ground sectors and excitation scans

`--sz` selects an allowed half-integer spin projection. The minimizing motif
has `M=N/2-|Sz|` entries, packed with spacing two as symmetrically about
`N/2` as the integer lattice allows. Odd rings have two placements unless
`M=0`. Returned dimensions still describe the **whole** Yangian multiplet,
not just the selected projection. Every reported gap references the global
ground energy, including when selecting a nonminimal spin sector.

`--levels COUNT` returns the lowest COUNT motifs, including ground motifs.
`--levels all` enumerates every motif and thus describes the entire finite
Hilbert-space spectrum with multiplicities. Equal energies remain separate
rows; COUNT can cut through a degeneracy. Ordering uses exact integer energy
coefficients, then momentum index, then the motif list.

The number of motifs grows as `F_(N+1)`. Before enumeration we check it
against `--max-motifs` (default 100,000). If the budget is insufficient,
exit status 2 reports the refusal and publishes no alleged lowest levels.
Ground-sector and specified-motif calculations do not require a spectrum scan.
The explicit size range `2<=N<=1,000,000` keeps energy numerators in signed
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
independently constructed spin-basis Hamiltonian through `N=8`, check every
magnetization-sector minimum through `N=9`, and verify that motif dimensions
sum to `2^N`. Native-precision tests cover even/odd formulas, exact ordering,
invalid labels, count overflow and budget refusal. Gaps are formed from exact
integer differences before floating conversion.

Wavefunctions, SU(2) decomposition of each Yangian multiplet, correlations,
open-chain variants and a thermodynamic spinon frontend are not implemented.
