# Finite four-spinon spectra

The experimental scanner attempts a bounded set of periodic XXX four-spinon
triplets, then attaches their physical raising weights to discrete spectral
lines. It uses the [selected finite-pair solver](experimental-xxx-four-spinon.md)
and its stabilized form factors. Every requested label remains in the result,
including unsupported geometry and numerical failures. Completing enumeration
does not establish a complete four-spinon spectrum: singular and collapsed
strings are outside this solver's supported geometry.

## A first line, then a momentum slice

From a build directory outside the source tree:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_APPS=ON
cmake --build . --target bethe-experimental-xxx-four-spinon
./experimental/bethe-experimental-xxx-four-spinon 12 \
  --holes=-3,-1,0,2 --branch 1 --csv line.csv
```

This solves the singlet ground state and one finite triplet afresh. Its absolute
state momentum index is 9, while the spectral momentum-transfer index is 3:
q=pi/2. The default channel is zz; use `--channel raising` for S-+.

To attempt all four-hole labels at that momentum, remove the explicit holes:

```sh
./experimental/bethe-experimental-xxx-four-spinon 16 \
  --momentum-index 4 --threads 2 --json slice.json \
  --csv-table spectrum=lines.csv --csv-table diagnostics=states.csv
```

Some requested weights can be unavailable. The tool then writes its accepted
lines and diagnostics and exits with status 2. A rejected line is absent from
the spectrum; it is not replaced by zero or estimated from nearby states.
A validated q=0 state is different: spin symmetry fixes its weight to exactly
zero. Such states remain accepted lines in the spectrum and diagnostics and
contribute zero to the moments; they do not cause exit status 2.
Without `--momentum-index`, the catalogue includes every momentum. The default
10,000-candidate budget rejects a larger catalogue before solving the ground
or opening exports. Explicit `--holes` and `--branch` avoid constructing a full
catalogue; they can also be checked against a supplied momentum index.

The executable is an opt-in build-tree tool and is not installed with the
supported frontends. Native and enabled wide precisions use the usual
[`--precision` option](command-line.md). MPFR scratch and workers inherit the
prepared singlet's finite precision; no ambient default is used.

## Labels, momentum and normalization

N is divisible by four and at least eight. Four distinct integer holes lie in
[-N/4,N/4], and the two-string branch is -1, 0 or 1. The catalogue has
3 binomial(N/2+1,4) requests before momentum filtering. These are quantum-number
requests, not a count of converged nonsingular physical states. At N=12 there
are 105 requests.

The Hamiltonian and Fourier convention are

```math
H=\sum_j\mathbf S_j\cdot\mathbf S_{j+1},\qquad
S_q^+=\frac{1}{\sqrt N}\sum_j e^{-iqj}S_j^+,\qquad
q=P_0-P_n\pmod{2\pi}.
```

For these N, the singlet has P0=0. A state's index is the hole sum minus its
branch modulo N; the Fourier transfer has the opposite sign. The
`--momentum-index` option and exported `momentum_index` always refer to the
transfer, not the absolute state momentum. The fixed-momentum catalogue solves
the constraint for the last hole, avoiding an all-momentum enumeration first.

The raising weight is N times the squared local matrix element. The spectral
function is 2 pi times the weighted delta-function sum. SU(2) relates the zz
weight to half the raising weight. No broadening or sum-rule rescaling is
applied. Every explicit label has multiplicity one: reflection partners and
energy-degenerate states remain distinct.

Stable text identities such as `N12-s4-h-3,-1,0,2-b1` identify canonical labels
across selections. Numeric `LehmannLine::state_id` values instead index records
within the returned scan, starting at one. Neither identity is an energy rank.

## Tables and missing states

| Table | Contents |
| --- | --- |
| `spectrum` (primary) | Accepted candidate ID, momentum transfer, gap, weight and multiplicity |
| `diagnostics` | Every requested ID and hole/branch labels; solve and factor status; acceptance, signed pair deviation, its logarithm, residual and pivot indicator; rejection reason |
| `moments` | Retained zeroth/first moments in all N bins and the full first-moment rule |
| `roots` | Prepared ground and converged excited standard ABA lambda roots |

`--diagnostics` and `--roots` show those tables on screen. Exporting a named
table requests its file without enabling the corresponding screen output.
Unknown table names are errors. CSV/TSV exports include the common provenance
preamble; JSON carries metadata and the tables selected for screen output.
Request named CSV/TSV exports to save diagnostics or roots without displaying them. The report records the channel,
normalization, precision, request count and missing-weight count.

For a momentum-filtered scan, unrequested bins of the retained-moment table
are zero because no lines were supplied there. Those zeros do not estimate
the physical spectrum. Ratios to the full first-moment rule measure the
accepted subset, not the missing part. A pivot margin or root residual is
a rejection diagnostic, not a certified form-factor error bar.

Rounded complex roots may lose an exponentially small pair deviation next to
one half. The separate `deviation` and `log_abs_deviation` columns preserve
the solver's representation; do not reconstruct it from `lambda_imag` alone.

## Library composition

```cpp
#include <bethe/experimental/xxx_four_spinon_scan.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace finite = bethe::experimental::xxx;
auto ground_state = bethe::heisenberg::ground_state<double>(12);
finite::PreparedSinglet<double> ground(12, ground_state);
auto requests = finite::four_spinon_catalogue(12, 3, 1000);
uni20::async::TbbScheduler scheduler(2);
auto pending = finite::scan_four_spinons_async(ground, requests, {&scheduler});
auto const& result = pending.get_wait(scheduler);
// Inspect result.records as well as result.lines; result.converged() only
// means every requested weight was obtained, not that the full DSF was found.
```

The synchronous `scan_four_spinons` observes that same final result. The Async
operation owns its prepared ground snapshot, canonicalized requests and
per-state scratch. The borrowed scheduler must outlive completion. Its bounded
continuous runners fuse solve and weight evaluation; no intermediate batch
barrier waits for the slowest preceding state. Final records and moment sums
follow canonical label order regardless of completion order. Duplicate labels
and mismatched ground/request sizes are rejected before submission.

Cancellation stops new work and lets active kernels finish before throwing
`ExecutionCancelled`. An unexpected exception likewise joins admitted work
before propagating. Neither outcome publishes partial results as a successful
scan. Dropping an Async handle does not cancel the operation. The implementation
waits on the returned result and never drains the whole scheduler with `run_all`.

Tests compare small catalogues against independent brute-force enumeration,
physical line weights against the selected-state arbitrary-precision controls,
and serial/concurrent output across native and runtime precisions. Controlled
held-state tests check final-result lifetime during failures and cancellation.
These checks validate the stated bounded interface, not general string coverage.

Source provenance: the momentum-constrained catalogue comes from
`tests/research_xxx_sector_scan.hpp` at research revision
`2242fd5d78ec2c6c81a64993eca5dd6eea6cd732`. This port replaces its historical
batch/backend choices with the supported continuous-work execution and the
reviewed finite-deviation evaluator. No research cache or archive is required.
