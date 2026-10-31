# Tutorial precision audit

Audit started 31 October 2026 against public main `4e776f6`. This records
precision choices and their evidence; it does not claim a single accuracy
bound for every model or every plotted point.

## Inventory

- All 654 checked-in tutorial CSV exports identify `fp64`. Their regenerators
  either request binary64 explicitly or use `capture_csv_tables`, which does so.
  This covers the ordinary finite-chain, gas, impurity, field-theory and
  dispersion examples. Their existing checks inspect successful statuses,
  model identities and the relevant analytical or numerical consistency tests.
- The two-spinon thermodynamic heat map and introductory point command both
  use binary64. The point command fixes the kernel tolerance at `1e-12`; its
  complete rerun returns the same S2 value as the five-format comparison below,
  differing from MPFR192 by less than 2e-15 relative. A successful two-spinon
  value does not imply that the four-spinon part of the same query succeeds.
- The 75 JSON Lines records in the one-magnon precision tutorial deliberately
  compare five formats. They are checked against an independent Decimal sine
  evaluation. The examples exhibit observable cancellation and root sensitivity,
  and explain the stable special-case formula; they are not ordinary defaults.
- Four-spinon refinement and the small six-spinon example originally used MPFR
  without a cheaper-format comparison. These are the main corrective work here.
  The separate resolved-slice tutorials must make their own checks on their much
  larger soft domains.

Most root-solving CSV examples use residual tolerances near a few multiples of
binary64 epsilon. Those are stopping criteria, not energy-error bounds. The
finite-volume sine-Gordon examples explicitly use a `1e-7` absolute tolerance
and report quadrature, source and tail diagnostics separately. Their physical
finite-size/asymptotic limitations remain distinct from arithmetic.

### What the ordinary tutorial checks establish

Binary64 is the ordinary plotting format, not a claim that every observable
has sixteen accurate digits. The model tutorials retain their own checks and
limitations. These examples show where to look before increasing precision:

| Examples | Checks on the saved calculations | Separate limitation or refinement |
| --- | --- | --- |
| [Finite XXX spectra](../tutorials/xxx-structure-factor.md#4-read-coverage-not-just-convergence) | State convergence and measured intensity/first-moment coverage | A converged selected family can omit spectral weight; broadening does not restore it. |
| [Hubbard dispersions](../tutorials/hubbard-doped.md#reproduce-and-extend), [Lieb–Liniger](../tutorials/lieb-liniger.md#5-reproduce-the-figures) | Background/row status, branch identities, units and limiting or symmetry relations | Refine the integral-equation mesh separately from arithmetic, especially near extreme densities or couplings. |
| [Lee–Yang](../tutorials/lee-yang.md#4-match-a-non-hermitian-mps-calculation), [sine-Gordon](../tutorials/sine-gordon.md#4-compare-with-mps-calculations) | Exported convergence diagnostics, mass/length scaling and available exact limits | Distinguish equation, source and tail errors from finite-volume approximations and lattice/MPS comparison errors. |
| [Richardson](../tutorials/richardson.md#5-reproduce-and-validate), [central spin](../tutorials/central-spin.md#5-reproduce-and-validate) | Regular-variable equations, constraints and energy reconstruction; optional independent small matrices | Check the physical observable near root collisions; an ill-conditioned root representation need not indicate a singular energy. |
| [Kondo](../tutorials/kondo.md#4-check-signs-units-and-numerical-errors) | Response bounds, derivatives, scale covariance and independent reference values | Integral/series and tail diagnostics are not rigorous enclosures; their units differ between magnetization and energy. |
| [Potts](../tutorials/potts.md#5-reproduce-and-validate), [non-Hermitian XXZ](../tutorials/xxz-nonhermitian.md#5-reproduce-and-validate) | Sector labels, dimensions, equations and independent small-matrix checks | Precision does not remove finite-size corrections or enlarge the selected state family. |
| [TASEP/ASEP](../tutorials/asep.md#5-reproduce-and-validate) | Analytic endpoints, rate/particle symmetries and an independent small generator | Follow the declared mode; a converged continuation is not an enumeration of every eigenvalue. |

These are existing model checks, not new wider-arithmetic comparisons for all
654 exports. For an observable requiring more accuracy than the tutorial plot,
repeat its calculation at fixed non-arithmetic controls before refining the
mesh, state family or domain separately. The spinon comparisons below measure
the precision and cost choices explicitly.

## Four-spinon starting point

The [compact comparison](../tutorials/data/xxx-four-spinon-precision.json)
retains all five formats, including binary64's numerical failure, and can be
checked or regenerated with `scripts/check_xxx_four_spinon_precision.py`.
It fixes q=1.5, omega=1.8, J=1, order4/T8, density rtol/atol=1e-8/1e-30,
64/80 residue terms and 16/20 tail orders. The two-spinon kernel tolerance is
explicitly 1e-12 in every format. One worker is used throughout.

| Arithmetic | Density comparison / guard bits | Result | Query CPU seconds |
| --- | ---: | --- | ---: |
| Binary64 | 24 / 53 | Arithmetic screen fails; no S4 or sum | 0.0069 to failure |
| x86 long double | 53 / 64 | Complete | 0.329 |
| Binary128 | 64 / 113 | Complete | 6.86 |
| MPFR, 128 bits | 64 / 128 | Complete | 51.0 |
| MPFR, 192 bits | 128 / 192 | Complete | 65.1 |

Hardware: Intel Core i7-12700KF, GCC 13.3 Release, source `0b496a0`, Uni20
`1cb2474`. Timings include the query and ordered output on a shared host; they
are not portable performance guarantees. Failure time is not the cost of a
complete binary64 calculation. The native 64-bit-significand S4 differs from
the MPFR192 value by about 1.14e-16 relative; MPFR128 and MPFR192 differ by
about 2.42e-38. These comparisons describe the **same truncated quadrature**,
not the unrestricted physical response.

The [neighbourhood audit](four-spinon-neighbourhood-precision-20261031.json)
extends this check to 54 points at the same order and cutoff. It independently
reproduces the review finding: excluding five physical zeros, long double
completes 26/49 points at `1e-8`, 47/49 at `1e-6`, and 49/49 at `1e-4`;
binary64 completes none at any of those tolerances. All 324 native outcomes
are retained. Five additional checks compare binary128/MPFR at three difficult
points: the complete MPFR128 references agree with accepted native values to
within 1.3e-16 relative, and MPFR192 confirms one reference. This is wider-value
evidence for those three points only. The recorded binary128 failure uses the
previous 53-bit lower comparison; this neighbourhood audit predates the
intermediate native comparison.

The later [binary128 neighbourhood record](four-spinon-native-neighbourhood-20261031.json)
retains all 98 queries at the same 49 supported points, with a 64/113-bit pair.
It completes 48/49 points at `1e-8` and 49/49 at `1e-6`; the strict failure at
q=1.5, omega=4.5 remains unavailable. All 97 completed estimates agree with
the matched MPFR128 S4 values within 1.3e-30 relative. These are clean-checkout
exports from `865d701` on a Xeon Gold 6254, with two workers per query. The older
MPFR exports retain their unavailable embedded revision and the review author's
source attribution. They request sectors 2,4 with different worker counts;
only the S4 value is compared, not total-query cost. Both sets retain complete
numerical controls and raw hashes. The checker reconstructs the coverage and
comparisons from the compact record; it does not claim to inspect raw archives:

```sh
python3 scripts/check_xxx_four_spinon_native_neighbourhood.py
```

The starting command now uses `1e-6`, with a margin below the integration
changes it demonstrates; its actual two-sector rerun is checked separately.
The strict five-format and refinement comparisons retain `1e-8`. Changing the
momentum/frequency can still require a stronger precision pair. The full native neighbourhood takes
17.5 CPU seconds and the five reference checks 81.8; these select only sector4,
so they are not directly comparable with the two-sector costs above.

The original 53/64-bit full-grid check is separate: orders 2/4 pass at T8 and T12, but
orders 8/16 fail the arithmetic screen at both cutoffs. Thus the cheap starting
command is not an adequate uniform setting for the refinement figure. The
[refinement comparison](four-spinon-refinement-precision-20261031.json) also
records binary128 failing at order 8 on both cutoffs with a 53-bit lower
comparison. MPFR128 raises that comparison to 64 bits and passes
all eight points, agreeing with the original MPFR192 S4 values to within
3.4e-38 relative. That MPFR128 grid cost 41 CPU minutes; the original MPFR192
grid cost 54. Sources and two-spinon kernel tolerance differ between those grids,
so this compares the unchanged S4 calculation and measured costs, not identical
two-sector benchmarks.

The [native follow-up](spinon-native-precision-20261031.json) uses the stronger
64/113-bit binary128 pair, with the same strict density tolerances and contour
controls. All eight quadratures now pass and agree with MPFR192 within 6.3e-31
relative. They cost 348 CPU seconds, or 5.80 CPU minutes. **The tutorial figure
and reproducer now use these binary128 data.** Earlier failed cases and MPFR
reference values remain in their original audit. No failed node is replaced by
a successful prefix or a looser acceptance threshold.

The follow-up also validates the exact small six-sector API spectrum and moment
commands in six formats. Native 80-bit agrees with MPFR384 to better than 4e-16
relative across values, replicate means/maxima, sampling errors and covariance;
binary128 agrees within 3e-30. Both binary64 failures are retained. These small
API examples now recommend native 80-bit on platforms providing that format.
Their tiny samples and finite cutoff still limit the physical estimates.

## Six-spinon starting point and difficult geometry

The [eight-format comparison](six-spinon-precision-20261031.json) fixes the
three-frequency R8 example, including all 384 draws and replicate identities.
Native x86 long double succeeds in 0.31 seconds wall / 0.38 CPU seconds and
agrees with the 384-bit reference means to within 2.7e-16 relative. MPFR128
also succeeds, costing 43.6 CPU seconds; MPFR384 costs 101.1. Binary64 fails
arithmetic screens and an exponent-range check. MPFR96 fails its 32/96-bit
arithmetic pair. The record preserves these failures rather than publishing
partial estimates. Independent reconstruction checks successful means,
sampling errors and work counts. The large sampling errors and omitted R8
domain remain; precision agreement does not remove them.

The subsequent [seed study](six-spinon-seed-precision-20261031.json) keeps all
40 predeclared native queries at the same R8/4x32 settings and `rtol=1e-5`.
Seeds 2 and 35 fail; the first failure, seed2, is replayed completely with
MPFR128/192 at unchanged tolerance. All stream/work identities match, and the
means agree within 4.35e-37 relative. A native `1e-3` replay agrees with that
strict192 reference within 3.25e-16 in the means and 7.08e-15 in replicate
means/maxima. Binary64 still fails, including exponent range. Independent
reconstruction checked all44 raw hashes, failure identities, successful means
and sampling errors. The printed seed42 example remains valid; its arithmetic
recommendation does not automatically transfer to another stream or sample
count. No passing seed is substituted for a failed one.

The [native binary128 follow-up](six-spinon-native-comparison-20261031.json)
repeats the exact three-frequency queries for seeds 42, 2 and 35 with the
64/113-bit pair and unchanged `rtol=1e-5`, `atol=0`. All three complete, costing
7.32, 6.34 and 7.02 CPU seconds respectively on the same i7-12700KF. Independent
reconstruction of all six native/reference exports verifies every compact
field, stream/work identity and replicate statistic. Means, sampling SE and
replicate means/maxima differ from MPFR192 by less than 4.9e-30 relative.
The corresponding checker is `scripts/check_xxx_six_spinon_native_precision.py`.
This supplies a measured native retry before MPFR for the two observed failures;
it does not extend the seed scan to forty binary128 queries or establish
sampling/domain convergence. Original failed native80 records remain intact.

A separate [geometry comparison](six-spinon-geometry-precision-20261031.json)
holds the density pair at 192/256 bits and probes R1040. Geometry at 256,
512 and 1024 bits is rejected before export: the soft endpoint rounds onto
pi. The same 64-draw stream succeeds at 1536 and 2048 bits, with matching
support/work identities and relative changes in replicate means and maxima
below 1e-450. Individual observation values were not retained in this probe.
This demonstrates a genuine precision requirement of the current coordinate
representation, without asserting that the physical estimate has that accuracy.
Only 64 draws were used; this is a geometry and arithmetic probe, not a
converged spectral calculation.

## Scope of the resolved slice

The small examples and large-domain resolved slice need separate
recommendations. Representative [four-spinon controls](four-spinon-slice-calibration-20261031.json)
and [six-spinon controls](six-spinon-slice-calibration-20261031.json) have been
checked on both sides of the two-spinon upper edge. They separate arithmetic,
contour, quadrature/sampling and domain comparisons. The sixteen-frequency
display calculation is still in progress; its final sampling quality and
figures remain to be assessed. The measured control changes are not bounds on
the omitted tails. This inventory does not replace model-specific accuracy
tests for every parameter in the library.
The [precision guide](../numerical-accuracy.md) gives the decision procedure and
explains the relative/absolute screen and demonstrated higher-precision case.
