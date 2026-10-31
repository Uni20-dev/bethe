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
- The two-spinon thermodynamic heat map also uses binary64. Its introductory
  single-point command did not need a different arithmetic recommendation.
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

## Four-spinon starting point

The [compact comparison](../tutorials/data/xxx-four-spinon-precision.json)
retains all five formats, including binary64's numerical failure, and can be
checked or regenerated with `scripts/check_xxx_four_spinon_precision.py`.
It fixes q=1.5, omega=1.8, J=1, order4/T8, density rtol/atol=1e-8/1e-30,
64/80 residue terms and 16/20 tail orders. The two-spinon kernel tolerance is
explicitly 1e-12 in every format. One worker is used throughout.

| Arithmetic | Density comparison / guard bits | Result | Query CPU seconds |
| --- | ---: | --- | ---: |
| Binary64 | 24 / 53 | Arithmetic screen fails; no S4 or sum | 0.0026 to failure |
| x86 long double | 53 / 64 | Complete | 0.337 |
| Binary128 | 53 / 113 | Complete | 6.69 |
| MPFR, 128 bits | 64 / 128 | Complete | 48.2 |
| MPFR, 192 bits | 128 / 192 | Complete | 64.0 |

Hardware: Intel Core i7-12700KF, GCC 13.3 Release, source `b573fb4`, Uni20
`1cb2474`. Timings include the query and ordered output on a shared host; they
are not portable performance guarantees. Failure time is not the cost of a
complete binary64 calculation. The native 64-bit-significand S4 differs from
the MPFR192 value by about 1.14e-16 relative; MPFR128 and MPFR192 differ by
about 2.42e-38. These comparisons describe the **same truncated quadrature**,
not the unrestricted physical response.

The native full-grid check is separate: orders 2/4 pass at T8 and T12, but
orders 8/16 fail the arithmetic screen at both cutoffs. Thus the cheap starting
command is not an adequate uniform setting for the refinement figure. The
[refinement comparison](four-spinon-refinement-precision-20261031.json) also
records binary128 failing at order 8 on both cutoffs: its lower comparison
remains at 53 bits. MPFR128 raises that lower comparison to 64 bits and passes
all eight points, agreeing with the original MPFR192 S4 values to within
3.4e-38 relative. The tutorial figure/reproducer now uses this measured
128-bit grid (41 CPU minutes total; the original used 54). Sources and two-spinon
kernel tolerance differ between those grids, so this is a comparison of the
unchanged S4 calculation and measured costs, not an identical two-sector
benchmark. No failed point is replaced by a successful prefix or a looser
acceptance threshold.

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

## Scope still being checked

The small examples and large-domain resolved slice need separate
recommendations. The resolved-slice pilots are still checking their larger
domains, quadrature and sampling requirements. This inventory does not replace
model-specific accuracy tests for every parameter in the library.
The [precision guide](../numerical-accuracy.md) gives the decision procedure and
explains the relative/absolute screen and demonstrated higher-precision case.
