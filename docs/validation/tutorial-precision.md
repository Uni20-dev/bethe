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

The native full-grid check is separate: orders2/4 pass at T8 and T12, but
orders8/16 fail the arithmetic screen at both cutoffs. Thus the cheap starting
command is not an adequate uniform setting for the refinement figure. The
existing MPFR192 reference grid remains available while the lower-cost MPFR
refinement checks are evaluated. No failed point is replaced by a successful
prefix or a looser acceptance threshold.

## Scope still being checked

The six-spinon owner is validating the identical three-frequency small-demo
streams across native and MPFR formats. Its finite R8 demonstration and the
large-domain resolved slice are different calculations and need separate
recommendations. Higher-order four-spinon MPFR comparisons are also in progress.
The [precision guide](../numerical-accuracy.md) gives the decision procedure and
explains the relative/absolute screen and demonstrated higher-precision case.
