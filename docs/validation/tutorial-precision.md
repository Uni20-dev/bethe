# Tutorial precision audit

This ledger links the tutorials' precision choices to measured comparisons.
Each record identifies its numerical controls, source revision and hardware;
agreement at those inputs does not establish one accuracy bound for every
model or plotted point. See the [precision guide](../numerical-accuracy.md)
for the decision procedure and the
[development study](../spinon-numerics-development.md) for implementation history.

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
- Four-spinon refinement uses binary128; the small six-spinon example uses
  native extended80 where available. Their measured comparisons are below.
  The resolved-slice tutorials have [separate recorded checks](#scope-of-the-resolved-slice)
  on their much larger soft domains.

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
compare successive precisions with non-arithmetic controls fixed, then refine
the mesh, state family or domain separately. The spinon comparisons below measure
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

### Changing the momentum or frequency

The [native neighbourhood record](four-spinon-neighbourhood-precision-20261031.json)
checks 49 supported points and five physical zeros at the same order and cutoff.
With a 53/64-bit pair, extended80 completes 26/49 supported points at `1e-8`,
47/49 at `1e-6`, and 49/49 at `1e-4`. Binary64 completes none at these
tolerances. All outcomes are retained. Wider references at three difficult
points agree with accepted native S4 values within 1.3e-16 relative; that
comparison does not certify every point at the relaxed tolerances.

The [binary128 neighbourhood record](four-spinon-native-neighbourhood-20261031.json)
uses a 64/113-bit pair at all 49 supported points. It completes 48/49 at `1e-8`
and 49/49 at `1e-6`; the strict failure at q=1.5, omega=4.5 remains unavailable.
All 97 successful estimates agree with matched MPFR128 S4 references within
1.3e-30 relative. The binary128 exports identify source `865d701`, a Xeon Gold
6254 and two workers per query. The MPFR exports lack an embedded revision;
the record retains the review author's source attribution. Worker counts and
sector selections differ, so the comparison concerns S4 values, not total-query
cost. Both sets preserve numerical controls and raw hashes. The checker
reconstructs the compact coverage and comparisons without reading raw archives:

```sh
python3 scripts/check_xxx_four_spinon_native_neighbourhood.py
```

The starting command uses `1e-6`, with a margin below the integration
changes it demonstrates; its actual two-sector rerun is checked separately.
The strict five-format and refinement comparisons retain `1e-8`. Changing the
momentum/frequency can still require a stronger precision pair. A successful
relaxed screen needs a separate comparison with a stricter reference.

### Refining the quadrature

The [refinement checks](four-spinon-refinement-precision-20261031.json) show
that extended80 passes orders 2/4 at T8 and T12, but fails the arithmetic screen
at orders 8/16. The cheap starting command therefore cannot support the whole
refinement figure. With the 64/113-bit pair, all eight binary128 quadratures
pass the strict `1e-8` screen and agree with MPFR192 within 6.3e-31 relative.
The [comparison record](spinon-native-precision-20261031.json) retains their
controls and total cost of 348 CPU seconds on an i7-12700KF. **The tutorial
figure and reproducer use these binary128 data.** Failed cases remain in the
records; no failed node is replaced by a successful prefix.

The same record validates the small six-sector API spectrum and moment commands
in six formats. Native extended80 agrees with MPFR384 within 4e-16 relative
across values, replicate means/maxima, sampling errors and covariance;
binary128 agrees within 3e-30. Both binary64 failures are retained. Those small
API examples recommend extended80 where available, with measured CPU costs
of 0.036 and 0.040 seconds respectively. Their tiny samples and finite cutoff
limit the physical estimates separately.

## Six-spinon starting point and difficult geometry

The [eight-format comparison](six-spinon-precision-20261031.json) fixes the
three-frequency R8 example, including all 384 draws and replicate identities.
Native x86 long double succeeds in 0.31 seconds wall / 0.38 CPU seconds and
agrees with the 384-bit reference means to within 3e-16 relative. MPFR128
also succeeds, costing 43.6 CPU seconds; MPFR384 costs 101.1. Binary64 fails
arithmetic screens and an exponent-range check. MPFR96 fails its 32/96-bit
arithmetic pair. The record preserves these failures rather than publishing
partial estimates. Independent reconstruction checks successful means,
sampling errors and work counts. The large sampling errors and omitted R8
domain remain; precision agreement does not remove them.

### Changing the sample stream

The [seed comparison](six-spinon-native-comparison-20261031.json) retains
complete R8/4x32 queries for seeds 42, 2 and 35. Extended80 fails the arithmetic
screen at seeds 2 and 35. Binary128 with a 64/113-bit pair completes all three
at unchanged `rtol=1e-5`, `atol=0`, costing 7.32, 6.34 and 7.02 CPU seconds
on the same i7-12700KF. Every stream/work identity matches its MPFR192
reference; means, sampling SE and replicate means/maxima agree within 4.9e-30
relative. The checker is `scripts/check_xxx_six_spinon_native_precision.py`.
These measurements support trying binary128 on a failed extended80 stream;
they do not establish sampling or cutoff convergence. The development page
links the larger seed scans and additional tolerance experiments.

### Separating geometry from density precision

The [endpoint-preserving geometry comparison](six-spinon-endpoint-geometry-20261101.json)
checks R1040 with 256-bit geometry and a 192/256-bit density pair against a
1536-bit reference on identical 64-draw streams. Sample identities and support
counts match, and the mean changes by less than 8e-74 relative. Its roughly
100% relative sampling SE remains; this is arithmetic evidence on those inputs,
not a production precision floor or a bound on omitted soft weight.
The [development account](../spinon-numerics-development.md#a-representation-can-create-its-own-precision-requirement)
explains the representation and preserves the earlier failed geometry checks.

## Scope of the resolved slice

The small examples and large-domain resolved slice need separate
recommendations. Representative [four-spinon controls](four-spinon-slice-calibration-20261031.json)
and [six-spinon controls](six-spinon-slice-calibration-20261031.json) have been
checked on both sides of the two-spinon upper edge. They separate arithmetic,
contour, quadrature/sampling and domain comparisons. The complete
[sixteen-frequency comparison](../tutorials/xxx-spinon-slice.md#read-the-reference-comparison)
retains the independent sector estimates and their partial sums. All six-spinon
queries are complete, including precision retries; observed sampling standard
errors range from 2.49% to 6.52%, with fifteen points meeting the initial 5% aim.
The [reference record](six-spinon-reference-20261102.json) preserves the
individual errors, influence diagnostics, arithmetic stages and measured costs.

The [smaller installed example](six-spinon-native-retry-tutorial-20261101.json)
demonstrates the same retry workflow at R=256 in about 68 CPU seconds, or
20 seconds with four workers on the measured Xeon Gold 6254. It retains its
15.43% sampling standard error and separate domain limitation. The full
reference is optional; readers can redraw its compact data without integrations.
Neither the sampling bars nor the measured control changes bound omitted tails.
This inventory does not replace model-specific accuracy tests for every
parameter in the library.
The [precision guide](../numerical-accuracy.md) gives the decision procedure and
explains the relative/absolute screen and demonstrated higher-precision case.
