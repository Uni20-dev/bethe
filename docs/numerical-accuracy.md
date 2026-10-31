# Choosing precision and numerical tolerances

[Command-line options](command-line.md#select-arithmetic-precision) ·
[Worked precision example](tutorials/xxx-precision.md)

Use the least expensive arithmetic that reliably resolves the observable you
need. More printed digits, a small equation residual, or a successful density
screen do not establish the accuracy of the final physical answer. First decide
which observable and absolute or relative accuracy matter for the comparison.

## Separate the sources of error

| Check | What to change | What it cannot settle |
| --- | --- | --- |
| Arithmetic and cancellation | Repeat the same calculation with wider arithmetic | Missing states, integration domains or mesh points |
| Root or integral-equation solve | Tighten stopping criteria and check the resulting observable | An ill-conditioned observable can amplify a small residual |
| Contour or series evaluation | Increase residue budgets, tail order or analytic remainder controls | Outer quadrature and sampling |
| Deterministic quadrature | Increase order or mesh resolution on the same domain | Weight outside that domain |
| Sampling | Increase independent work and inspect influential observations | A biased proposal implementation or omitted domain |
| Integration tails | Extend the domain separately, or validate a derived tail treatment | Arithmetic within the retained domain |

Start with an inner numerical budget comfortably below the change you wish to
resolve. Then demonstrate that tightening it does not change that observable
appreciably. There is no useful reason to demand thirty-digit inner accuracy
for a plotting example limited by percent-scale quadrature or sampling, unless
conditioning requires those working digits to obtain the modest final accuracy.

## What the spinon density tolerances mean

The W4/W6 recovery screen compares its numerical diagnostics with

```math
\mathrm{atol}+\mathrm{rtol}\,|W|.
```

For `--density-relative-tolerance 1e-8 --density-absolute-tolerance 1e-30`,
the relative term dominates when $`|W|\gg10^{-22}`$. This requests roughly
eight-digit relative agreement there, not thirty digits. Below that crossover,
the absolute floor dominates. It permits a large relative discrepancy in a
sufficiently small density. Choose that floor in the units and scale of the
**density being checked**, not the final plotted integral.

Jacobian factors, integration weights and repeated evaluations intervene between
W and the observable. The integrated arithmetic/contour diagnostics are useful
screening quantities, not certified total error bounds. Neither tolerance
accounts for omitted tails, quadrature error or Monte Carlo sampling error.
The two-spinon transition integral has its own `--tolerance`.

### Why evaluate a density at two precisions?

Cancellation and special-function evaluation can lose digits even when a
floating-point format has ample range. For each integration configuration, the
controller independently evaluates the same form-factor density at a lower and
a higher precision, using the same contour controls. The adapters convert the
retained geometry coordinates separately into each arithmetic context; they do
not merely round a previously computed density, nor recompute the outer geometry.
The higher-precision density supplies the result. The lower value is a numerical
check, not a second physical contribution or a value to average with the first.

Their absolute difference, both kernels' arithmetic indicators and the change
in the scalar prefactor form the arithmetic screen. That screen plus the higher
kernel's contour remainder must fit within the tolerance above. The calculation
also repeats at a larger contour budget: both controls must pass and their
values must agree within their combined indicators. These are stability checks,
not rigorous error bounds or estimates of the outer integration error.

| Selected arithmetic | Lower / higher density significand bits |
| --- | --- |
| Binary64 | 24 / 53 |
| x86 80-bit `long double` | 53 / 64 |
| Binary128 with an intermediate native `long double` | 64 / 113 on x86 |
| Binary128 without an intermediate native format | 53 / 113 |
| MPFR with N bits | max(24, N-64) / N |

Native selection uses actual significand widths, not the operating system or
processor name. If `long double` has 53 bits, it uses a 24/53 pair; if it has
113 bits, it uses 53/113. Binary128 uses `long double` only when its width lies
strictly between binary64 and binary128. The effective widths are exported.
`--density-guard-bits` can separately set the MPFR density precision while
preserving wider geometry; its lower comparison is 64 bits below that guard,
with the same minimum of 24.

A failed pair does not by itself prove that its higher-precision value is
inaccurate: the lower comparison may be limiting. Keep the failure and inspect
its reason. Try the available native sequence (`fp64`, `long-double`, `fp128`)
before MPFR, and repeat the complete query. Increasing precision does not repair
a contour-budget failure or geometry that has already lost information.

A different tolerance can be justified by the accuracy actually needed. Keep
the quadrature/domain fixed, compare accepted values against stricter wider
references at representative difficult points, and retain the failed cases.
Changing the tolerance until every point passes supplies no such validation.

## A fast ordinary example and a demanding one

The [four-spinon tutorial](tutorials/xxx-four-spinon.md) measures precision at
fixed momentum, frequency, domain, quadrature and contour controls. Its modest
order-4/T8 example succeeds in x86 `long double` and agrees with MPFR at 128 and
192 bits far beyond the accuracy needed for that demonstration. This validates
that point and those controls. Other momenta/frequencies and finer orders
encounter nodes that the same native precision pair cannot certify. The
tutorial retains those failures and distinguishes the cheap starting command
from the more demanding refinement study. The platform matters: check the
significand bits and exponent range with `bethe-info` before transferring a
`long-double` recommendation.

The [one-magnon example](tutorials/xxx-precision.md) shows actual information
loss. A computed root residual can be zero while subtracting the extensive
energy erases the small energy shift. It compares five arithmetic formats with
an independent exact dispersion. A stable expression for the observable avoids
one subtraction, but cannot repair a rounded root; the root itself must then
be recalculated. For that special one-root problem, the exact sine formula
avoids both difficulties. This is why a stable reformulation should be
considered alongside additional working precision.

## Validated small API examples

The exact six-sector [spectrum](xxx-spinon-sectors.md#six-spinons-inspect-the-sampling-strategy)
and [moment](xxx-spinon-moments.md#refine-quadrature-and-sampling-separately)
commands were repeated with the same seed, sample streams, cutoff and contour
controls in six arithmetic formats. Binary64 fails the arithmetic screen in
both examples. Native 80-bit, binary128 and MPFR128/192/384 complete them.
The native 80-bit estimates, replicate values and maxima, sampling errors and
moment covariance agree with the 384-bit reference to better than
$`4\times10^{-16}`$ relative. Binary128 agrees to better than
$`3\times10^{-30}`$. These checks support the economical native commands for
those specific examples, not arbitrary seeds or extended domains.

| Arithmetic | Spectrum CPU seconds | Moments CPU seconds |
| --- | ---: | ---: |
| Native 80-bit | 0.036 | 0.040 |
| Binary128 | 0.611 | 0.644 |
| MPFR128 | 3.78 | 3.81 |
| MPFR384 | 8.45 | 8.32 |

These are exported application CPU times on an Intel Core i7-12700KF with two
workers. The compact record also retains elapsed time. The very small samples
and rapidity cutoff limit the physical estimates far more than their arithmetic;
these are execution examples, not precision spectral or sum-rule calculations.

The same follow-up repeats all eight four-spinon refinement points in binary128
with the stronger native 64/113-bit pair. Every point passes the unchanged
$`10^{-8}`$ relative density screen and agrees with the independently retained
MPFR192 quadrature within $`6.3\times10^{-31}`$ relative. Total cost is
348 CPU seconds, about seven times less than the recorded MPFR128 grid.
The earlier binary128 failures used a 53-bit lower comparison.

The [compact audit](validation/spinon-native-precision-20261031.json) retains
all 20 queries, including failed values as null, complete replicate/covariance
tables, controls, effective bit counts, source revisions, raw hashes and costs.
Check it or regenerate it from installed executables with:

```sh
python3 scripts/check_xxx_native_density_precision.py --check
# Needs native binary128 and MPFR support; retains expected binary64 failures.
python3 scripts/check_xxx_native_density_precision.py \
  --solver bethe-xxx-structure-factor-thermo --moments bethe-xxx-moments-thermo \
  --raw-dir native-precision-exports
```

Regeneration keeps all requested queries. A failure cannot be replaced by an
estimate over successful samples. Higher precision is still needed where
geometry, cancellation or exponent range exceed the native formats; the next
example illustrates why those requirements must be measured separately.

## Deep-domain geometry needs its own precision

The six-spinon [deep-domain probe](validation/six-spinon-geometry-precision-20261031.json)
shows another reason for extra precision: representing an endpoint. The current
geometry constructs $`p_{\min}=2\arctan(e^{-R})`$ and
$`p_{\max}=\pi-p_{\min}`$. Keeping the second endpoint distinct from pi requires
roughly $`R/\ln 2`$ significand bits, before further accuracy margins. At
$`R=1040`$, 1024 bits fail this check while 1536 and 2048 bits produce matching
replicate means and maxima with the same 192/256-bit density pair. The final
observable does not need hundreds of accurate digits; its current coordinates
need enough bits to retain a very small soft distance. A reformulation that
stores that distance separately could avoid this subtraction, but would also
need checked momentum sums, support decisions and solved roots. Removing the
endpoint guard would not make the lower-precision calculation reliable.

## A reproducible decision procedure

1. Fix the state or integration domain, all non-arithmetic controls and, for
   sampling, the exact seed and stream identities. Save every outcome.
2. Compare the observable at increasing supported precisions. Record the
   effective bits, statuses, diagnostics and CPU/wall time. A failed run is
   evidence; it has no usable partial estimate.
3. Repeat representative difficult points, including edges or near-colliding
   roots. Recompute geometry and roots when checking their precision; changing
   only a density evaluator cannot repair coordinates already rounded away.
4. With an adequate arithmetic choice, independently refine the solver,
   quadrature, contour/series, sampling and domain controls that apply. Preserve
   their distinct meanings when reporting uncertainties.
5. Use the cheapest validated setting for ordinary examples. Keep expensive
   reference runs optional, and repeat the checks when parameters or domains
   change. Displayed digits identify a numerical record; they are not an
   accuracy certificate.

The tutorial [precision audit](validation/tutorial-precision.md) records the
scope and measured comparisons behind the current recommendations.
