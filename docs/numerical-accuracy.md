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

The paired precision check is deliberately conservative. With binary64 geometry,
the density is compared at 24 and 53 significand bits. With x86 `long double`,
it is compared at 53 and 64 bits; binary128 uses 53 and 113. MPFR normally uses
N-64 and N bits, with a minimum lower precision of 24 bits. Therefore a failed
pair does not by itself prove that the higher-precision value is inaccurate:
the lower comparison may be the limiting calculation. Keep the failure, inspect
its reason, and repeat with a genuinely stronger pair rather than discarding
the difficult nodes or weakening the screen to obtain success.

## A fast ordinary example and a demanding one

The [four-spinon tutorial](tutorials/xxx-four-spinon.md) measures precision at
fixed momentum, frequency, domain, quadrature and contour controls. Its modest
order-4/T8 example succeeds in x86 `long double` and agrees with MPFR at 128 and
192 bits far beyond the accuracy needed for that demonstration. Finer orders
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
