# Resolving the four-spinon contribution

<!-- BEGIN GENERATED PAGE LINK -->
[Read this page on the documentation site](https://uni20-dev.github.io/bethe/docs/tutorials/xxx-four-spinon/)
<!-- END GENERATED PAGE LINK -->

[Two-spinon thermodynamic spectrum](xxx-structure-factor-thermo.md) ·
[Selected-sector options](../xxx-spinon-sectors.md)

The two-spinon response is a useful starting point for an iMPS comparison, but
it does not exhaust the local-spin spectrum. This tutorial adds four spinons
at one interior spectral point and examines two different numerical questions:
how well the quadrature resolves its chosen domain, and how much weight changes
when that domain is extended.

We use the infinite spin-$`1/2`$ XXX chain at zero temperature and field,
$`H=J\sum_j\mathbf S_j\cdot\mathbf S_{j+1}`$, with $`J=1`$. The point is
$`q=1.5`$ radians/site and $`\omega=1.8J`$. Output is the unbroadened longitudinal
density $`S^{zz}(q,\omega)`$; frequency integrals use $`d\omega/(2\pi)`$.

## Calculate both contributions

With an MPFR-enabled [build or installation](../building.md), run:

```sh
bethe-xxx-structure-factor-thermo --sectors 2,4 \
  --momentum 1.5 --omega 1.8 --precision bits-192 \
  --four-order 4 --four-cutoff 8 \
  --density-relative-tolerance 1e-8 --density-absolute-tolerance 1e-30 \
  --threads 2 --csv point.csv --csv-table diagnostics=point-checks.csv
```

From a build directory, prefix the program with `./`. The three output rows
are the two-spinon density, the four-spinon estimate and `sum:2+4`. They already
have compatible physical normalization: add the contributions directly.
The sum is still a partial response; six and higher sectors are absent.

This deliberately modest rule gives approximately

| Contribution | Density at this point |
| --- | ---: |
| Two spinons | 0.8731304551 |
| Four spinons, order 4 and cutoff 8 | 0.1869408679 |
| Their partial sum | 1.0600713230 |

These extra digits identify the calculation, not its accuracy as an
unrestricted four-spinon response. The refinement study below shows much
larger integration changes than the arithmetic diagnostics.

## Change order and soft cutoff separately

![Four-spinon order/cutoff study and measured CPU cost](figures/xxx-four-spinon-convergence.svg#gh-light-mode-only)
![Four-spinon order/cutoff study and measured CPU cost](figures/xxx-four-spinon-convergence-dark.svg#gh-dark-mode-only)

Each curve holds the soft cutoff fixed and increases the Gauss order per
panel. Comparing points vertically at the same order changes the integration
domain instead. The two-spinon baseline is unchanged throughout this study.

`--four-cutoff T` is a **log-radius cutoff** in the support-adapted four-spinon
coordinates. Increasing it includes configurations closer to a soft endpoint.
It is not a cutoff on frequency, a number of spinons, or a bound on every
absolute rapidity. See the [integration geometry](../xxx-four-spinon-integration.md)
for its definition.

The higher-order points can stabilize within one curve while the curves
remain separated. That is convergence within a truncated domain, not evidence
that the soft contribution has been exhausted. At order 16, changing the cutoff
from 8 to 12 changes the four-spinon density by about 0.02647; changing the
order from 8 to 16 at cutoff 12 changes it by about 0.00011. A small order
change alone would therefore give a misleading accuracy impression here.

The frontend can make the four comparisons in one command. Add to the first
calculation:

```sh
--four-refine-order 8 --four-extended-cutoff 12
```

It evaluates coarse/fine orders on both domains. The displayed four-spinon
value uses the fine order on the extended domain. `quadrature_difference`
and `cutoff_change` report separate comparisons; neither is an automatic
estimate or correction of the remaining tail. Check the individual diagnostic
columns as well as the density and its status.

## Arithmetic accuracy is a different check

The measured runs use 192-bit geometry/output arithmetic, with density
comparisons at 128 and 192 bits. The two contour controls use 64/80 terms and
tail orders 16/20. These choices and the acceptance tolerances are retained in
the exported metadata.

At the initial order-4/cutoff-8 point, integrated contour and arithmetic
screening diagnostics are below $`10^{-34}`$. They are far smaller than the
order and cutoff changes. Tightening the density tolerance cannot resolve
weight that lies outside the domain or between insufficient quadrature nodes.
A larger production calculation should also repeat the complete geometry at
higher precision, independently of these density-only checks.

The absolute density tolerance matters near a zero of the density, where a
relative-only check may be ill-conditioned. Failed nodes remain explicit
failures: no incomplete prefix is used to form a successful contribution.
If a requested sector fails, its sum row is `incomplete` rather than a sum
of the surviving sectors.

## Cost and practical scope

The first point used 33 seconds wall time and 66 CPU seconds with two workers
on an Intel Core i7-12700KF. The right-hand panel shows measured CPU cost for
the complete queries, including the inexpensive two-spinon reference. These
measurements come from a shared host, with at most two simultaneous two-worker
queries; they are not portable timing promises. The largest point, order 16
and cutoff 12, used about 11.4 minutes wall time and 22.8 CPU minutes.

Doubling order increases the number of nodes in both integration directions.
Start with a few representative points and separate the order and domain
checks before requesting a dense spectral grid. More scheduler workers help
with independent node evaluations; they do not replace numerical refinement.

For integrated intensity or the first frequency moment, use
[`bethe-xxx-moments-thermo`](../xxx-spinon-moments.md). It evaluates the moment
measure directly rather than integrating a plotted frequency grid. Its
fixed-momentum and global four-spinon domains have different cutoff controls;
copying the number 8 from this example does not copy the physical domain.

## Try a small six-spinon slice

The same program can evaluate six spinons using independent sampling
replicates. This is a small **execution and reporting example**, at a different
momentum from the four-spinon study above:

```sh
bethe-xxx-structure-factor-thermo --sectors 6 \
  --momentum 1.57079632679489661923132169163975 \
  --omega-min 1.8 --omega-max 2.2 --points 3 \
  --precision bits-384 --six-cutoff 8 --replicates 4 --samples 32 \
  --seed 42 --threads 2 --density-relative-tolerance 1e-5 \
  --diagnostics --replicate-output --quiet --json slice.json
```

A fresh run of the installed program gave:

| $`\omega/J`$ | Six-spinon density | Sampling standard error |
| --- | ---: | ---: |
| 1.8 | 0.00003123 | 0.00002380 |
| 2.0 | 0.00018661 | 0.00012265 |
| 2.2 | 0.00017248 | 0.00011482 |

These errors are 66–76% of the estimates. Four replicates do not establish
reliable confidence intervals, and the absolute-rapidity cutoff 8 also omits
substantial soft weight. The table demonstrates fresh calculation and separate
error reporting; it is not a resolved six-spinon curve. More samples address
sampling variation, while a larger cutoff addresses a different omission.
The two cutoff types must not be equated: `--six-cutoff` bounds absolute
rapidities, whereas `--four-cutoff` bounds the soft log-radius coordinate.

The run used 62 seconds wall time and 101 CPU seconds on the same CPU model,
with two workers. Its [compact record](../validation/six-spinon-cli-20261031.json)
contains controls, source revisions and timings; the replicate exports allowed
an independent recomputation of the means and sampling errors. See
[six-spinon sampling](../xxx-six-spinon-structure-factor.md) before extending it.
To include lower sectors at these **same** spectral points, use `--sectors 2,4,6`
and refine their numerical controls separately. No normalization step is needed.

## Reproduce and extend the study

The checked [compact data](data/xxx-four-spinon-convergence.json) retain all
eight finite-domain estimates, separate diagnostics, work counts, timings,
compiler/hardware information and exact Bethe/Uni20 source revisions.
With NumPy and Matplotlib available, reproduce both light and dark figures:

```sh
python3 scripts/plot_xxx_four_spinon_tutorial.py
```

To recalculate the points from fresh inputs, put the program on `PATH`:

```sh
python3 scripts/plot_xxx_four_spinon_tutorial.py \
  --solver bethe-xxx-structure-factor-thermo --threads 2 \
  --raw-dir four-spinon-exports
```

Use a new raw-output directory; existing exports are not overwritten. This
runs all eight calculations and can take tens of minutes. The script validates
controls, statuses, provenance, the two-spinon reference and physical partial
sums before replacing the compact dataset or drawing a figure. No target
four-spinon fraction is inserted. Literature is available through
`--references` and [CITATIONS.md](../../CITATIONS.md).
