# Resolving six spinons with importance sampling

<!-- BEGIN GENERATED PAGE LINK -->
[Read this page on the documentation site](https://uni20-dev.github.io/bethe/docs/tutorials/xxx-six-spinon/)
<!-- END GENERATED PAGE LINK -->

[Two-spinon spectrum](xxx-structure-factor-thermo.md) ·
[Four-spinon quadrature](xxx-four-spinon.md) ·
[Selected-sector options](../xxx-spinon-sectors.md)

Six spinons contribute both inside and above the two-spinon band. Their large
kinematic bandwidth does not mean that most of their intensity lies at high
energy. To compare their response with an iMPS calculation, we need the
spectral density, including contributions from configurations containing very
soft spinons.

We use the infinite spin-$`1/2`$ antiferromagnetic XXX chain at zero temperature
and field, $`H=J\sum_j\mathbf S_j\cdot\mathbf S_{j+1}`$, with $`J=1`$.
The output is the unbroadened longitudinal density $`S^{zz}_6(q,\omega)`$.
Frequency integrals use $`d\omega/(2\pi)`$. The physical normalization is
already included; no multiplier is fitted to a sum rule or another calculation.

## Why sampling needs a physical proposal

Momentum and energy conservation determine two momenta after four rapidities
are chosen. The integrator includes every momentum branch and both labelled
solutions. Its proposal mixes a background distribution with distributions
that place one or two pairs near soft configurations. The full mixture
probability enters each weight, including overlapping ways to generate it.

This is ordinary Monte Carlo integration with independent replicate streams.
Changing the proposal changes the distribution of work and statistical error,
while preserving the target integral. Concentrating draws in one region is
useful only if all required regions retain support and the complete proposal
probability is used. The [geometry guide](../xxx-spinon-sampling.md) describes
the measure and balancing between different choices of the solved pair.

## Start with an inexpensive calculation

This small calculation checks the workflow before extending the rapidity domain
or requesting a resolved spectrum:

```sh
bethe-xxx-structure-factor-thermo --sectors 6 \
  --momentum 1.57079632679489661923132169163975 \
  --omega-min 1.8 --omega-max 2.2 --points 3 \
  --precision long-double --six-cutoff 8 --replicates 4 --samples 32 \
  --seed 42 --threads 2 --density-relative-tolerance 1e-5 \
  --density-absolute-tolerance 0 --density-terms 64 --density-tail-order 16 \
  --diagnostics --replicate-output --json small-slice.json
```

### Why evaluate the form factor at two precisions?

For each supported spinon configuration, the program evaluates the same
form-factor density twice: once at a lower precision and once at a higher
precision, called the **guard precision**. Both evaluations start from the
same geometry coordinates, converted separately into their respective formats,
and use the same contour controls. The lower result is a fresh calculation,
not a rounded copy of the higher result.

This tests sensitivity to roundoff. The formula contains special functions
and sums whose terms can nearly cancel; choosing a floating-point format alone
does not tell us how many useful digits survive. The difference between the
two results contributes to the arithmetic diagnostic, together with the
kernel's roundoff indicators and a check of its scalar prefactor. The program
uses the **higher-precision value**, not an average of the two.

For example, `long-double` on the tested x86 system compares a 53-bit `double`
evaluation with a 64-bit `long double` evaluation. The table below gives both
precisions for every measured run. With MPFR, the higher density precision
normally follows `--precision bits-N`; `--density-guard-bits` can lower it
independently of geometry precision. The comparison uses 64 fewer bits, with
a minimum of 24, so a 192-bit density guard gives a 128/192-bit pair.

The arithmetic diagnostic plus the estimated contour remainder must fit within
$`\mathrm{atol}+\mathrm{rtol}\,|W_{\mathrm{high}}|`$.
The program also repeats this check with a larger contour budget and checks
agreement between the passing results. Agreement is a numerical stability
check, not a rigorous error bound. Disagreement can originate in the lower
precision even when the higher result is accurate: it means the requested
check has not passed, not that the higher format is necessarily inadequate.
Neither check measures Monte Carlo sampling error or the weight excluded by
the rapidity cutoff.

### Compare cost and arithmetic stability

On an x86 machine with native extended80 `long double`, this complete query
took 0.38 CPU seconds. The same domains, streams and numerical controls gave:

| Arithmetic | Density comparison / guard bits | Outcome | CPU seconds |
| --- | ---: | --- | ---: |
| `fp64` | 24 / 53 | Recovery failed | 0.04 |
| `long-double` | 53 / 64 | Complete | 0.38 |
| `fp128` | 53 / 113 | Complete | 6.92 |
| `bits-96` | 32 / 96 | Recovery failed | 6.41 |
| `bits-128` | 64 / 128 | Complete | 43.59 |
| `bits-192` | 128 / 192 | Complete | 55.73 |
| `bits-384` | 320 / 384 | Complete | 101.14 |

All successful runs retain the same replicate identities and work counts. The
largest relative change in the densities from extended80 to 384-bit MPFR is
$`2.7\times10^{-16}`$. The [comparison record](../validation/six-spinon-precision-20261031.json)
contains the commands, effective precision, diagnostics, sampling errors,
hardware and timings, including an additional 256-bit check. These tiny
arithmetic differences do not measure the accuracy of the physical spectrum:
the observed sampling errors remain 66–76%, and cutoff 8 omits substantial
soft weight.

The two failing settings illustrate why the guard precision alone is not
enough. Recovery compares two arithmetic precisions. Binary64 uses a binary32
comparison, which fails arithmetic checks and encounters a subnormal-density
failure in one replicate. MPFR at 96 bits uses a 32-bit comparison and also
fails arithmetic checks. Neither failed query is replaced by an estimate over
its surviving replicates. Native extended80 is the cheapest checked setting
for this particular query; `bits-128` is a measured MPFR alternative. Check
`bethe-info` for the actual native format on your platform, and repeat the
precision comparison when changing the domain, grid or sample streams.

### Changing the seed can change the precision requirement

More samples explore more configurations, including ones with stronger
cancellation. With the same three frequencies, cutoff and 32 samples per
replicate, 38 of the 40 seeds from 1 through 40 complete in extended80 at
relative tolerance `1e-5`. Seeds 2 and 35 fail an arithmetic screen. This
small, fixed scan illustrates the limitation; it is not a failure probability
for other calculations. The [seed comparison record](../validation/six-spinon-seed-precision-20261031.json)
retains all 40 outcomes, including failures.

Seed 2 provides a reproducible example. Use `--seed 2` and distinct output
filenames. First vary only `--precision`, keeping relative tolerance `1e-5`;
the final two rows separately test `--density-relative-tolerance 1e-3`:

| Arithmetic | Relative tolerance | Complete three-point query | CPU seconds |
| --- | ---: | --- | ---: |
| `long-double` | `1e-5` | Arithmetic screen fails | 0.31 |
| `bits-128` | `1e-5` | Complete | 37.81 |
| `bits-192` | `1e-5` | Complete | 48.54 |
| `long-double` | `1e-3` | Complete | 0.31 |
| `fp64` | `1e-3` | Arithmetic and exponent-range failures | 0.03 |

Increasing MPFR precision from 128 to 192 bits changes the three means by
at most $`4.4\times10^{-37}`$ relative, with matching replicate identities
and work counts. This repairs the failure without changing the tolerance,
domain or sample stream. A failed paired-precision check does not by itself
show that the higher-precision value was inaccurate: it shows that the
requested agreement with its lower-precision comparison was not established.

The `1e-3` native run is a separate tolerance experiment. Its means agree with
the unchanged, strict 192-bit reference within $`3.3\times10^{-16}`$ relative;
the largest relative change among replicate means and maxima is
$`7.1\times10^{-15}`$. Such a tolerance can be economical for this coarse
demonstration after checking it against a stricter reference. These results
do not certify it for a larger sample, domain or resolved spectrum. Binary64
also encounters an exponent-range limit, which loosening the tolerance does
not repair.

When your query fails, preserve that result and rerun the **same** frequencies,
seed and sample count at higher precision. Compare successful wider runs
before increasing the sampling budget. Do not search for a passing seed or
average only successful replicates. Choose an inner tolerance for the accuracy
you need, then validate it against a stricter calculation; a passing screen
alone does not establish the error of the spectral integral.

## Extend the soft domain

An MPFR-enabled [build or installation](../building.md) supports the explicit
precision used here. This decimal input denotes a point extremely close to
$`q=\pi/2`$; its exact string is retained for comparing runs.

```sh
q=1.5707963267948966192313216916397514420985846996875529104874722961539082031431045
bethe-xxx-structure-factor-thermo --sectors 6 \
  --momentum "$q" --omega 2 --exchange 1 --channel zz \
  --precision bits-512 --density-guard-bits 256 --six-cutoff 256 \
  --replicates 16 --samples 128 --seed 7304 --threads 4 \
  --density-terms 256 --density-tail-order 24 \
  --density-relative-tolerance 1e-6 --density-absolute-tolerance 0 \
  --diagnostics --replicate-output --json point.json
```

From a build directory, prefix the executable with `./`. The sample count is
per replicate, so this requests 2,048 draws. These settings define a finite
integration domain and numerical checks; they do not promise a particular
relative sampling error. The point is a starting calculation to inspect before
running a frequency grid.

Here geometry and output retain 512 bits, while the density comparison uses
192/256 bits. Replaying all 2,048 draws against the original 448/512-bit density
checks preserves every replicate's work counts and changes the mean by less
than $`2\times10^{-75}`$ relative. In these runs, CPU cost fell from 1,850 to
1,150 seconds (about five minutes wall time with four workers). The
[comparison record](../validation/six-spinon-density-precision-20261031.json)
retains the actual point and its replicate summaries. Its sampling error is
still 15.4%; lowering arithmetic cost does not improve that error or recover
the omitted rapidity tail.

The `spectrum` table contains the density. The `diagnostics` table separates
sampling error from contour and arithmetic diagnostics. The `replicates` table
retains each mean, sample count, supported-root count and largest observation.
An unsuccessful replicate prevents an estimate over the remaining replicates.

The accompanying analysis script independently checks all replicate identities
and recomputes the reported mean and sampling standard error:

```sh
python3 scripts/plot_xxx_six_spinon_tutorial.py \
  --raw point.json --momentum "$q" --output point-check.json
```

Run the script from the source directory. It needs only Python's standard
library for this check and refuses to overwrite an existing analysis file.
It also reports the fraction of the accumulated weight supplied by the
largest observation and by the largest replicate.

## Read the errors before increasing the grid

For $`B`$ complete, equally sized replicate means $`\bar S_b`$, the reported
sampling error is

```math
\bar S=\frac1B\sum_b\bar S_b,\qquad
\mathrm{SE}(\bar S)^2=\frac{1}{B(B-1)}
\sum_b(\bar S_b-\bar S)^2.
```

This estimate measures the scatter that the run observed. A few large weights
can dominate a small run, and an unsampled rare region does not appear in its
error bar. Examine replicate scatter and influential observations, increase
independent work, and check a different seed. Calibrate proposal choices on
separate streams before using them for a validation calculation.

These controls answer different questions:

| Change | What the comparison tests |
| --- | --- |
| More independent draws/replicates | Sampling scatter and sensitivity to rare large weights |
| Larger `--six-cutoff` | Weight outside the current rapidity domain |
| More residue terms and higher tail order | Evaluation of the auxiliary contour sums |
| Higher `--precision` | Energy-surface arithmetic; also density if no guard override is set |
| Higher `--density-guard-bits` at fixed geometry | Density arithmetic alone, with MPFR |

The density check uses $`\mathrm{atol}+\mathrm{rtol}\,|W_6|`$. The relative
term sets the scale for nonzero matrix elements; the absolute term supplies a
floor near zero. For example, an absolute tolerance of $`10^{-30}`$ would not
request thirty-digit accuracy everywhere. These examples use a zero absolute
floor and explicit relative tolerances. Passing this inner check does not
bound sampling, integration or omitted-domain error.

The rapidity cutoff $`R`$ restricts all six rapidities to $`|\beta_j|\le R`$,
where $`p_j=2\arctan(e^{\beta_j})`$. Large $`|\beta_j|`$ means a spinon near
momentum zero or $`\pi`$, with very little energy. Increasing $`R`$ can therefore
add weight inside the bright continuum; it does not merely extend the
high-frequency tail. A small sampling error does not bound this omitted weight.

Compare cutoffs with common precision and sufficiently strong contour controls.
When sample streams are shared, use their paired replicate differences rather
than assuming independent errors for the difference. Keep a failed calculation
as a failed record. A numerical failure or missed sampled support is not a
physical zero.

For two exports differing only in cutoff, the script checks their controls and
computes the paired difference and its standard error:

```sh
python3 scripts/plot_xxx_six_spinon_tutorial.py \
  --raw lower-cutoff.json --compare-raw higher-cutoff.json \
  --momentum "$q" --output cutoff-check.json
```

A small difference does not bound the weight beyond the larger cutoff. Both
estimates can omit the same soft region.

## Distinguish contour failures from arithmetic failures

Larger domains can require both stronger contours and wider arithmetic. The
[slice calibration record](../validation/six-spinon-slice-calibration-20261031.json)
keeps fifteen bounded queries, their commands and costs, two failed queries,
paired cutoff comparisons and exact replays of the problematic configurations.

At $`\omega=2J`$ and $`R=1040`$, a complete 2048-draw query fails with
1040/1056 residue terms and tail orders 32/36. Raising the density guard from
256 to 384 bits leaves the original contour failures. Two stronger contour
choices pass the failed configurations and agree within
$`5\times10^{-61}`$ relative. The cheaper passing choice uses 1560/1576 terms
and tail orders 48/52; a complete replay then passes.

At $`\omega=3J`$, those contours pass their truncation check, but one sampled
configuration fails the 192/256-bit arithmetic comparison: its relative
indicator is $`2.73\times10^{-6}`$ against a $`10^{-6}`$ target. Keeping the
contours and tolerance fixed, the 256/320-bit pair passes; a 320/384-bit replay
agrees within $`1.3\times10^{-91}`$ in the checked densities. The complete
1024-draw query then passes at 256/320 bits, as does the $`\omega=4J`$ anchor.
These cases motivate the wider density pair for the large-domain slice.

Neither repair changes the sample stream or replaces a failed result with an
average of surviving observations. The checks establish usable numerical
controls at the tested points. They do not establish sampling convergence or
bound the omitted soft domain; larger runs still check every observation.

Increasing the budget to 32 replicates of 768 draws at $`\omega=2J`$ exposed
two further contour failures. The complete query was rejected: its other
replicates do not supply a usable spectral estimate. Replaying the two failed
configurations showed that wider arithmetic again left the contour problem
unresolved. Stronger contour controls repaired both cases without changing the
precision, tolerance or sample identities.

We also tested four constructed geometries close to the rapidity cutoff;
these are additional stress tests, not Monte Carlo samples. The
[recovery record](../validation/six-spinon-contour-recovery-20261031.json)
retains 22 bounded checks. At 3120/3136 residue terms and tail orders 52/56,
both failed configurations and all four constructed cases pass. Their largest
relative contour diagnostic is below $`8.6\times10^{-11}`$, compared with the
unchanged $`10^{-6}`$ target. Separate increases in arithmetic precision and
residue count preserve the checked densities. These are the controls in the
[large-domain reproduction recipe](xxx-six-spinon-reproduction.md); the
replacement whole-query calculation is still in progress. Passing these
selected configurations does not replace checking every production observation.

## When the geometry needs higher precision

Extending the rapidity domain can demand much more precision even with a loose
sampling target. In a separate fixed-stream check at $`R=1040`$, we held the
density comparison at 192/256 bits and varied only the geometry precision.
At 256, 512 and 1024 bits the program rejected the surface with
`surface cutoff unresolved; increase geometry precision`. At 1536 and 2048
bits, all 64 draws completed with matching supported-root counts and agreeing
replicate estimates. The [geometry comparison record](../validation/six-spinon-geometry-precision-20261031.json)
retains these controls, outcomes and costs.

Very soft momenta approach zero or $`\pi`$ exponentially in the rapidity.
The current surface representation needs enough bits to distinguish its domain
boundaries: $`\pi-2\arctan(e^{-R})`$ must remain distinguishable from $`\pi`$.
This alone requires roughly $`R/\ln 2`$ bits, before allowing for accuracy in
the subsequent calculation. Raising density precision cannot repair an unresolved surface:
rebuild the geometry with higher `--precision`. With MPFR,
`--density-guard-bits 256` keeps the density checks at 192/256 bits while the
geometry precision increases. That separation substantially reduces cost.
The 64-draw check establishes feasibility and arithmetic agreement only; it
does not resolve the spectrum or its omitted soft tail.

For a new calculation, hold the domain, contour budgets and sample identities
fixed while comparing precisions. Then vary the non-arithmetic controls
separately. A contour-tail failure calls for a contour study; a large sampling
error calls for more independent samples and examination of influential weights.
Keep failed configurations visible, and retain the diagnostic that identifies
what needs improvement.

## Keep the plotted data independently checkable

The compact `bethe-six-spinon-slice-v1` format retains the full requested
frequency grid and stores numerical metadata once per run, together with
spectrum, diagnostic and replicate rows. It includes each raw export's SHA256,
source revisions, cutoff, proposal controls and timing summary. Local command
paths are omitted. Raw observations and archives are not needed to redraw a
figure, and are not included in the compact file.

For the single-point example above:

```sh
python3 scripts/xxx_six_spinon_tutorial_data.py \
  --raw point.json --momentum "$q" --cutoff 256 \
  --frequency 2 --output six-compact.json
```

Repeat `--raw` for separate runs and `--frequency` for every requested node,
in increasing order. Decimal input frequencies are matched uniquely within the
export's recorded arithmetic-roundoff allowance. The reader returns the input
as `omega` and retains the rounded export coordinate as `computed_omega`; the
stored source rows remain unchanged. There is no interpolation, and ambiguous
matches are errors. The packer rejects duplicate frequencies and missing
nodes. An unavailable node needs an explicit record with `omega`, a `failed`
or `unavailable` status, a reason, and `raw_sha256` (or `null` if no export
exists); supply that list using `--unavailable records.json`.

Both tutorial plotters can use
`read_compact_slice(document, momentum=q, exchange='1')` from the same module.
It recomputes the estimates through the sampling validator and returns ordered
points with their original controls. Unavailable points have `s6: null` and
`sampling_se: null`; they must appear as gaps, not zeroes or connecting lines.
These checks preserve the evidence for an estimate; they do not certify its
sampling or cutoff convergence.

To plot a compact slice, install the optional plotting dependencies and run:

```sh
python3 -m pip install -r docs/tutorials/requirements.txt
python3 scripts/plot_xxx_six_spinon_tutorial.py \
  --compact six-compact.json --momentum "$q" --figure six-slice.svg
```

This writes `six-slice.svg` and `six-slice-dark.svg`, each with linear and
logarithmic density panels. The points are unbroadened; joining lines only
guide the eye. Error bars show the sampling SE, excluding cutoff uncertainty.
Unavailable frequencies interrupt the curve. The two-, four- and six-spinon
upper boundaries are marked when they fall within the requested frequency
range. A single-point file is useful for checking the renderer, but a physical
line shape requires a resolved frequency grid.

The optional [large-domain reproduction recipe](xxx-six-spinon-reproduction.md)
records the installed-program commands and predetermined streams for a
sixteen-frequency slice. Its numerical production is still in progress; the
recipe is not a claim of completed or converged data.

The [six-spinon library interface](../xxx-six-spinon-structure-factor.md)
provides independent geometry and recovery precision when needed. To obtain
zeroth and first moments, use the [direct moment calculation](../xxx-spinon-moments.md)
and its own domain checks; a sparse display grid is not a substitute for a
frequency integral.
