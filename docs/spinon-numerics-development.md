# From a reliable spinon calculation to an efficient one

The spinon code illustrates a common route in numerical software: establish a
reliable reference, measure where it spends time, then change the representation
and execution strategy while checking the same inputs. Some apparent precision
requirements turned out to belong to the implementation rather than the physics.

This page records that development. For current commands, start with the
[six-spinon tutorial](tutorials/xxx-six-spinon.md) or
[four-spinon tutorial](tutorials/xxx-four-spinon.md). Dataset reproduction pages
retain the settings used for particular reference calculations. The measurements
below describe their recorded source revisions and hardware; they are not a
table of current performance guarantees.

## Establish a reference before optimizing

The six-spinon response combines several numerical problems. Four sampled
rapidities determine two remaining momenta through energy and momentum
conservation. The density contains special functions and auxiliary contour
sums. Proposal probabilities and the energy-shell Jacobian convert that density
into a weighted observation. Independent replicate means then give the spectral
estimate and its sampling standard error.

Each layer can fail differently. We kept the random-stream identities, input
parameters and failures, and used precision-generic C++ implementations to
repeat the same calculation at different precisions. A comparison must include
support classifications and unsuccessful observations, not just values returned
by both implementations. Otherwise an optimization can appear accurate by
silently dropping difficult configurations that carry appreciable weight.

Arbitrary precision supplied a useful reference. It did not remove contour
truncation, finite-domain error or sampling error. In particular, close agreement
between two precisions says little about whether a Monte Carlo estimate has
sampled its rare, large contributions adequately.

## Separate the precision of different tasks

The four-spinon quadrature supplied one example. At $`q=\pi/2`$,
$`\omega=2J`$, log-radius cutoff 768 and order 4, retaining 1536-bit geometry
but changing the density comparison from 1472/1536 to 192/256 bits reduced
measured CPU cost by a factor of 4.37. The result changed by less than
$`7.90\times10^{-76}`$ relative. The
[calibration record](validation/four-spinon-slice-calibration-20261031.json)
preserves both calculations. This checks density precision at one quadrature;
it does not determine the minimum geometry precision.

Initially, a large geometry precision also made every form-factor evaluation
expensive. Geometry and density do not necessarily need the same number of
bits. Separating them lets the geometry retain its precision while the density
uses a cheaper comparison pair.

For an $`R=256`$, 2,048-draw example, changing the density pair from 448/512 bits
to 192/256 bits while retaining 512-bit geometry reduced measured CPU time from
1,850 to 1,150 seconds. All replicate work counts agreed, and the mean changed
by less than $`2\times10^{-75}`$ relative. Its sampling error remained 15.4%.
The [comparison record](validation/six-spinon-density-precision-20261031.json)
keeps both costs and the numerical results.

This was a useful optimization, but it did not explain why the geometry needed
so many bits in the first place.

## The lower comparison can limit an accurate result

The density screen compares two independently evaluated precisions. Initially,
binary128 compared against binary64 even where a native intermediate format
was available. Some four-spinon refinement quadratures failed that 53/113-bit
pair while MPFR128, with its 64/128-bit pair, passed. The higher precision alone
therefore did not explain which queries succeeded. The
[refinement record](validation/four-spinon-refinement-precision-20261031.json)
retains those failures and the MPFR controls.

Using an intermediate 64-bit-significand `long double` for binary128's lower
comparison let all eight quadratures pass at the unchanged tolerance. Their
S4 values agreed with MPFR192 within $`6.3\times10^{-31}`$ relative and used
348 CPU seconds, compared with about 41 CPU minutes for the recorded MPFR128
grid. The [native comparison](validation/spinon-native-precision-20261031.json)
records the controls and source identities. Those are costs of recorded runs,
not a controlled comparison of every aspect of the two executables.

This does not justify weakening a screen until a query passes. The lower
evaluation is an independent arithmetic check; a failed pair does not identify
which value lost accuracy. Same-input wider references distinguish that case
from a bad higher-precision result. The
[neighbourhood study](validation/four-spinon-neighbourhood-precision-20261031.json)
and [64/113-bit comparison](validation/four-spinon-native-neighbourhood-20261031.json)
retain failures, tolerance studies and wider references across momenta and
frequencies. A successful ordinary point is not a uniform precision guarantee.

## A representation can create its own precision requirement

Very soft spinons have momenta exponentially close to zero or pi. At rapidity
cutoff $`R`$, an endpoint distance is

```math
d_R=2\arctan(e^{-R}).
```

The original surface formed $`\pi-d_R`$. Distinguishing this number from pi
requires approximately $`R/\ln 2`$ significand bits, before allowing for later
calculations. At $`R=1040`$, an archived 64-draw check failed with 256-, 512-
and 1024-bit geometry; 1536 and 2048 bits completed and agreed. The
[original geometry record](validation/six-spinon-geometry-precision-20261031.json)
documents that genuine limitation of the original representation.

Keeping the endpoint distance separate from integer multiples of pi removed
the destructive subtraction. The same 64 draws then completed with 256-bit
geometry and the same 192/256-bit density pair. The mean changed by less than
$`8\times10^{-74}`$ relative to the 1536-bit reference, with the same sample
identities and supported-root counts. Its roughly 100% relative sampling error
was unchanged. See the [endpoint-preserving replay](validation/six-spinon-endpoint-geometry-20261101.json).

The lesson is to inspect what the bits represent before concluding that the
problem requires arbitrary precision throughout. A well-conditioned physical
quantity can be assembled from poorly represented intermediate quantities.
Near-degenerate solved pairs still need conditioning checks: changing the
representation does not remove genuine sensitivity to inputs.

## Significand width and exponent range are different limits

Preserving a small endpoint distance avoids cancellation, but the distance
itself can underflow in native double. A scaled representation stores a bounded
factor and a logarithmic scale separately. Endpoint and energy contributions
can then be compared at compatible scales, with explicit bounds for omitted
tiny terms. A projected numerical zero must not become a physical zero.

The same issue appears inside the density kernel. In a fixed workload of
84 raw evaluations, native double returned 40 values and failed 44 at a
Gamma-residue seed whose magnitude exceeded its exponent range. Four failed
pilot configurations carried 21.6% of that pilot's reference weighted sum.
Discarding those failures would have biased the calculation substantially.

Native extended80 returned all 84 densities, with maximum relative change
$`7.08\times10^{-13}`$ against the wider reference. Binary128 returned all 84
with maximum relative change $`7.08\times10^{-28}`$. The recorded CPU times were:

| Density arithmetic | CPU seconds |
| --- | ---: |
| MPFR, 256 significand bits | 56.23 |
| Native extended80, 64 significand bits | 0.59 |
| Binary128, 113 significand bits | 14.18 |

These are raw kernel timings on an i7-12700KF with GCC 13.3, excluding geometry,
recovery checks and integration. The roughly 95-fold MPFR-to-extended80 ratio
is not an end-to-end speedup. Some accurate raw values still fail the requested
arithmetic screens; checking a higher-precision result on every observation
also costs time. The [precision study](validation/six-spinon-precision-floor-20261101.json)
records failures and screens alongside values and timings. A useful native-first
policy must measure the cost and frequency of its higher-precision fallback.

## Let failures identify the calculation that needs strengthening

Contour truncation and roundoff require different remedies. At $`\omega=2J`$
and $`R=1040`$, wider density arithmetic did not repair two contour failures in
a production sample. Increasing residue counts and tail orders did. Exact
replays and four constructed boundary cases checked the repair without changing
the failed configurations. The [recovery record](validation/six-spinon-contour-recovery-20261031.json)
retains the unsuccessful and successful calculations.

Conversely, a configuration at $`\omega=3J`$ passed its contour check but failed
the 192/256-bit arithmetic comparison. The 256/320-bit pair passed at the same
contour controls and tolerance, and a 320/384-bit comparison agreed within
$`1.3\times10^{-91}`$ in the checked densities. These examples are in the
[slice calibration record](validation/six-spinon-slice-calibration-20261031.json).

Small-domain seed scans made the same point at lower cost. Extended80 completed
38 of 40 three-frequency queries with the tested controls; binary128 completed
the two failed seeds without relaxing the tolerance, and agreed with MPFR
controls. The [native comparison](validation/six-spinon-native-comparison-20261031.json)
and [broader scan](validation/six-spinon-native-seed-scan-20261031.json) retain
the input identities and outcomes. These counts describe those scans, not a
universal failure probability.

The compact native-comparison record can be checked without rerunning the
solver, from the source directory:

```sh
python3 scripts/check_xxx_six_spinon_native_precision.py
```

Increasing the sample count may expose a difficult configuration absent from
the calibration sample. Recovery must retain that observation and recompute it,
or report an incomplete calculation. Neither a different seed nor an average
over surviving observations repairs the original estimate.

## Parallelize work more finely than the statistical replicates

A replicate is a useful unit of statistical analysis, but its runtime varies
with the sampled configurations. Scheduling each complete replicate as one
task left workers idle while the last few slow replicates finished. The first
successful large query used 38.0 CPU hours in 143 minutes with 24 workers,
averaging about sixteen occupied cores. That measurement alone does not assign
every idle interval to scheduling; profiling identified the uneven completion
of whole-replicate tasks as a source of underutilization.

The calculation now splits replicate streams into sample blocks and submits
them through Uni20's asynchronous scheduler. Sample identities and their final
accumulation order remain fixed, independently of task completion order.
This exposes more parallel work without changing the statistical meaning of a
replicate. Very small tasks would instead spend too much time on scheduling
and allocation, so the block size remains an explicit tuning parameter.
The [parallel-calculation guide](parallel-calculations.md#six-spinon-monte-carlo)
describes the interface and cancellation behavior. Whole-query throughput is
the useful performance measure, rather than a brief snapshot of full CPU use.

## Preserve the evidence as the code improves

For every change, retain enough information to compare the same calculation:
input identities, numerical controls, source revisions, failure reasons and
measured costs. Compact records can carry estimates, replicate summaries and
hashes of the original exports; users need not download per-observation
archives to redraw a figure.

The high-precision calculations remain valuable references even when they are
no longer economical production settings. The tutorials should teach the
current supported method, while the records and this development account
explain how it was established. The next useful benchmark is the complete
checked native calculation, including its fallback work: that is what determines
how much finer a spectrum or how many more samples the same compute budget buys.
