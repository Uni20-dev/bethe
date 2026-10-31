# Resolving six spinons with importance sampling

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

## Calculate an interior point

An MPFR-enabled [build or installation](../building.md) supports the explicit
precision used here. This decimal input denotes a point extremely close to
$`q=\pi/2`$; its exact string is retained for comparing runs.

```sh
q=1.5707963267948966192313216916397514420985846996875529104874722961539082031431045
bethe-xxx-structure-factor-thermo --sectors 6 \
  --momentum "$q" --omega 2 --exchange 1 --channel zz \
  --precision bits-512 --six-cutoff 256 \
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

Four controls answer different questions:

| Change | What the comparison tests |
| --- | --- |
| More independent draws/replicates | Sampling scatter and sensitivity to rare large weights |
| Larger `--six-cutoff` | Weight outside the current rapidity domain |
| More residue terms and higher tail order | Evaluation of the auxiliary contour sums |
| Higher `--precision` | Resolution of the energy surface and arithmetic in the density |

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
in increasing order. The packer rejects duplicate frequencies and missing
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

The [six-spinon library interface](../xxx-six-spinon-structure-factor.md)
provides independent geometry and recovery precision when needed. To obtain
zeroth and first moments, use the [direct moment calculation](../xxx-spinon-moments.md)
and its own domain checks; a sparse display grid is not a substitute for a
frequency integral.
