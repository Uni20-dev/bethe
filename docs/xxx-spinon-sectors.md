# Selected thermodynamic spinon spectra

[XXX family](xxx-family.md) · [Two-spinon formula](xxx-structure-factor-thermo.md)

`bethe-xxx-structure-factor-thermo` keeps its inexpensive two-spinon default.
Use `--sectors` to select distinct entries from `2,4,6`. Four spinons use
support-adapted deterministic quadrature; six use explicit importance sampling.
The output contains the requested contributions and, for more than one sector,
a row named `sum:2+4`, `sum:2+4+6`, or the corresponding selected sum.
No normalization correction is required. For direct frequency integrals, use
[the separate moments frontend](xxx-spinon-moments.md).

For worked calculations, follow [four-spinon refinement](tutorials/xxx-four-spinon.md),
[six-spinon sampling](tutorials/xxx-six-spinon.md) and the
[combined spectral comparison](tutorials/xxx-spinon-slice.md).

## Start with a small point

```sh
bethe-xxx-structure-factor-thermo --sectors 2,4 \
  --momentum 1.5 --omega 1.8 --precision long-double \
  --four-order 2 --four-cutoff 2 --threads 2 --diagnostics
```

This deliberately small finite-domain calculation checks the workflow. It is
not a converged four-spinon response. To compare quadrature order and domain
size independently, add two controls:

```sh
bethe-xxx-structure-factor-thermo --sectors 2,4 \
  --momentum 1.5 --omega 1.8 --precision fp128 \
  --four-order 4 --four-cutoff 8 \
  --four-refine-order 8 --four-extended-cutoff 12 \
  --density-relative-tolerance 1e-8 --density-absolute-tolerance 1e-30 \
  --tolerance 1e-12 \
  --threads 2 --csv spectrum.csv --csv-table diagnostics=checks.csv
```

This evaluates coarse/fine orders on both domains. The displayed four-spinon
value uses the fine order and extended domain. The order difference and signed
cutoff change are separate fields, not a total error bar or an automatic tail
correction. `--four-max-nodes` limits each geometry control; exceeding it leaves
an explicit failed contribution.

The four component quadratures pass with the native 64/113-bit density pair
on the tested x86 system, and agree with a wider MPFR reference. See the
[measured native precision check](numerical-accuracy.md#validated-small-api-examples).
On a platform without an intermediate native format, inspect the actual pair
and use MPFR if the complete query cannot pass.

The four-spinon cutoff is a **log-radius cutoff** in the
[support-adapted integration coordinates](xxx-four-spinon-integration.md).
The six-spinon cutoff restricts **each absolute rapidity**. Equal numerical
cutoffs do not define equal domains. Both are partial integrations, and their
sum is labelled accordingly. The global-moment common-shift reduction has yet
another domain: individual rapidity gaps with the common shift unrestricted.

## Six-spinons: inspect the sampling strategy

A small fresh execution check is:

```sh
bethe-xxx-structure-factor-thermo --sectors 6 \
  --momentum 1.57079632679489661923132169163975 --omega 2 \
  --precision long-double --six-cutoff 8 --replicates 3 --samples 8 \
  --seed 42 --threads 2 --diagnostics --replicate-output
```

Three replicates of eight configurations are not enough to claim a reliable
spectrum. Inspect support counts, replicate scatter and influential samples,
then increase independent work and vary the proposal and cutoff. A pilot that
misses support returns `no_support_sampled`, with no density, rather than zero.
A failed replicate prevents an estimate over the survivors.

For this exact point, seed and controls, 80-bit `long double` agrees with the
384-bit MPFR reference in the estimates, replicate means, maxima and sampling
error to better than $`4\times10^{-16}`$ relative. This validates the arithmetic
of the example, not its sampling accuracy or omitted tail. Other seeds or
domains can require wider arithmetic; follow the
[precision procedure](numerical-accuracy.md#why-evaluate-a-density-at-two-precisions).

`--samples-per-task` (default 16) sets the maximum observations per Async task.
Smaller blocks can keep more of `--threads` busy when observations have uneven
costs, including within a single replicate. This scheduling option does not
change the stream, replicate count or sample-order accumulation, and is recorded
in output metadata. See [parallel calculations](parallel-calculations.md#six-spinon-monte-carlo).

The physical mixture remains explicit:

| Option | Default | Meaning |
| --- | --- | --- |
| `--six-cutoff` | 8 | Absolute rapidity bound |
| `--six-soft-fraction` | 0.5 | Soft component fraction of the native marginal |
| `--six-background-weight` | 0.5 | Background joint-proposal weight |
| `--six-one-pair-weight` | 0.25 | One-soft-pair joint-proposal weight |
| `--six-two-pairs-weight` | 0.25 | Two-soft-pair joint-proposal weight |
| `--six-one-pair-scale`, `--six-one-pair-gap` | 4, 0.5 | One-pair center and gap scales |
| `--six-two-pairs-scale`, `--six-two-pairs-gap` | 8, 1 | Two-pair center and gap scales |
| `--replicates`, `--samples` | 8, 32 | Independent replicates and configurations per replicate |
| `--seed` | 0 | Reproducible stream identity |
| `--uniform-bits` | 52 | Exact open midpoint grid, separate from arithmetic precision |

Mixture weights must be positive and sum to one. The integrator evaluates the
complete mixture PDF, including overlapping channels; changing mixture weights
does not change the physical normalization. Every point uses distinct replicate
IDs, starting at `grid_index * replicates`. This avoids silently reusing samples
between grid points. The seed and stream rule are recorded in metadata.

See [six-spinon spectral sampling](xxx-six-spinon-structure-factor.md) and
[the proposal geometry](xxx-spinon-sampling.md) for the full physical construction.
The [Monte Carlo engine](monte-carlo.md) manages execution and statistics; it
does not choose or adapt the physical importance distribution.

## Precision and density recovery

By default, `--precision` controls retained geometry, output and the higher density check.
The lower comparison uses binary32 for binary64 and binary64 for binary80.
Binary128 uses native `long double` when its actual significand width lies
strictly between binary64 and binary128; otherwise it uses binary64. MPFR uses
64 fewer bits, with a minimum of 24. Both precisions and their contour controls
are recorded. The [accuracy guide](numerical-accuracy.md#why-evaluate-a-density-at-two-precisions)
explains why a failed pair need not imply an inaccurate higher-precision value.
Raising precision does not replace order, sampling or cutoff studies.

The six-spinon surface keeps endpoint distances and integer multiples of pi
separately when summing free momenta. A cutoff of 1040 therefore does not by
itself require a 1500-bit significand to distinguish a momentum from pi.
Endpoint distances and energies now retain their common exponential factor as
`log_scale = -abs(beta)`. The cutoff distance is compared logarithmically.
Ordinary projections below the square of the arithmetic epsilon carry an
absolute omission bound into the support and conditioning screens; a projected
zero is not a physical zero. This lets binary64 geometry handle cutoffs whose
endpoint distances are outside its exponent range. Unresolved boundaries still
fail explicitly.

Near a degenerate solved pair, arithmetic can still amplify roundoff. The
surface checks conditioning of the solved rapidities and log coarea Jacobian
against the square root of the input epsilon. MPFR retries an unresolved
surface at twice and then four times the requested precision, with results
rounded back to the requested output precision; native formats report the
failure. The input rapidities and query values remain the same represented
numbers. Refinement cannot recover input digits already rounded away. These
checks are arithmetic indicators, not interval certification or an error
bound on the spectral density. Exhausted geometry remains a failed sample.
The internal surface result records working bits and retry count; CLI
precision metadata describes the requested base and output precision.

The [fixed 64-draw validation](validation/six-spinon-endpoint-geometry-20261101.json)
at cutoff 1040 completes at 256 geometry bits with the same sample counts as
the archived 1536-bit calculation. Its geometry-only replay at 128 and 256
bits matches support classifications against 2048 bits without any internal
refinement. This checks those arithmetic paths; the small pilot does not
establish sampling or cutoff convergence.

Before the scaled-endpoint stage, the
[precision-floor probes](validation/six-spinon-precision-floor-20261101.json)
extended this to 1359 geometry solves over seven queries and prescribed deep and
nearly degenerate configurations. The raw density comparison uses 42 labelled
configurations, including rapidities of magnitude 1039. Relative to 512-bit
geometry and 320-bit density, the largest error in a weighted contribution was:

| Arithmetic | Largest relative change | Returned labels |
| --- | ---: | ---: |
| MPFR, 53 significand bits | 2.53e-6 | 42/42 |
| MPFR, 64 significand bits | 1.19e-7 | 42/42 |
| MPFR, 80 significand bits | 1.72e-12 | 42/42 |
| Native extended80, 64 significand bits | 7.07e-13 | 40/42 |
| Native binary128, 113 significand bits | 7.08e-28 | 42/42 |

These raw comparisons recombine returned densities and log weights in
200-digit Decimal arithmetic. MPFR uses bounded internal geometry retries;
extended80 rejects the nearly degenerate geometry. That case regenerates
the energy at each precision, so its error includes input rounding. The geometry
then required an explicitly representable cutoff endpoint distance, excluding
native binary64. MPFR with 53 bits has a much wider exponent range.

The complete 64-draw pilot also passes the existing recovery controller with
96-bit geometry and an explicitly configured 80/96-bit density pair. Its mean
changes by 9.32e-26 relative to the reference; a 64/80-bit pair rejects two draws.
This library policy differs from the CLI's default 64-bit precision gap.
The raw density alone reaches 1e-6 accuracy at 40 MPFR bits in this finite probe,
but that leaves almost no margin and does not satisfy the recovery screens.

Extreme mixed-sign configurations still fail low-precision error screens even
with larger contour budgets. Three representative stress configurations pass
an offline replay of the screens at 320/384 bits with controls (3120,128) and
(3136,132); their values change by at most 1.71e-23 from the original controls.
This supports a low-precision first attempt with fallback, not relaxing the
screens or claiming a universal precision floor. No production defaults change.

Native binary64 returns 40 of 84 raw density evaluations (two controls for each
of 42 labelled configurations). The largest observed relative difference among
those successful evaluations is about 1.52e-12. The remaining evaluations fail
at the Gamma-residue seed because of exponent range. In the original pilot, four of 20 labels fail this way and
carry 21.6% of its reference weighted sum, so discarding failures is not viable.
For the same 84-evaluation raw density workload, native extended80 takes 0.59
CPU seconds, binary128 14.18 seconds and 256-bit MPFR 56.23 seconds on the
recorded host. These single-run timings exclude recovery checks and geometry.

The [scaled-endpoint replay](validation/scaled-spinon-endpoints-20261101.json)
is the first stage toward wider native-binary64 coverage. It repeats all 1359
geometry solves against the frozen 512-bit reference. Native binary64 resolves
1358 cases with no support mismatches: 52 inside and 1306 outside. The maximum
absolute log-weight difference is 1.35e-11. The nearly degenerate case remains
an explicit conditioning failure; native extended80 rejects the same case.
The deep reference fixtures, bounded-projection checks and unresolved-boundary
tests exercise this behavior independently of the replay.

The next stage is to retain residue amplitudes and density prefactors in scaled
form. Recovery and weighted observations must then compare errors and combine
scales before materializing the final contribution. This should preserve
conditioning checks and higher-precision fallback. The geometry stage alone
does not make the deep-domain density representable in binary64, change the
recovery policy, or establish that pure binary64 covers every configuration.

Build `bethe_xxx_spinon_precision_probe` explicitly to run these opt-in probes.
Its `geometry` mode replays the fixed designs; `density` reads lines containing
an identifier, log weight and six rapidities; `checked bits-96 80` runs the
complete pilot through recovery. The store retains the raw inputs, runners,
outputs and receipts. Regenerate the compact record using
`scripts/check_xxx_six_spinon_precision_probe.py STORE_DIRECTORY --pilot-reference FILE`.
`FILE` is the archived `geometry-1536.json` pilot output (1536-bit geometry,
192/256-bit density). It is a separate required input; its SHA256 appears in
the summary. Neither its path nor its parent directory is inferred from
`STORE_DIRECTORY`. To check relocation using temporary copies in host-local
storage, run
`scripts/check_xxx_spinon_archive_relocation.py STORE_DIRECTORY --pilot-reference FILE --scratch-directory SCRATCH_DIRECTORY`.
Choose an existing host-local directory for the temporary copies.
The later geometry-only record uses
`scripts/check_xxx_scaled_spinon_endpoints.py STORE_DIRECTORY --reference-directory DIRECTORY`.
The reference directory must contain the frozen PR #93 `geometry-bits-512.jsonl`
and `geometry-receipt.json`; their hash relationship is checked before analysis.
It need not be adjacent to the stage-1 archive. Its relocation check is
`scripts/check_xxx_spinon_archive_relocation.py STORE_DIRECTORY --reference-directory DIRECTORY --scratch-directory SCRATCH_DIRECTORY`;
this also verifies rejection of an altered reference file.

For MPFR, `--density-guard-bits B` lets the W4/W6 density use a smaller guard
precision while retaining the geometry and output precision. It requires sector
4 or 6 and $`25\le B\le N`$ for `--precision bits-N`; native formats reject this
override. The lower comparison uses `max(24, B-64)` bits. Both arithmetic checks
and both contour budgets still have to pass. For example, a query may use
768-bit geometry while its density passes at 192/256 bits:
use `--precision bits-768 --density-guard-bits 256`. Validate this separation
against a higher density precision on the same complete query before using it
for a larger calculation. The two-spinon kernel has its own `--tolerance`.

The density controller checks two increasing contour budgets at two distinct
precisions. The defaults are 64/80 residue terms and tail orders 16/20, with
relative tolerance `1e-6`, absolute tolerance zero and a 60-second budget per
configuration. Use `--density-terms`, `--density-tail-order`,
`--density-relative-tolerance`, `--density-absolute-tolerance` and
`--density-seconds` to set these controls. The second budget uses 16 additional
terms and four additional tail orders. An exhausted screen leaves a failure;
it does not redraw a difficult sample.

The two-spinon contribution retains its separate `--tolerance`,
`--max-evaluations` and `--max-quadrature-order` controls. Its reported kernel
error is not an error estimate for a higher-sector partial sum.

## Retrying six-spinon observations at wider precision

`--six-precision-stages fp80,fp128,reference` tries each complete six-spinon
observation in increasing arithmetic precision. It requires sector 6 and an
MPFR reference selected with `--precision bits-N` or `digits-N`. Omit an
unavailable native stage explicitly; an unsupported requested format is an
error. The list must be ordered, contain no duplicates and end in `reference`.
Without this option, the fixed-precision calculation is unchanged.

For example, this small interface check requests eight draws, with the same
default finite domain and density controls at every stage:

```sh
bethe-xxx-structure-factor-thermo --sectors 6 --momentum 1.5 --omega 2 \
  --precision bits-192 --six-precision-stages fp80,fp128,reference \
  --replicates 2 --samples 4 --seed 42 --threads 2 \
  --diagnostics --replicate-output --quiet --json staged-point.json
```

Eight draws do not resolve the physical spectrum. Increase the sampling and
domain controls independently; the arithmetic stages do not remove their
uncertainties. Validate a proposed policy against a wider calculation with
the same sample identities before using it for reference data.

Each retry starts from the **same exact uniform sample**, rebuilding the
proposal, energy surface, density and weighted observation at the next stage.
Geometry conditioning, exponent-range failures and exhausted density screens
can trigger a retry. Invalid inputs, inconsistent controls, cancellation and
time limits remain terminal failures. No observation is discarded or redrawn.
The reference point determines continuum boundaries before native evaluation.

Every native stage uses its own paired density precisions and retains the
requested tolerances and contour controls. `--density-guard-bits` applies to
the MPFR reference; its guard must be at least as wide as every native stage.
`--precision` sets the reference geometry, accumulation and output precision,
and also the arithmetic of selected two/four-spinon sectors. Metadata records
these separately from each six-spinon stage. A 192-bit output scalar therefore
does not claim that all observations were evaluated to 192-bit accuracy.

## Tables and failures

The `spectrum` table contains `q`, `omega`, `sector`, `density`,
`estimated_error`, `evaluations` and `status`. `estimated_error` applies only
to the two-spinon kernel; higher-sector and sum rows leave it empty.
Evaluation counts mean transition-integrand evaluations for sector 2 and
checked-density attempts for sectors 4/6; the sum row uses zero rather than
combining unlike work units.

The `continuum` table gives unrestricted kinematic edges for each selected
sector. It does not describe cutoff-induced holes or certify sampled support.

Optional detail tables are available:

- `diagnostics`: two-spinon kernel error; four-spinon order and cutoff changes;
  contour and arithmetic diagnostics; six-spinon sampling error; work count
  and reason. Available when sectors include 4 or 6.
- `replicates`: six-spinon replicate IDs, completed sample and supported-root counts,
  density-attempt counts, largest absolute observation, density and
  numerical diagnostics, failed sample identity and reason. Available with 6.
- `recovery`: stage attempts, accepted observations, retries, inclusive supported
  roots and density attempts, plus failure buckets with status and reason.
  Available with `--six-precision-stages`; unused stages have explicit zero rows.

`--diagnostics` shows diagnostics; `--replicate-output` shows replicates and,
with staged precision, recovery. These selections also include the tables in
a JSON document. A named CSV/TSV export, such as `--csv-table recovery=retry.csv`,
requests that table without showing it on screen. Unknown or unavailable table
names fail before exports are opened. Complete recovery-aware analysis needs
both detail options, as in the example above.

Physical replicate work counts only the final observation. Recovery work counts
all retained attempts, including failed stages and a terminal failed sample;
it excludes speculative work discarded after ordered failure handling. Rows
are keyed by momentum, frequency, replicate and stage. `kind=stage` holds totals;
`kind=failure` groups identical failure tags and reasons. Failure buckets have
absent numeric `accepted`/`retried` fields. Empty enum strings mean no tag;
the literal `none` is a distinct enum value. Failed runs keep these records and
an absent spectral estimate, rather than reporting a mean of survivors.

The primary spectrum streams during calculation. Detail rows are retained only
when requested, then written after the primary table. Their explicit combined
buffer limit is `--detail-row-limit`, default 100000 rows; it applies even with
`--no-retain`, which disables each output table's own row retention. Split a
large grid or raise this limit when requesting detailed records. A spectrum-only
stream does not build those detail buffers.
For recovery, the preflight bound allows a distinct failure reason for every
draw at every stage, plus the stage-total rows. It is checked before opening
files, even if the eventual successful run uses far fewer rows.

A failed requested contribution makes its partial-sum density absent, even if
other sectors succeeded. Thresholds also have no finite sum when one term has
no finite point value. The exit code is 2 for incomplete numerical estimates,
1 for invalid input or delivery failures, and 0 for successful estimates and
identified thresholds. Success at finite cutoffs is not a convergence claim.

Channel conversion is internal: `--channel raising` multiplies longitudinal
values and all their diagnostics by two. Exchange conversion uses frequency
omega/J and density 1/J. No broadening, fitted normalization or archive is needed.

`--threads` bounds the shared Uni20 scheduler. Higher-sector configurations and
replicates execute within bounded admission windows. Independent sectors start
before waiting; the final point result joins them and drains started work after
exceptions. Grid points are collected in order. The two-spinon-only path retains
its existing grouped parallel grid for inexpensive evaluations.
