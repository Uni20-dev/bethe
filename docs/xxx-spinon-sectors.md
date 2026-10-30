# Selected thermodynamic spinon spectra

[XXX family](xxx-family.md) · [Two-spinon formula](xxx-structure-factor-thermo.md)

`bethe-xxx-structure-factor-thermo` keeps its inexpensive two-spinon default.
Use `--sectors` to select distinct entries from `2,4,6`. Four spinons use
support-adapted deterministic quadrature; six use explicit importance sampling.
The output contains the requested contributions and, for more than one sector,
a row named `sum:2+4`, `sum:2+4+6`, or the corresponding selected sum.
No normalization correction is required. For direct frequency integrals, use
[the separate moments frontend](xxx-spinon-moments.md).

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
  --momentum 1.5 --omega 1.8 --precision bits-192 \
  --four-order 4 --four-cutoff 8 \
  --four-refine-order 8 --four-extended-cutoff 12 \
  --threads 4 --csv spectrum.csv --csv-table diagnostics=checks.csv
```

This evaluates coarse/fine orders on both domains. The displayed four-spinon
value uses the fine order and extended domain. The order difference and signed
cutoff change are separate fields, not a total error bar or an automatic tail
correction. `--four-max-nodes` limits each geometry control; exceeding it leaves
an explicit failed contribution.

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
  --precision bits-384 --six-cutoff 8 --replicates 3 --samples 8 \
  --seed 42 --threads 2 --diagnostics --replicate-output
```

Three replicates of eight configurations are not enough to claim a reliable
spectrum. Inspect support counts, replicate scatter and influential samples,
then increase independent work and vary the proposal and cutoff. A pilot that
misses support returns `no_support_sampled`, with no density, rather than zero.
A failed replicate prevents an estimate over the survivors.

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

`--precision` controls retained geometry, output and the higher density check.
The lower comparison uses binary32 for binary64, binary64 for binary80 or
binary128, and 64 fewer bits for MPFR, with a minimum of 24. Both precisions
and their contour controls are recorded. Low native precision can fail the
recovery screen, especially for six spinons; choose wider precision in that
case. Raising precision does not replace order, sampling or cutoff studies.

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

## Tables and failures

The `spectrum` table contains `q`, `omega`, `sector`, `density`,
`estimated_error`, `evaluations` and `status`. `estimated_error` applies only
to the two-spinon kernel; higher-sector and sum rows leave it empty.
Evaluation counts mean transition-integrand evaluations for sector 2 and
checked-density attempts for sectors 4/6; the sum row uses zero rather than
combining unlike work units.

The `continuum` table gives unrestricted kinematic edges for each selected
sector. It does not describe cutoff-induced holes or certify sampled support.

Two optional detail tables are available:

- `diagnostics`: two-spinon kernel error; four-spinon order and cutoff changes;
  contour and arithmetic diagnostics; six-spinon sampling error; work count
  and reason. Available when sectors include 4 or 6.
- `replicates`: six-spinon replicate IDs, completed sample and supported-root counts,
  density-attempt counts, largest absolute observation, density and
  numerical diagnostics, failed sample identity and reason. Available with 6.

`--diagnostics` and `--replicate-output` show these on screen and in a JSON
document. A named CSV/TSV export requests that table without showing it on
screen. Unknown or unavailable table names fail before exports are opened.

The primary spectrum streams during calculation. Detail rows are retained only
when requested, then written after the primary table. Their explicit combined
buffer limit is `--detail-row-limit`, default 100000 rows; it applies even with
`--no-retain`, which disables each output table's own row retention. Split a
large grid or raise this limit when requesting detailed records. A spectrum-only
stream does not build those detail buffers.

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
