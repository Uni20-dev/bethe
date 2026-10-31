# Direct thermodynamic spinon moments

[XXX family](xxx-family.md) · [Selected spectra](xxx-spinon-sectors.md) ·
[C++ two-spinon moments](xxx-two-spinon-moments.md)

`bethe-xxx-moments-thermo` integrates moments directly, without first sampling a
frequency plot. Choose `--momentum Q` for a fixed momentum in radians/site, or
`--global` to integrate over momentum too. The default is the two-spinon sector;
`--sectors 2,4,6` includes the requested contributions and their partial sum.
Every contribution already has its physical normalization.

## Definitions and reference denominators

For the longitudinal channel, the output moments are

```math
m_n^{(r)}(q)=\int_0^\infty\frac{d\omega}{2\pi}\,
\omega^n S_r^{zz}(q,\omega),\qquad
I_n^{(r)}=\int_0^{2\pi}\frac{dq}{2\pi}\,m_n^{(r)}(q).
```

Here $`r`$ is the spinon sector and $`n=0,1`$. The full-response reference
values, summed over all sectors, are

```math
I_0=\frac14,\qquad
I_1=\frac{2J}{3}\left(\ln2-\frac14\right),\qquad
m_1(q)=I_1(1-\cos q).
```

The `fraction` column divides the computed contribution by the corresponding
full-response denominator. No multiplier is fitted to a sum rule. A fixed
momentum has no supplied exact zeroth-moment denominator: `sum_rule` and
`fraction` are absent for that row. In particular, $`m_0(q)`$ is not generally
$`1/4`$.

`--channel raising` multiplies both moments and their reference denominators
by two; covariance scales by four and fractions stay unchanged. Changing
`--exchange J` leaves the zeroth moment unchanged and scales the first moment
by $`J`$.

## Start with a finite-domain calculation

```sh
bethe-xxx-moments-thermo --global --diagnostics
```

The two-spinon calculation uses relative rapidity $`g=\beta_2-\beta_1\gt0`$,
with a default cutoff of 48 and order 16 in panels of width at most two. The
common rapidity shift is integrated over its full range. Compare higher orders
at the same cutoff, then extend the cutoff; the transition-kernel diagnostic
alone does not estimate either outer-quadrature error or the omitted tail.

A small two-plus-four execution check is:

```sh
bethe-xxx-moments-thermo --momentum 1.5 --sectors 2,4 \
  --precision long-double --two-order 4 --two-gap-cutoff 4 \
  --four-order 2 --four-cutoff 2 \
  --density-relative-tolerance 1e-4 --density-absolute-tolerance 1e-24 \
  --threads 2 --diagnostics
```

This deliberately short calculation checks the interface, not convergence.
The nonzero absolute density tolerance permits a meaningful check near a
vanishing density, where a relative-only test can fail. It is an explicit
accuracy control, recorded in metadata; it does not skip failed nodes.

The cutoff definitions differ:

| Sector and calculation | Integration domain |
| --- | --- |
| Two, fixed or global | Relative gap below `--two-gap-cutoff`; common shift unrestricted |
| Four, fixed momentum | Soft log-radius controlled by `--four-cutoff` |
| Four, global | Three positive ordered gaps bounded by `--four-gap-cutoffs A,B,C`; common shift unrestricted |
| Six, fixed or global | Every absolute rapidity lies in `(-R,R)`, with `--six-cutoff R` |

Equal numbers in those options do not describe equal domains. The reported
sum combines the requested finite-domain estimates; it contains no automatic
rapidity-tail correction.

The same `--four-cutoff` also need not select identical domains in a spectrum
and a fixed-momentum moment: the spectral radius bound depends on frequency.
See the [domain comparison](xxx-four-spinon-thermodynamics.md#keep-numerical-checks-separate) before
using frequency-integrated spectral data to check a direct moment.

## Refine quadrature and sampling separately

For four-spinons, `--four-refine-order` requests four controls: coarse and fine
orders on the original and extended domains. Pair it with
`--four-extended-cutoff` at fixed momentum, or with
`--four-extended-gap-cutoffs A,B,C` for global moments. The result uses the
fine order on the extended domain. `quadrature_difference` is the fine/coarse
change there; `cutoff_change` is the signed change between domains at fine
order. They are separate comparisons, not a certified total error bar.

Six-spinons use the same explicit mixture and recovery controls as the
[selected-spectrum frontend](xxx-spinon-sectors.md). One configuration supplies
both moment observations, and their covariance is retained. For example:

```sh
bethe-xxx-moments-thermo \
  --momentum 1.57079632679489661923132169163975 --sectors 6 \
  --precision bits-384 --six-cutoff 8 --replicates 2 --samples 4 \
  --seed 42 --threads 2 --diagnostics --replicate-output --covariance
```

These eight draws are only an execution check. Increase independent work,
inspect influential observations and change proposal weights at fixed domain.
Separately extend the rapidity cutoff: a small replicate error does not bound
the missing tail. Changes in arithmetic precision and contour budgets provide
further checks, distinct from both sampling and domain truncation.

With MPFR, the shared `--density-guard-bits` option can reduce W4/W6 density
cost while preserving a higher geometry/output precision. Its bounds, lower
comparison precision and validation procedure are described in the
[spectrum precision controls](xxx-spinon-sectors.md#precision-and-density-recovery).
It applies to both fixed-momentum and global moments; the recovery checks and
finite-domain interpretation are unchanged.

## Tables and failures

The frontend declares four exportable tables:

| Table | Contents | Availability / screen option |
| --- | --- | --- |
| `moments` | Sector, moment, value, denominator, fraction and status | Always; primary output |
| `diagnostics` | Kernel, contour, arithmetic, quadrature, cutoff and sampling diagnostics | Always / `--diagnostics` |
| `replicates` | Six-spinon replicate means, sample counts, maxima and failures | Sector 6 / `--replicate-output` |
| `covariance` | The 2-by-2 covariance of six-spinon moment means | Sector 6 / `--covariance` |

Use, for example, `--csv moments.csv --csv-table diagnostics=checks.csv`.
Named table exports do not enable their screen display. Unknown or unavailable
tables are errors, detected before opening export files. JSON/CSV/TSV exports
record the numerical controls and run metadata.

A failed requested sector has no value, and its partial sum is marked
`incomplete`; successful sectors remain visible. No failed observation or
replicate is removed to manufacture a survivor-only estimate. Numerical
incompleteness exits with status 2; invalid arguments exit with status 1.

At zero momentum the response vanishes and the first-moment fraction is
undefined ($`0/0`$). At $`q=\pi`$ the unrestricted static response diverges;
a finite-cutoff zeroth moment must not be interpreted as its converged value.
The result carries `unrestricted_static_divergence` when the momentum equals
$`\pi`$ at working precision.

For library composition, use the owned moment queries' `evaluate_async`
methods with a shared execution policy and wait for the final returned result.
Independent sectors share the scheduler without draining unrelated work.
See [Async execution](async-execution.md),
[four-spinon integrals](xxx-four-spinon-thermodynamics.md) and
[six-spinon moments](xxx-six-spinon-moments.md).
