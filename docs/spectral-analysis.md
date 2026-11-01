# Reading spectra and calculating spatial correlations

[Spectral quadrature and time transforms](spectral-quadrature.md) ·
[Command-line exports](command-line.md)

The standalone Python modules in `scripts/` provide precision-preserving table
readers and linear spatial/cutoff operators. The reader uses only Python's
standard library. Numerical operators use NumPy binary64, with the same
[dependency environment](spectral-quadrature.md)
as the time-transform tools. They do not require Python bindings or a private
calculation archive.

From a source checkout, run the analytic example:

```sh
python -m pip install -r scripts/requirements-analysis.txt
python scripts/example_spectral_analysis.py
```

It calculates on-site, nearest- and next-nearest-neighbour correlations from
a positive test structure factor, then compares nested cutoffs using correlated
replicates. To inspect a fresh application export as well:

```sh
bethe-xxx-structure-factor 8 --format json --roots > finite.json
python scripts/example_spectral_analysis.py --tables finite.json
```

## Read the values that were exported

Add the checkout's `scripts` directory to your Python import path, or run your
analysis script from that directory:

```python
from spectral_tables import read_tables, group_rows, index_rows

export = read_tables("finite.json")
spectrum = export.tables["spectrum"]
print(export.sha256, dict(spectrum.metadata), spectrum.summary)
```

`read_tables` accepts combined `--format json` output. For a single
`--json-table TABLE=FILE` export, pass `table_name="TABLE"`. An optional
`sha256=...` checks the exact bytes read. `expected_metadata` maps table names
to required metadata fields, for example
`{"spectrum": {"Channel": "zz"}}`; a missing table or different convention
is an error. Choose fields from that application's exported schema rather
than assuming every model has the same columns or normalization.

The reader validates column types, integer ranges, nullability, row widths and
precision encodings. Duplicate JSON keys and duplicate column identifiers are
errors. It retains every row, including failed states, and keeps `first_row`
for a sink attached partway through a calculation. A missing final `summary`
remains absent. Envelope status `complete` means the JSON output finished;
it does **not** mean all states converged or all physical families were included.

Fixed and runtime-precision real cells become `RealValue` objects:

| Access | Meaning |
| --- | --- |
| `value.text` | Original numerical spelling, including signed zero |
| `value.precision` | Working binary bits, or the `"exact"` tag |
| `value.decimal()` | Decimal value without passing through binary64 |
| `value.fraction()` | Exact rational represented by the exported text |
| `value.finite_decimal()` | Decimal value, rejecting nonfinite diagnostics |
| `value.float64()` | Explicit binary64 conversion; rejects overflow, nonfinite values and a nonzero value rounding to zero |

Exact MPFR rationals such as `1/3` are retained as written. Use `fraction()`;
`decimal()` rejects a nonterminating rational representation rather than
silently selecting a decimal rounding precision. Optional missing values stay
`None`; `nan`/`inf` diagnostics are retained but rejected by finite conversion.
Wide integers and native half-integers are read exactly. None of these checks
establishes that a value is accurate to its advertised working precision.
Decimal **arithmetic**, as distinct from parsing, still needs an explicit
`decimal.localcontext()` precision appropriate to your calculation.

Tables, rows and metadata are immutable. For finite-state exports with the
corresponding columns, `index_rows(table, ["candidate_id"])` requires unique
state identities. For roots use `index_rows(roots, ["candidate_id", "root_index"])`;
`group_rows(roots, ["candidate_id"])` retains all roots for each state.
Repeated calculations must keep their run/precision identity and source hashes.
No duplicate is silently replaced by a preferred recovery result, and no
failure is converted to a zero spectral line. Family labels, root-coordinate
conventions and separate tiny-deviation columns remain the producer's contract.

## Fourier transform the zeroth frequency moment

For a full-zone convention, the equal-time correlator is

```math
C(r)=\int_{-\pi}^{\pi}\frac{dq}{2\pi}\,e^{iqr}m_0(q),\qquad
m_0(q)=\int_0^\infty\frac{d\omega}{2\pi}\,S(q,\omega).
```

`spatial_weights(q, separations, measure)` returns the matrix multiplying
the supplied `m0` values. Momenta are in radians and separations are lattice
integers. Supply the momentum quadrature measure explicitly, including its
normalization. For a periodic N-site chain this is `1/N` at each of all N
allowed momenta. There is no additional frequency factor for Lehmann weights
already integrated over `dω/(2π)`.

For a reflection-even spectrum, `even=True` uses `cos(q*r)` and half-zone
weights must include `dq/pi`. This is a physical symmetry assumption made by
the caller. No endpoint values are inserted, and an incomplete momentum range
returns only its measured contribution. In particular, on-site intensity,
neighbour correlations and frequency first moments are different observables.

```python
from spectral_spatial import spatial_weights
from spectral_time_transform import statistics

W = spatial_weights(q, [0, 1, 2], momentum_measure, even=True)
result = statistics(W, m0_mean, m0_covariance_of_mean)
```

`statistics` retains the full input covariance, including shared momentum or
sector samples. To compare two observables, pass the difference of their
operator rows before propagating covariance. To transform a time signal at
each momentum, first obtain it with the frequency quadrature operators and
compose the two linear operations; their measures remain separate.

## Cumulative domains without changing the estimator

`cumulative_weights(coordinates, cutoffs, measure)` implements the inclusive
mask `coordinate <= cutoff`. A rapidity span `max(beta)-min(beta)` is different
from a depth `max(abs(beta))` or an energy cutoff. Compute the intended
coordinate explicitly and retain the original domain in the result's provenance.
The operator compares binary64 coordinates; resolve delicate boundary decisions
at the source precision before reducing to this representation.

For raw Monte Carlo observations, include the **original** sample denominator
in `measure`, usually `1/N`. Physical-zero observations still count toward N.
Do not divide by the number of observations left after applying the mask.
For already integrated bin means, supply unit measure. An unrecorded overflow
bin or missing observation cannot be filled with zero.

All cutoffs derived from the same observations are correlated. Apply the
operator to complete replicate vectors or their full covariance of the mean;
subtract operator rows to calculate a shell. Sum means and covariance matrices
only for explicitly independent domains. Correlated domains instead need their
joint observations/covariance. Failed or incomplete replicates require recovery
or an explicit incomplete result, not an estimate based only on successes.

These are sampling errors, separate from frequency/momentum discretization,
finite rapidity domains and numerical evaluation errors. The example and tests
check analytical Fourier coefficients, the Fourier sign, sample-level covariance,
nested cutoffs, unchanged denominators, failures, decimal values and source
identity. Frozen catalogue overlays, machine paths, shell manifests, empirical
finite-size fits, archive collectors and publication plotting stay in research
adapters; they are not part of these general tools.
