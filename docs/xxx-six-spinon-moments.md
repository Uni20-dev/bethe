# Direct six-spinon moments

`bethe/xxx_six_spinon_moments.hpp` calculates the longitudinal six-spinon
intensity and first moment directly, without integrating a displayed frequency
grid. `ThermodynamicSixSpinonMoments<Geometry,Low,High>` reuses the same physical
W6, bounded recovery and soft-pair proposal as the
[spectral-point calculation](xxx-six-spinon-structure-factor.md). It uses the
shared [Monte Carlo engine](monte-carlo.md) to accumulate both moments and their
covariance from each configuration.

Both fixed-momentum and global moments are available. The installed moments
front end remains integration work. The calculation uses a finite
absolute-rapidity cutoff, not a claim of unrestricted convergence.

## Measure and importance proposal

Let q4 be the complete four-variable joint proposal and g its native marginal.
A draw chooses one of the fifteen unordered pairs uniformly, samples that pair
from g, and samples the other four rapidities from q4. Its actual density is

```math
Q_6(\boldsymbol\beta)=\frac1{15}\sum_{i\lt j}
g(\beta_i)g(\beta_j)q_4(\boldsymbol\beta_{\ne i,j}).
```

Every term includes all overlapping channels of q4. The mixture contains the
background and the requested soft-pair/cluster components. Seven open uniform
coordinates select the native pair, joint component and continuous coordinates.
The complete density determines the weight; the selected component alone does
not. This importance strategy remains an explicit part of the physical adapter.

For labelled rapidities restricted by max(abs(beta)) less than R, a contribution
to the order-r global moment is

```math
Y_r=\frac{W_6(\boldsymbol\beta)}{(2\pi)^2Q_6(\boldsymbol\beta)}
\prod_j\mathrm{sech}\,\beta_j\,
\left[\frac{\pi J}{2}\sum_j\mathrm{sech}\,\beta_j\right]^r,
\qquad r=0,1.
```

The Jacobian converts rapidities to momenta. W6 already contains its single
state-counting factorial; there is no additional permutation correction.
Both moment measures are included. The unrestricted first-moment momentum law
is not imposed on this global integral.

The same cutoff box is used for spectra and these global moments, allowing
matched-domain comparisons. It excludes soft configurations throughout the
spectrum. A gap cutoff with an analytically integrated common shift is a
different domain and needs its own comparison. The common-shift reduction is a
useful alternative; its five-gap importance proposal requires separate validation.

### Fixed-momentum geometry

`detail/xxx_spinon_moment_proposal.hpp` supplies the direct fixed-momentum
geometry used by `fixed_moments(k)`. It chooses uniformly among thirty labelled
charts: eliminate momentum i, draw a distinct rapidity j
from g, and draw the other four from q4. It solves the momentum constraint in
every section and retains the solution only inside the same absolute cutoff
box. A draw with resolved empty support contributes zero and is not redrawn.
Numerically ambiguous support fails explicitly.

The chart intensity against the constrained labelled momentum measure is

```math
h_{i,j}(\boldsymbol\beta)=
q_4(\boldsymbol\beta_{\ne i,j})g(\beta_j)
\prod_{a\ne i}\cosh\beta_a,\qquad i\ne j.
```

Its reciprocal mixture weight is 30 divided by the sum of all thirty
intensities. Eliminating a momentum has unit Jacobian; the cosh factors convert
the five free proposal densities to momentum densities. The moment
observation is this weight times W6 and the sampled energy to power r, divided
by 2*pi. No further factorial or normalization adjustment is needed. Chart
selection is actually sampled, so the measure also applies to nonsymmetric
labelled observables.

Geometry checks cover the full mixture and label orbit, all chart choices,
momentum conservation, analytic periodic box-spline volume and labelled-momentum
integrals. R=1040 conservation is checked at 1800 and 2048 bits; insufficient
384-bit geometry is rejected. These establish the geometric measure, not the
efficiency of this proposal for physical W6 or a converged moment estimate.

## Interface and result

Construct the model with a `SpinonJointProposal<Geometry>`, a
`SpinonRecoveryPolicy` and positive exchange J, as in the spectral example.
`global_async(options, execution)` returns an owned Async result;
`global(options, execution)` waits using the same scheduler. A retained
`UniformSample<7>` can also be passed to `global_observation` for exact replay.
The scheduler remains borrowed until completion.

For fixed momentum, use `model.fixed_moments(k).evaluate_async(options, execution)`
or the synchronous `evaluate` wrapper. The query owns the model state through
completion. `observation(UniformSample<6>, stop)` replays one exact sample.
The result has the same six integral columns and covariance, with m0(k) and
m1(k) replacing I0 and I1. It records the requested momentum and an optional
`first_moment_denominator` and `first_moment_fraction`. The denominator is the
full-response first moment, `(2J/3)(log(2)-1/4)(1-cos(k))`; the numerator is
integrated from the sampled energy and is not forced to follow that law.

`SixSpinonMomentStatus` distinguishes an estimate, an exact physical zero,
no sampled support, numerical failure and cancellation. At k=0 or 2*pi both
moments are exactly zero without sampling; their first-moment fraction is
undefined. At k=pi, `unrestricted_static_divergence` warns that the finite-R
zeroth moment is not a converged static structure factor. No exact zeroth-moment
denominator is invented at other momenta. A complete sample stream that finds
no support has no aggregate estimate; its zero observations remain in the
replicate records for diagnosis.

The result records J, R, proposal parameters and recovery settings. Its
`integral` contains the complete numerical-engine result, with columns:

| Column | Quantity |
|---:|---|
| 0 | Global intensity I0 |
| 1 | Global first moment I1 |
| 2, 3 | Sampled contour diagnostics for I0 and I1 |
| 4, 5 | Sampled arithmetic diagnostics for I0 and I1 |

The 6 by 6 `integral.covariance` is the covariance of the mean, in row-major
order. Entries 0 and 7 give the two moment variances; entry 1 gives their
covariance. Retained replicate means and maxima expose sampling scatter and
unusually influential observations. All numerical quantities remain at the
explicit geometry/accumulation precision.

`fraction`, when present, divides the absolute moments by the full-response
sum rules I0=1/4 and I1=(2J/3)(log(2)-1/4). These are reference denominators,
not fitted multipliers. Fractions share the same cutoff and sampling limitations
as the absolute moments. No fraction or mean is returned after an unrecovered
sample failure or cancellation. Increasing J leaves I0 unchanged and scales
I1 by J; both fractions are invariant under that exchange change.

`xxx_spinon_observation.hpp` shares the checked density weighting between
spectra and moments. Observable adapters supply logarithms of all physical
weights; the helper imports the same rapidities into the declared recovery
contexts and returns values, contour diagnostics and arithmetic diagnostics
separately. It supplies no missing measure or normalization factor.

## Validation and practical limits

Independent controls check the complete fifteen-component PDF, all native-pair
choices, permutation invariance, the labelled momentum volume, total energy and
a nonsymmetric labelled observable. Deep-soft log weights are tested at R=1040
with 160/256-bit MPFR; that geometric test is not a physical R=1040 integration.

An independent 85-decimal-digit Python construction checks a complete physical
I0/I1 observation and a fixed-k m0/m1 observation, including all mixture terms
and moment measures. Physical
integration tests check serial/concurrent agreement, covariance, exchange
scaling, owned lifetimes and explicit failure/cancellation. The optional
`tests/reference_xxx_six_spinon_moments.py` regenerates the physical reference.
Use `--fixed` for the fixed-k reference. The compact
[sampling record](validation/six-spinon-moments-20261031.json) retains the
following pilot settings, replicate summaries, timings and comparisons.

An accurate density does not ensure an accurate integral. Production work must
vary sample budgets, independent streams, proposal parameters, uniform-grid
resolution, geometry precision and rapidity cutoff as appropriate. The numerical
density diagnostics exclude sampling and omitted-tail uncertainty. The R=8
settings used for interface tests and small pilots do not resolve the full
thermodynamic six-spinon weight.

A fresh R=8, J=1 execution pilot used 8 replicates of 16 draws at 384-bit
geometry precision, with 128/192-bit density controls and a declared 224/320-bit
recovery stage. All 128 configurations completed, requiring 256 density-control
attempts. It gave I0=0.00012085 with sampling standard error 0.00007627 and
I1=0.00016657 with standard error 0.00011297. These large errors (about 63% and
68%) show that the pilot budget is insufficient for precise integration. It
used 48.73 seconds wall time and 97.03 seconds CPU time on two workers. The tiny
sampled density diagnostics do not reduce those sampling errors or quantify
the omitted rapidity region. These are execution measurements, not estimates
of the unrestricted sector fractions.

### Compare proposals on new samples

At k=pi/2 and R=8, a 512-draw calibration suggested emphasizing the two-pair
component. A separate comparison froze both proposals, then used 16 new
replicates of 256 draws each (seed 42, replicate IDs 100--115). The same identities
were used between proposals, so their difference has a paired sampling error.
The calibration used IDs 0--7 and was not reused as validation data.

The joint-mixture weights for background, one pair and two pairs were
(0.5, 0.25, 0.25) and (0.1, 0.2, 0.7). Both used native soft fraction 0.5,
one-pair scale/gap 4/0.5 and two-pair scale/gap 8/1, with the precision and
density policy above. These are explicit physical proposal choices, not
adaptation performed by the Monte Carlo engine.

| Proposal | m0(k) | m1(k) | CPU seconds |
|---|---:|---:|---:|
| Original mixture | 0.0001106 +/- 0.0000253 | 0.0002509 +/- 0.0000545 | 1577 |
| More two-pair sampling | 0.0000982 +/- 0.0000128 | 0.0002223 +/- 0.0000245 | 1488 |

Each calculation used two workers on an Intel Core i7-12700KF; the two processes
ran concurrently. Timings describe this pilot, not a portable performance target.

Errors are between-replicate sampling standard errors. The two means agree
within 0.51 and 0.54 paired standard errors. The second proposal's observed
variance times CPU cost is 0.24 and 0.19 of the original for m0 and m1.
Its largest individual draws nevertheless account for 6.8% and 5.3% of the
corresponding sums, compared with 18.8% and 17.7% for the original proposal.
The comparison supports this proposal for further finite-R pilots; it does
not establish a universally efficient mixture or exclude unseen rare regions.
No library default is changed by this measurement.

The earlier 512-draw calculations at pi/4, pi/2 and 3pi/4 also compared the
independently measured first-moment fractions. All pairwise differences are
less than 0.51 paired standard errors, but individual fractional uncertainties
are 33--36%. There is no resolved momentum dependence in this small sample,
and no numerical verification of the unrestricted momentum law. A finite
absolute-rapidity box can itself violate that law. The record also includes
a separate 1024-draw global integral; it is not reconstructed from these three
momenta or from a spectral grid.
