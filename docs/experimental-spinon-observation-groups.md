# Shared samples and annular spinon observations

`bethe::experimental::xxx::SixSpinonObservationGroup<Geometry>` prepares a
bounded group of six-spinon spectral requests from one exact uniform identity.
It uses the public importance proposal, energy-shell solver and density
recovery. A request names momentum, frequency, momentum section and a rapidity
annulus. The preparation retains both labelled roots and distinguishes a
resolved support zero, a masked zero and unresolved geometry.

This checkpoint supplies numerical preparation and observation evaluation.
Persistent group records and their integration with the
[density-control store](experimental-spinon-recovery-store.md) are still being
implemented. This interface alone does not make an integration restartable.

## Measure and masks

The proposal generates four free rapidities in its full outer box. For a
complete six-rapidity configuration, the mask retains

```math
R_{\mathrm{inner}}\le\max_j|\beta_j|\lt R_{\mathrm{outer}}.
```

The proposal is unchanged when a mask is applied. In particular, a compact
observation and a remainder observation use the same free coordinates, full
mixture density and original sample denominator. Their sum reconstructs the
unsplit observation, and their covariance must be retained when combining
uncertainties. Drawing independently from the two regions or normalizing only
by retained samples would be a different estimator.

The complete fifteen-chart balance factor and exchange-energy Jacobian come
from the public proposal and surface implementation. Physical density recovery
already includes its single labelled-momentum factorial. No additional factorial
or Fourier factor is introduced by this adapter. One requested section is only
that section's contribution; sum sections 0, 1 and 2 for the full six-spinon
point at the chosen cutoff.

## Prepare once, evaluate checked roots

Construct a `SpinonJointProposal<G>` and a vector of
`SpinonSectionRequest<G>{id, momentum, omega, inner_cutoff, outer_cutoff, section}`.
Then use:

```cpp
namespace x = bethe::experimental::xxx;
x::SixSpinonObservationGroup<G> group(proposal, exchange, requests);
bethe::MonteCarloStream<5> stream(seed, replicate, uniform_bits);
auto prepared = group.prepare("sample-0", stream.next());
auto observation = group.observation(prepared, density_recovery);
```

Here `density_recovery` is an existing `XXXSpinonDensityRecovery<Low, High>`.
Every geometry scalar has one common precision, large enough to retain the
largest density guard precision. Runtime precision is explicit throughout;
the ambient MPFR default is neither read nor changed.

`prepared.identity()` retains the exact integer words and declared uniform grid.
`prepared.draw()` exposes the shared free coordinates and proposal component.
`prepared.sections()` retains each unique solved surface, both labelled roots,
geometry working precision and refinement count, balanced log weights or a
failure diagnostic. `prepared.rows()` maps every request, in its original order,
to a surface and mask decision. A preparation belongs to its original group or
a copy of that group; combining it with an independently constructed group is
an error. Copies share immutable owned data, and moved-from preparations cannot
be used.

Identical momentum/frequency/section requests share geometry even when their
masks differ. Each selected geometric root is evaluated once and reused across
the masks which contain it. Work metadata counts unique evaluated roots and
density controls, not a duplicated count for every output column.

For K requests, a successful observation has 3K columns: K physical values,
K contour diagnostics, then K arithmetic diagnostics. These are compatible
with the public vector-valued Monte Carlo and replicate-covariance collectors.
The diagnostics are not rigorous error bars or estimates of omitted support.
An unsuccessful observation has no value vector; do not insert zeros or compute
a mean over surviving observations.

## Numerical and API limits

- A group contains 1–32 uniquely named requests. Momentum is strictly between
  zero and 2 pi; the zero-momentum inelastic response is outside this adapter.
  Frequency and exchange are positive and finite. Sections are 0, 1 or 2.
- Masks satisfy `0 <= inner < outer <= proposal.cutoff()`. Both solved roots
  have the same maximum rapidity magnitude, so a mask retains both or neither.
- Exact lower thresholds require a limiting calculation. Frequencies below
  the lower boundary, above the sector upper boundary, or at the vanishing
  outer upper edge produce resolved zeros. Ambiguous internal surface edges
  remain geometry failures.
- A mask boundary within the geometry's `sqrt(epsilon)` screening scale is
  unresolved, rather than silently assigned to one region. This conservative
  check is not an interval certificate.
- The root callback seam `observe(prepared, callback, stop)` exists for bound
  persistence. It receives surface/root indices, rapidities, the balanced log
  weight and cancellation token, and must return a checked weighted observation
  with three finite nonnegative fields. Exceptions propagate. Cancellation or a
  failed root prevents a complete group value.

The adapter owns no scheduler, filesystem or stream controller. Statistical
independence between sample identities remains the caller's sampling contract;
correlations between columns from the same identity must be preserved.

## Current checks

Fresh deterministic samples exercise the original momentum and energy equations,
both root permutations, retained support zeros, mask partitioning and covariance
reconstruction with all 24 draws in the denominator. Native and MPFR observations
are compared with the public full-point evaluator after summing sections.
Those comparisons test adapter equivalence, not a new independent vertex-operator
derivation. Separate tests retain failed densities, cancellation and ambiguous
masks, and reject mixed precision and preparation ownership.
