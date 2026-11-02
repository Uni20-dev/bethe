# Shared samples and annular spinon observations

`bethe::experimental::xxx::SixSpinonObservationGroup<Geometry>` prepares a
bounded group of six-spinon spectral requests from one exact uniform identity.
It uses the public importance proposal, energy-shell solver and density
recovery. A request names momentum, frequency, momentum section and a rapidity
annulus. The preparation retains both labelled roots and distinguishes a
resolved support zero, a masked zero and unresolved geometry.

The numerical adapter can be used independently. The POSIX
`SpinonObservationStore<G, Low, High>` adds bounded persistent groups around
the [density-control store](experimental-spinon-recovery-store.md). Both are
experimental build-tree interfaces, enabled by `BETHE_BUILD_EXPERIMENTAL`;
they are not installed with the supported library.

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
  the lower boundary, above the requested momentum section's upper boundary,
  or at its vanishing upper edge produce resolved zeros. Each section's maximum
  occurs at six equal momenta. Other momentum sections can still contribute at
  that frequency, so this check uses the requested section, not the full band's
  maximum. Other ambiguous internal surface edges remain geometry failures.
- A mask boundary within the geometry's `sqrt(epsilon)` screening scale is
  unresolved, rather than silently assigned to one region. This conservative
  check is not an interval certificate.
- The root callback seam `observe(prepared, callback, stop)` exists for bound
  persistence. It receives surface/root indices, rapidities, the balanced log
  weight and cancellation token, and must return a checked weighted observation
  with three finite nonnegative fields. Exceptions propagate. Cancellation or a
  failed root prevents a complete group value.

The numerical adapter owns no scheduler, filesystem or stream controller. Statistical
independence between sample identities remains the caller's sampling contract;
correlations between columns from the same identity must be preserved.

## Persistent chunks and replay

Include `bethe/experimental/spinon_observation_store.hpp`. Supply a directory,
group, ordered vector of `SpinonObservationIdentity{id, uniforms}`, recovery
policy, producer revision and provenance. A chunk holds at most 64 samples and
64 selected labelled roots across those samples. Split larger work into explicit
chunks before opening a store. The prepared requests still have the 32-row
limit. Scalar strings are bounded by 64 KiB and each snapshot by 16 MiB.

The binding retains proposal parameters, exchange and normalization, request
order, uniform bit counts and integer words, prepared geometry, mask decisions,
both roots and balanced weights. Reopening reconstructs preparation and requires
an exact match. Changing precision, proposal, executable binding, request order
or geometry therefore fails visibly. This is a same-configuration restart, not
a migration format across numerical implementations. The producer identifier
must identify the actual executable/kernel revision; the library cannot infer
that provenance for a downstream application.

Finite binary rapidities are written as their **exact terminating decimals**
before entering the density store. A shorter decimal that round-trips at the
geometry precision can fall on the opposite side of a lower-precision rounding
tie. Exact decimals preserve the original binary input for every declared
density control. Extremely large exponents/precision exceeding the bounded
representation are rejected before filesystem creation.
Exactly coincident rapidities require a separate limiting calculation and are
also rejected during density-binding preflight. Distinct stored coordinates
that collide only after conversion to a lower precision retain the density
driver's explicit failure/recovery result.

`run(index)` resumes missing controls and atomically commits a complete weighted
observation or explicit numerical failure. Cancellation and unexpected exceptions
leave the observation uncommitted; completed controls remain recoverable.
`load(index)` returns an optional observation: an absent record is unfinished,
while a present record with absent `values` is a numerical failure. No missing
row is filled with zero. A failed root may leave later roots unattempted; these
remain pending rather than becoming invented failures or zeros.

On replay, the density driver revalidates control histories, and the common
weighting function rederives the observation and compares it with the saved
record. A cached acceptance flag or value is insufficient. Completed density
histories execute no new kernels, although preparation, validation and weighting
are repeated. Metadata counts the retained control history, not newly executed
work. The optional `run(index, progress, stop)` callback observes each newly
completed control after publication; it receives section/root indices and the
recovery result, outside the storage mutex.

The filesystem contains a chunk lock/snapshot and a `controls` subdirectory
with its own lock/snapshot. Temporary replacement files may also exist: at most
six regular files, independent of sample/root count. An all-zero chunk needs no
density store. Locks retain a stable inode; a second owner is rejected. Atomic
rename protects process restarts, not power-loss durability. No fsync or NFS
durability claim is made. A publication failure poisons that store instance;
close it and reopen rather than continuing with uncertain in-memory state.
The underlying recovery store preserves original per-root deadlines, including
downtime for unfinished work. Chunks are not cryptographically authenticated
and do not detect deliberate rollback of the whole dataset.

## Async execution and statistics

`evaluate_spinon_observation_chunk_async(shared_store, collect, execution)`
uses the existing bounded observable grid. It owns the shared store through
all admitted work, collects in input order, and drains siblings before exposing
an exception or cancellation. The borrowed scheduler must outlive completion.
Changing `max_in_flight` on restart changes neither identity nor accumulation
order. Direct concurrent `run` calls on different indices are supported; a
second active call for the same observation is rejected. Callbacks must allow
concurrent execution where the relevant interface documents it.

The returned count includes explicit failed observations. It is not a successful
integration count. Feed all successful columns together into the existing
replicate/covariance tools, preserving the original denominator and independence
assumptions. If any required observation fails, do not report a survivor-only
mean. Stream partitioning, replicate identity and cluster scheduling stay with
the caller; the store does not infer them from a filename or seed.

The small example generates twelve fresh shared samples and nine output columns
(three sections, each with whole/compact/remainder masks). It requires a native
`long double` wider than binary64. From your configured build directory:

```sh
cmake --build . --target bethe-experimental-observation-store-example
./experimental/bethe-experimental-observation-store-example /path/to/chunk --stop-after-control
./experimental/bethe-experimental-observation-store-example /path/to/chunk
./experimental/bethe-experimental-observation-store-example /path/to/chunk
```

The first command that runs the example exits with status 2 after saving one
new control. The second resumes, and the third replays committed observations.
Use a new directory for a new numerical policy or producer revision. This tiny
calculation demonstrates restart and joint statistics; it is not a converged
six-spinon spectrum. No archive or cluster access is required.

## Current checks

Fresh deterministic samples exercise the original momentum and energy equations,
both root permutations, retained support zeros, mask partitioning and covariance
reconstruction with all 24 draws in the denominator. Native and MPFR observations
are compared with the public full-point evaluator after summing sections.
Those comparisons test adapter equivalence, not a new independent vertex-operator
derivation. Separate tests retain failed densities, cancellation and ambiguous
masks, and reject mixed precision and preparation ownership.

Storage tests use native binary64/long-double recovery and MPFR80/128 recovery
with MPFR160 geometry. They include cross-precision rounding ties, process death
after a completed control, replay without kernels, valid-checksum altered results,
publication failure, and held Async work on cancellation/error. The
[validation record](../experimental/provenance/spinon-observation-chunks.json)
separates fresh CMake checks from focused configured-provider builds.

A small serial cost check of the twelve-sample example on local storage gave
median times of 12.0 ms for direct evaluation, 16.2 ms with new persistence and
2.50 ms for reopening/replay. All three produce identical numerical observations;
the first two execute four controls and replay executes none. These are five
interleaved trials of a tiny interior example, including geometry and open-time
validation. They do not measure production throughput or network filesystem cost.
