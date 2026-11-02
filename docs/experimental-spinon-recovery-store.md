# Restartable spinon density chunks

`bethe::experimental::xxx::SpinonRecoveryStore<Low, High>` saves completed
two-, four- and six-spinon density controls in a bounded input chunk. Reopening
the same chunk resumes an unfinished point at its next missing control. A
completed result is revalidated and returned without evaluating more density
kernels. It uses the supported [recovery driver](spinon-recovery-progress.md),
including its precision, contour and cross-control agreement checks.

This is an opt-in POSIX component. It does not submit jobs, choose a numerical
policy or estimate a spectrum. The caller supplies a directory and retains
responsibility for its storage location and lifetime.

## A small restart example

From a build directory, configure the source and run the supplied example:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_EXAMPLES=ON
cmake --build . --target bethe-experimental-recovery-store-example
./experimental/bethe-experimental-recovery-store-example ./recovery-example --stop-after-control
./experimental/bethe-experimental-recovery-store-example ./recovery-example
./experimental/bethe-experimental-recovery-store-example ./recovery-example
```

The first command exits with status 2 after saving one control. The second
finishes that two-spinon point and calculates fresh four- and six-spinon points.
The third reports zero new controls for every point. Numerical rejection also
returns status 2; the diagnostic distinguishes it from deliberate cancellation.
The example uses binary64 with a wider `long double` guard, and exits with
status 77 on platforms where `long double` is not wider. These small interior
points demonstrate restart behavior, not sufficient precision for arbitrary
soft configurations or a production integration.

The source is [experimental/examples/spinon_recovery_store.cpp](../experimental/examples/spinon_recovery_store.cpp).
Link `bethe::experimental` to use the header in your own build-tree consumer.
The example does not require GoogleTest. Experimental headers and executables
are not added to the supported installation.

## Bind inputs and numerical policy

The constructor takes a `RecoveryChunkBinding` containing:

- `producer`: the executable/kernel revision or another identifier that changes
  when the numerical implementation changes;
- `provenance`: the caller's request provenance;
- `policy`: the complete `SpinonRecoveryPolicy`, including its original time
  budget, tolerances, precisions and both contour controls at every stage;
- `points`: an ordered list of unique, nonempty IDs and original decimal
  rapidity strings.

Both provenance fields must be nonempty. The library compares their exact
contents but cannot infer whether they identify the correct executable.
Input order and decimal spelling are part of the binding, even for equivalent
spellings such as `1` and `1.0`. Native precision and exponent ranges are also
bound, so a change in the platform's `long double` cannot silently reuse a file.
MPFR stages retain explicit precision and never change its ambient default.

A chunk contains 1–64 points and at most 16 precision stages. Serialized data
is limited to 16 MiB, each string to 64 KiB and precision to 65536 bits. Invalid
or coincident rapidities are rejected before computation. The format is
versioned experimental storage, with no promise of compatibility with older
research-worker files.

Each stored scalar carries its own precision, binary exponent and round-trip
decimal fraction. This preserves native subnormals as well as mixed-precision
MPFR controls. Restore checks framing, checksum and input binding, then asks
the numerical driver to validate and reassess the complete ordered history.
An accepted flag or cached estimate is never sufficient evidence. Earlier
contradictions and failed controls survive restart; a failed control is not a
zero density.

## Publication and process failure

The directory contains `owner.lock`, `snapshot.bin`, and at most one temporary
`snapshot.pending`. One snapshot contains the whole bounded chunk. A nonblocking
exclusive `flock` rejects a second owner; the lock inode remains in place after
close or process death.

Publication writes and closes the temporary file, then renames it over the
snapshot in the same directory. After process death, reopening sees the old or
new complete snapshot. A stale temporary file is ignored. A control interrupted
before publication may execute again: recovery does not promise exactly-once
kernel execution. The checksum detects accidental corruption; it is neither
authentication nor protection against deliberate rollback.

There is no `fsync` protocol and no guarantee of durability across power loss,
kernel failure or loss of the storage server. Tests use a local POSIX filesystem.
Before deploying on a network or parallel filesystem, establish its locking
and rename behavior; those filesystems are not certified by these local tests.
Do not replace or unlink a chunk's lock file while it is in use.

An I/O or publication exception poisons that store instance. Close it after
its callers finish, repair the storage problem and reopen it; do not continue
with its unpublished in-memory history. Previously published controls remain
available. An observer exception likewise propagates, with its preceding
publication available for restart.

## Time, cancellation and concurrent calls

The original point deadline is published before admitting its first kernel.
For an unfinished point, elapsed time includes publication, downtime and time
waiting between calls. Restart never resets the budget. Backward wall-clock
movement is rejected rather than interpreted as extra available time. This
uses system wall time across processes, not CPU time.

An already finalized accepted or rejected result can be replayed later without
acquiring a new deadline. Consumed time is frozen at final assessment: publishing
that finalized snapshot and replaying it later cannot expire an accepted result.
Waiting for another point to release the publication mutex after final assessment
also leaves that result unchanged. Publication and observer time for the last control still count, because they
precede final assessment. Cancellation leaves a resumable prefix when the
numerical history is not already terminal. Cancellation and budget checks
occur between kernels; they cannot interrupt one kernel or an in-progress
filesystem operation.

`run(index, observer, stop_token)` publishes each completed control before
calling the observer outside the store mutex. The observer receives a borrowed
const view; copy anything retained beyond the callback. `load(index)` returns
an optional reassessed history without admitting kernels. Distinct points may
run concurrently, but duplicate active calls for one index are rejected and
snapshot publication is serialized. Changing concurrency on restart does not
change point identities.

The store owns no scheduler and must outlive every call. An Async adapter must
drain all admitted sibling work before propagating failure or destroying the
store. Full weighted-observation records, common-sample covariance and annular
selection are separate integration responsibilities.

## Validation

Tests cover partial and completed replay, independent point concurrency,
mixed-precision histories, preserved failures and contradictions, input/policy
mismatches, corruption, duplicate ownership, publication failure and forced
process exits immediately before and after atomic replacement. Deadline tests
include downtime and slow initial publication. These establish the restart
contract; they do not independently validate the underlying physical formulas.

A local one-core cost check of the example's three points took median 6.07 ms
through the recovery driver and 8.16 ms with publication; completed replay took
0.456 ms and executed zero new controls. These include context construction and
use binary64 with a binary80 guard. Eight alternating repetitions follow one
warmup, with identical accepted values. This measures about 2.09 ms of storage
overhead for very small kernels on one local filesystem, not network-filesystem
or production throughput. The [validation record](../experimental/provenance/spinon-recovery-store.json)
records the configurations, frozen source provenance and remaining limits.
