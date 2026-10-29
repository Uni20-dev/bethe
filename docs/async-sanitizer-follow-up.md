# Async dependency sanitizer findings

**Historical investigation, 26–27 October 2026; resolved upstream.**
For the current sanitizer configuration, symbol preparation and suppression
policy, use [testing and CI](testing.md#continuous-integration). The counts and
revisions below describe the original checkpoints.

The integration checkpoint used Uni20 `04163b9b84bb71830b2305a2a346d526ecb946ab`, including the
sanitizer boundaries in [Uni20 PR #97](https://github.com/Uni20-dev/uni20/pull/97)
and [PR #99](https://github.com/Uni20-dev/uni20/pull/99),
and notification-based waits in [Uni20 PR #98](https://github.com/Uni20-dev/uni20/pull/98),
with fetched oneTBB 2022.3.0. The original 26 October 2026 integration at Uni20 `1e7c4a96` exposed
two dependency findings under Clang 19, Debug (`-O1 -g1`),
address/undefined-behaviour sanitizers and leak detection. They reproduce without Bethe headers,
numerical kernels or callbacks. The maintainer accepts Uni20's process-lifetime
diagnostic metadata, the existing oneTBB instrumentation exclusion, and the
narrow MPFR cache suppression policy described below. The sanitizer preset
selects these boundaries explicitly; they require no Bethe scheduler workaround.

With the merged PR #99 policy wired through CTest, all **1,316 routine sanitizer
checks pass**, with zero skips, in the 27 October local run (180.69 seconds).
This check reuses existing numerical binaries: the dependency update changes
CMake, tests, and documentation, with no numerical or scheduler-header changes.
Fresh all-target builds remain part of GitHub CI. The merged upstream cache and
sanitizer-option regressions also pass.

## Uni20 debug DAG allocations

Constructing and destroying `uni20::async::Async<int> ready(7)` reports one
56-byte `NodeInfo` allocation on this x86-64 build. With `UNI20_DEBUG_DAG`
enabled, `EpochQueue::initialize_node` calls `NodeInfo::create`; the class
explicitly documents process-lifetime allocations that are never deallocated.
LeakSanitizer therefore reports these even after the Async value, buffers and
tasks are gone. This is distinct from an outstanding Bethe calculation or
retained result payload. The current Debug preset enables this diagnostic path;
normal Release builds disable it.

The lifetime policy is intentional and remains unchanged. Uni20's
`tests/sanitizers/lsan.supp`
matches allocation stacks in the debug-metadata header
`uni20/async/async_node.hpp`. Clang with reduced debug information may omit the
`NodeInfo` qualification on inlined `create<T>` frames, so a qualified function
name alone is insufficient.
Leak detection stays enabled; unrelated leaks must still fail. A complementary
run with `UNI20_DEBUG_DAG=OFF` can check without retained diagnostic pointers.

## oneTBB sanitizer suspend/resume path

A direct `task_arena`/`task_group` call to `task::suspend`, immediately resumed
with `task::resume`, reports:

```text
co_context.h:260:31: runtime error: reference binding to null pointer of type 'thread_data'
```

The sanitizer build uses oneTBB's thread-backed resumable-task path. In
`swap_coroutine`, `governor::set_thread_data(*prev_coroutine.my_thread_data)`
binds through a null pointer. The nearby destruction path wakes a suspended
coroutine without installing thread data. That explains a plausible teardown
route; the direct reproducer establishes the dependency boundary, rather than
proving every possible cause.

Uni20 has had a policy of excluding the compiled oneTBB implementation from
instrumentation since its first TBB scheduler. Its old source-property
implementation did not reach the dependency directory; the correction applies
a private compile option to the actual TBB build target. Bethe now sets
`UNI20_SANITIZE` in addition to its own global flags to activate that policy.
An upstream fix remains useful,
but is not a prerequisite for testing with this accepted boundary.
The provider retains its thread-backed resumable tasks so that disabling its
instrumentation does not switch instrumented application code through ASan's
incompletely supported `makecontext`/`swapcontext` path.

## Reproduce

From a configured `ci-sanitizers` build directory:

```sh
cmake --build . --target bethe_async_dependency_probe -j2
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
DEBUGINFOD_URLS= \
LSAN_OPTIONS=suppressions=/path/to/uni20/tests/sanitizers/lsan.supp:print_suppressions=0 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.Uni20DebugNodeLifetime
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
DEBUGINFOD_URLS= \
LSAN_OPTIONS=suppressions=/path/to/uni20/tests/sanitizers/lsan.supp:print_suppressions=0 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.TbbSuspendResume
```

The probe is opt-in; the regular Async tests already cover these paths and
remain enabled in CI. The pinned Uni20 revision contains the target-scoped
exclusion, so no local dependency override is required. The exclusion does not
constitute sanitizer validation of oneTBB itself; Bethe and Uni20 code, including
provider headers instantiated by them, remain instrumented.

For CTest runs, the sanitizer preset enables
`BETHE_TEST_UNI20_LSAN_SUPPRESSIONS`. CMake resolves the file from the source
directory of the actual `uni20_core` target and sets `LSAN_OPTIONS` for both
discovered GoogleTests and CLI/installation checks. This handles fetched,
locally overridden, and parent-supplied Uni20 sources without assuming a build
directory name. Bethe does not maintain its own suppression file. Set this
option to `OFF` to run an unsuppressed audit.

The preset also clears `DEBUGINFOD_URLS`. An inherited Ubuntu debuginfod URL
caused CLI timeouts while `llvm-symbolizer` waited for remote HTTPS connections;
the numerical processes were waiting for symbolization. Local debug information
remains available, and disabling those downloads does not disable sanitizer
checks or leak detection.
This matches [Uni20's documented test-environment policy](https://github.com/Uni20-dev/uni20/blob/04163b9b84bb71830b2305a2a346d526ecb946ab/docs/development/testing.md):
keep the setting in test commands or presets, without changing library behavior
or the user's shell configuration. System-library stack traces may have less
detail when their debug symbols are not installed locally.

## MPFR thread-cache reports and suppression policy

Before PR #99, with the TBB boundary corrected and only the debug-metadata
suppression selected, the 27 October 2026 validation against Uni20 `30af0808`
passed 1,311 of 1,316 tests, with no skips, in 170.98 seconds.
The five failures reported MPFR/GMP cache allocations after threads exit:

- `RootSolverMpfrShared.AlternatingPrecisionsAndThreadedScans`
- `XxxMpfr.RootOnlyConcurrentPrecisionsAndThreadedExcitations`
- `AsyncMpfr.ConcurrentPeriodicAndOpenScansOwnTheirPrecision`
- `RootSolverPrecision.cli`
- `XXXPbcPrecision.cli`

The numerical and precision assertions in the three unit tests passed before
LeakSanitizer reported allocations at process exit. PR #99's policy suppresses
the observed internal C frames `mpfr_atan`, `mpfr_cos`, and `mpfr_cache`, with
start and end anchors. The anchors matter: a substring rule also matches Uni20's templated wrapper
`mpreal_access::apply<&mpfr_atan>` and can hide leaked owning result storage
allocated before the C function is called. Upstream regression tests verify
that real scalar and ordinary heap leaks still fail. A function-based exclusion
can nevertheless hide other internal allocations within the named MPFR frames;
it is not a blanket claim about MPFR's ownership correctness.

Suppression statistics are disabled so accepted allocations do not disturb CLI
tests that require empty stderr; actual sanitizer diagnostics still fail.

The full runs also exposed concurrent use of the process-wide, single-threaded
`DebugScheduler` in the XXX MPFR root and form-factor tests' serial reference
calculations. Each concurrent reference callback now owns a separate
`DebugScheduler`; the references must not share the default scheduler.

A standalone oneTBB/MPFR control reproduces this with `mpfr_atan` at 333 bits
after `task::suspend`/`resume`, even after clearing every `mpfr_t` value.
Calling `mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE)` on the executing thread removes
the report. That control establishes the cache-lifetime cause; freeing caches
after every operation is not the proposed implementation.

[MPFR's memory-handling contract](https://www.mpfr.org/mpfr-current/mpfr.html#Memory-Handling)
asks users to free local caches before a thread terminates. Clearing only the
application thread's cache after a scheduler wait does not clear caches created
on other threads. Exit-time reports alone do not require a runtime cleanup
mechanism: persistent worker threads reuse their caches. Thread turnover and
the sanitizer provider's auxiliary resumable-task threads need separate
consideration. Suppression changes test reporting, not runtime cache lifetime;
these allocations are distinct from execution failures or numerical errors.

## GitHub runner findings

The sanitizer runner also needs `llvm-19`, which supplies
`llvm-symbolizer-19`. Installing Clang and compiler-rt with
`--no-install-recommends` did not supply a usable external symbolizer: CI warned
about its invalid path and emitted address-only stacks. Source-based LSan
suppressions cannot match those stacks. The workflow now installs the package,
checks its version and sets `ASAN_SYMBOLIZER_PATH` explicitly. This addresses
symbolization, without expanding the suppression.
The subsequent GitHub run at `086bef4` confirms 1,311/1,316 sanitizer tests pass,
with only the same five MPFR cache reports remaining.

The wide Release job exposed a separate scheduler assertion in
`TbbScheduler::ExecutionScope::~ExecutionScope`: `current_` is null instead of
the active scope. Both concurrent `AsyncMpfr` tests aborted; the other 1,823
routine cases passed. Repeating those tests locally with a two-CPU affinity
reproduced the assertion on iteration seven, despite the earlier full local
suite passing. This failure affected the previous pin, `3c0e7138`, independently
of MPFR cache retention.

## Notification-based scheduler correction

Uni20 PR #98 removes execution-scope accounting and the wait watchdog. Pending
waits instead require readiness notifications; both reader and writer waits
supply them. Bethe's callers need no API changes. The separate test correction
giving each concurrent serial reference its own `DebugScheduler` remains
necessary.

Before the merge, targeted GCC 13 Release checks against these Uni20 headers
and existing build configurations and provider libraries passed:

- All three `AsyncMpfr` tests, repeated 50 times with two-CPU affinity: 150 passes.
- The corrected concurrent XXX reference test, repeated 30 times: 30 passes.
- The complete `XxxMpfr` test executable: nine passes.
- Uni20's updated TBB module plus the direct suspend/resume regression:
  50 passes and one NUMA-topology skip. The regression fails with the previous
  scheduler's execution-scope assertion.

The complete local routine suites were then rebuilt and run against the fetched
merge commit, `30af0808`, with no local Uni20 source override:

| Configuration | Passed / total | Elapsed time |
| --- | ---: | ---: |
| GCC 13 Release, native scalars | 1,286 / 1,286 | 70.68 s |
| GCC 13 Release, binary80/binary128/MPFR | 1,825 / 1,825 | 372.33 s |
| Clang 19 ASan/UBSan, native scalars and MPFR | 1,311 / 1,316 | 170.98 s |

All three suites have zero skips and include CLI and installation tests. The
wide build reuses an optimized MPLAPACK 3.0.0 installation; GitHub CI builds the
pinned provider from source. The five scheduled long-running numerical tests
are outside these routine suites. The sanitizer failures are exactly the five
cache reports listed above: the scheduler correction resolves the execution
failure without changing their handling or disabling leak detection.
