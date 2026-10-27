# Async dependency sanitizer findings

Bethe pins Uni20 `30af0808380cdc7d467fa584edd6b3bec4703285`, including the
sanitizer boundaries in [Uni20 PR #97](https://github.com/Uni20-dev/uni20/pull/97)
and notification-based waits in [Uni20 PR #98](https://github.com/Uni20-dev/uni20/pull/98),
with fetched oneTBB 2022.3.0. The original 26 October 2026 integration at Uni20 `1e7c4a96` exposed
two dependency findings under Clang 19, Debug (`-O1 -g1`),
address/undefined-behaviour sanitizers and leak detection. They reproduce without Bethe headers,
numerical kernels or callbacks. The maintainer accepts Uni20's process-lifetime
diagnostic metadata and the existing oneTBB instrumentation exclusion. The
sanitizer preset now selects those boundaries explicitly; neither requires a
Bethe scheduler workaround.

## Uni20 debug DAG allocations

Constructing and destroying `uni20::async::Async<int> ready(7)` reports one
56-byte `NodeInfo` allocation on this x86-64 build. With `UNI20_DEBUG_DAG`
enabled, `EpochQueue::initialize_node` calls `NodeInfo::create`; the class
explicitly documents process-lifetime allocations that are never deallocated.
LeakSanitizer therefore reports these even after the Async value, buffers and
tasks are gone. This is distinct from an outstanding Bethe calculation or
retained result payload. The current Debug preset enables this diagnostic path;
normal Release builds disable it.

The lifetime policy is intentional and remains unchanged. `tests/lsan.supp`
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
LSAN_OPTIONS=suppressions=/path/to/bethe/tests/lsan.supp:print_suppressions=0 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.Uni20DebugNodeLifetime
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
DEBUGINFOD_URLS= \
LSAN_OPTIONS=suppressions=/path/to/bethe/tests/lsan.supp:print_suppressions=0 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.TbbSuspendResume
```

The probe is opt-in; the regular Async tests already cover these paths and
remain enabled in CI. The pinned Uni20 revision contains the target-scoped
exclusion, so no local dependency override is required. The exclusion does not
constitute sanitizer validation of oneTBB itself; Bethe and Uni20 code, including
provider headers instantiated by them, remain instrumented.

The preset also clears `DEBUGINFOD_URLS`. An inherited Ubuntu debuginfod URL
caused CLI timeouts while `llvm-symbolizer` waited for remote HTTPS connections;
the numerical processes were waiting for symbolization. Local debug information
remains available, and disabling those downloads does not disable sanitizer
checks or leak detection.
This matches [Uni20's documented test-environment policy](https://github.com/Uni20-dev/uni20/blob/bc14c127d28cfd10da1621f1805910f044c81683/docs/development/testing.md):
keep the setting in test commands or presets, without changing library behavior
or the user's shell configuration. System-library stack traces may have less
detail when their debug symbols are not installed locally.

## Remaining MPFR thread-cache reports

With the TBB boundary corrected and the debug-metadata suppression selected,
the 27 October 2026 validation against the fetched pin (no local Uni20 source
override) builds all frontends and tests.
The full routine sanitizer suite passes 1,311 of 1,316 tests, with no skips,
in 170.98 seconds. The five failures report MPFR/GMP cache allocations after threads exit:

- `RootSolverMpfrShared.AlternatingPrecisionsAndThreadedScans`
- `XxxMpfr.RootOnlyConcurrentPrecisionsAndThreadedExcitations`
- `AsyncMpfr.ConcurrentPeriodicAndOpenScansOwnTheirPrecision`
- `RootSolverPrecision.cli`
- `XXXPbcPrecision.cli`

These allocations are not suppressed. The numerical and precision assertions
in the three unit tests pass before LeakSanitizer reports the allocations at
process exit. Suppression statistics are disabled in the preset so accepted
debug metadata does not disturb CLI tests that require empty stderr; actual
sanitizer diagnostics still fail those tests.

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
consideration. These known reports remain unsuppressed while the cache policy
is reviewed; they are distinct from execution failures or numerical errors.

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
