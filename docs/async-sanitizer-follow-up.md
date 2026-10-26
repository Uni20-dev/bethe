# Async dependency sanitizer findings

Bethe pins Uni20 `3c0e713837eccf0918f0b8f95b5b7109fd5bbc11`, including
[Uni20 PR #97](https://github.com/Uni20-dev/uni20/pull/97), with fetched oneTBB
2022.3.0. The original 26 October 2026 integration at Uni20 `1e7c4a96` exposed
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

## Remaining MPFR cache cleanup

With the TBB boundary corrected and the debug-metadata suppression selected,
the 27 October 2026 validation against the fetched pin (no local Uni20 source
override) builds all frontends and tests and passes both dependency-only probes.
The full routine sanitizer suite passes 1,311 of 1,316 tests, with no skips,
in 191.70 seconds. The five failures report MPFR/GMP cache allocations after threads exit:

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

The full run also exposed concurrent use of the process-wide, single-threaded
`DebugScheduler` in the XXX MPFR test's serial reference calculations. Each
concurrent callback now owns a separate `DebugScheduler`; its numerical
assertions pass on rechecking, leaving the cache report above. The concurrent
references must not share the default scheduler.

A standalone oneTBB/MPFR control reproduces this with `mpfr_atan` at 333 bits
after `task::suspend`/`resume`, even after clearing every `mpfr_t` value.
Calling `mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE)` on the executing thread removes
the report. That control establishes the cache-lifetime cause; freeing caches
after every operation is not the proposed implementation.

[MPFR's memory-handling contract](https://www.mpfr.org/mpfr-current/mpfr.html#Memory-Handling)
asks users to free local caches before a thread terminates. Clearing only the
application thread's cache after a scheduler wait does not clear caches created
on other threads. A follow-up should choose the ownership of thread-exit cleanup
without introducing repeated cache clearing into numerical hot paths or changing
working precision. This remains a separate sanitizer blocker for the MPFR cases.
