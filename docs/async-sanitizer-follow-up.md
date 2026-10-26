# Async dependency sanitizer findings

The 26 October 2026 integration was checked against Uni20
`1e7c4a96b5cabf352cb1b5bb7efc2c2065498d26` and fetched oneTBB 2022.3.0.
Clang 19, Debug (`-O1 -g1`), address/undefined-behaviour sanitizers and leak
detection expose two dependency findings. They reproduce without Bethe headers,
numerical kernels or callbacks. They are unresolved; no suppression or disabled
sanitizer check is added to Bethe's CI.

## Uni20 debug DAG allocations

Constructing and destroying `uni20::async::Async<int> ready(7)` reports one
56-byte `NodeInfo` allocation on this x86-64 build. With `UNI20_DEBUG_DAG`
enabled, `EpochQueue::initialize_node` calls `NodeInfo::create`; the class
explicitly documents process-lifetime allocations that are never deallocated.
LeakSanitizer therefore reports these even after the Async value, buffers and
tasks are gone. This is distinct from an outstanding Bethe calculation or
retained result payload. The current Debug preset enables this diagnostic path;
normal Release builds disable it.

Uni20 needs an explicit diagnostic ownership policy compatible with leak
checking. Do not disable leak detection for the whole Bethe suite to hide it.

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
proving every possible cause. This needs a provider fix or supported provider
update, not a Bethe scheduler workaround.

## Reproduce

From a configured `ci-sanitizers` build directory:

```sh
cmake --build . --target bethe_async_dependency_probe -j2
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.Uni20DebugNodeLifetime
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
  ./tests/bethe_async_dependency_probe --gtest_filter=DependencyOnly.TbbSuspendResume
```

Both commands currently exit unsuccessfully. The probe is opt-in; the regular
Async tests already cover the affected paths and remain enabled in CI. After
upstream fixes, rerun both probes and the unchanged sanitizer suite before
treating the integration as ready to merge.
