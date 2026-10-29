# Testing and continuous integration

[Build configuration](building.md) · [Development guide](development.md)

## Run and extend the tests

CTest runs the GoogleTest numerical cases, CLI checks, and optional Python
citation-maintenance checks. Tests reuse a parent-provided GoogleTest target,
prefer an installed GoogleTest 1.12+ package, or fetch v1.17.0. GoogleTest is
only required when `BETHE_BUILD_TESTS=ON`; library consumers do not inherit it.

Numerical cases are instantiated separately for fp64, native `long double`,
and optional fp128. Names include the precision, so failures can be isolated.
From the build directory:

```sh
ctest -N
ctest -R 'Hubbard.*fp128' --output-on-failure
./tests/bethe_heisenberg_tests --gtest_filter='*/long_double.NativePrecisionSix'
```

The tests serve different purposes. Small independent exact-diagonalization
oracles check sectors, quantum numbers, and multiplicities in fp64. Analytic
energies, root equations, and Jacobians separately check the selected arithmetic
without narrowing to double. A successful fp64 ED comparison alone does not
establish fp128 accuracy.

New models should use the precision type list and comparison helpers in
`tests/test_support.hpp`.
Uni20's `EXPECT_FLOATING_EQ`/zero-ULP comparisons test exact scalar round trips,
including native fp80 where available. Solver and analytic checks instead use
native-precision `EXPECT_REAL_NEAR` (strict absolute error less than the supplied
tolerance), not GoogleTest's double-based `EXPECT_NEAR`. Keep fatal assertions
before indexing roots or dereferencing optional gaps. Use parameter traces
inside scans and keep expensive ED loops grouped into logical cases.

CTest gives each test one BLAS/OpenMP thread; it can still run independent
cases concurrently with `ctest --parallel`.

Expensive end-to-end numerical checks carry the `long-running` label. From a
build directory, choose routine checks, long checks alone, or the complete suite:

```sh
ctest --label-exclude '^long-running$' --output-on-failure
ctest --label-regex '^long-running$' --output-on-failure
ctest --output-on-failure
```

The long group includes expensive binary128 and MPFR continuum, thermal and
precision-convergence checks. The full Hubbard doped Fredholm reference retains
its 80/160/256-bit accuracy assertions there; an independent 80-bit reference,
mixed-precision Async checks and failure paths remain in routine validation.
Long checks normally have a 30-minute timeout; individual cases can declare a
larger budget. Lieb–Liniger thermodynamic CLI checks run separately for each
precision and for physical corrections, retaining every original assertion.
New costly cases should be declared with `LONG_RUNNING_TESTS` alongside their
target's `bethe_discover_tests` call, using GoogleTest names.

### Continuous integration

Pull requests and pushes to `main` run three Linux configurations, defined in
[`CMakePresets.json`](../CMakePresets.json):

| Preset | Coverage |
| --- | --- |
| `ci-native` | GCC 13 Release; fp64 and native `long double`; OpenBLAS; optional scalars disabled |
| `ci-sanitizers` | Clang 19; AddressSanitizer and UndefinedBehaviorSanitizer; native scalars and MPFR; CPU kernels without BLAS/MPLAPACK |
| `ci-wide` | GCC 13 Release; binary80, binary128 and MPFR; OpenBLAS and source-built, optimized MPLAPACK |

Each builds **all frontends and test targets** and runs the routine CTest suite,
including CLI/export checks and the installed-program relocation test.
Runtime-precision Lieb–Liniger, nested Newton and complex form-factor tests
run in both `ci-wide` and the MPFR CPU sanitizer configuration, exercising
MPLAPACK and CPU LU through normal dispatch. These jobs do not claim macOS, Windows or GPU
coverage. Documentation has its own workflow.

The sanitizer preset sets `UNI20_SANITIZE` as well as the application's compiler
flags, activating Uni20's oneTBB instrumentation boundary. The compiled oneTBB
implementation is excluded; Bethe and Uni20 remain instrumented. Leak detection
stays enabled. The preset enables `BETHE_TEST_UNI20_LSAN_SUPPRESSIONS`, which
selects `tests/sanitizers/lsan.supp` from the actual Uni20 dependency checkout
for every CTest case, including CLI subprocesses. This also works with local
source overrides and custom build directories. The upstream file covers
intentionally retained `NodeInfo` metadata and known MPFR internal cache frames;
ordinary owning scalar leaks remain detectable. Bethe keeps no separate copy
of these rules. See [the dependency sanitizer notes](async-sanitizer-follow-up.md).
Suppression statistics are quiet so successful commands retain their normal
stderr behavior; actual sanitizer errors still print and fail the process.
The preset clears `DEBUGINFOD_URLS` so symbolization uses locally available
debug information without waiting for remote downloads during CLI checks.
Remainder operations also use an MPFR integer pool whose private frame needs
matching debug symbols. The sanitizer CI job prepares these with Uni20's
`scripts/prepare-sanitizer-symbols.py` before running tests, verifies the loaded
library's build ID and records `sanitizer-symbols/symbols.json` in its artifacts.
For local sanitizer runs, follow Uni20's
[symbol preparation instructions](https://github.com/Uni20-dev/uni20/blob/main/docs/development/testing.md)
and set `ASAN_SYMBOLIZER_PATH` to the resulting wrapper. Downloads happen only
during preparation; neither the solvers nor sanitizer reports fetch symbols.

Weekly CI (Monday 02:00 Asia/Taipei) additionally runs `ci-long-running` on the
wide-precision build. A manual **C++ CI** run includes it by default, controlled
by the **Include long-running numerical validation** input. Run this group when
changing the relevant solvers, quadratures or precision machinery; ordinary
pushes and pull requests do not run it. Results are uploaded separately from
the routine suite. From a configured and built `ci-wide` tree:

```sh
ctest --preset ci-long-running
```

The workflow fetches Bethe's pinned Uni20, limits compilation and concurrent
tests to two jobs, and gives numerical kernels one BLAS/OpenMP thread.
Compiler outputs are cached; CMake/build directories are configured afresh.
JUnit results, capability reports and configure/build/test logs are uploaded
even on failure. MPLAPACK is built in Release mode without fast-math; the
first uncached wide-precision build is expected to take longer.

With the named compiler, CMake 3.28+, Ninja, Python and the relevant numerical
development packages installed, reproduce a configuration from the source root:

```sh
cmake --preset ci-native
cmake --build --preset ci-native
ctest --preset ci-native
```

Presets default to `build/<preset>` and do not need `ccache` locally. For a
shared source checkout, keep the build on local storage instead:

```sh
cmake --preset ci-native -B /path/to/local/bethe-ci-native
cmake --build /path/to/local/bethe-ci-native --parallel 2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ctest --test-dir /path/to/local/bethe-ci-native --parallel 2 \
    --label-exclude '^long-running$' --output-on-failure
```
