# Building Bethe and using Uni20

[Back to the overview](../README.md)

Start with the pinned dependency for a reproducible build. Use the sibling
checkout when developing both projects, and enable binary128 only when you
need the extra precision.

CMake 3.28+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
and Uni20's numerical dependencies are required. Uni20 is pinned to a tested
commit and fetched automatically unless a parent already supplies `uni20_core`.
The current pin is `04163b9b`, including Uni20 PR #99's sanitizer policy and
PR #98's notification-based scheduler waits, alongside recoverable square
and QR least-squares solves, reusable LU/log-determinants, binary80 provider
integration, explicit-precision MPFR/MPC scalars and conversions,
uniform scalar/backend configuration, typed run metadata, runtime-precision
data tables, and token-preserving help.
Local overrides must provide these APIs too; configuration rejects checkouts
without the recoverable-solve/least-squares/log-determinant headers and scalar-capability API.
Application builds also enable Uni20's
optional CLI11 dependency; library-only builds do not require it.
Scheduler-backed state scans use `uni20_async`. The parallel frontends and
parallel regression tests additionally link oneTBB, supplied by Uni20.

## Build the pinned version

Create a build directory wherever convenient, change into it, and give CMake
the path to the source checkout. Replace `/path/to/bethe` below with that path;
it can also be relative to the build directory (for example, `../bethe`).

```sh
mkdir bethe-build
cd bethe-build
cmake /path/to/bethe -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
ctest --output-on-failure
```

The build directory may be inside or outside the source checkout. Equivalently,
configure from any working directory with explicit source and build paths:

```sh
cmake -S /path/to/bethe -B /path/to/bethe-build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/bethe-build --parallel
ctest --test-dir /path/to/bethe-build --output-on-failure
```

Keep generated files separate from the source. When sharing a checkout between
hosts, use machine-local build storage and do not share CMake caches.

## Run the programs

From the build directory, run an executable directly:

```sh
./bethe-xxx-pbc 4
```

The model guides and tutorials use bare command names such as `bethe-xxx-pbc`.
These work from any calculation directory when the executables are on `PATH`,
whether supplied by an installation or your own build. To use the current
build for the rest of this shell session:

```sh
export PATH="$PWD:$PATH"
bethe-xxx-pbc 4
```

Run the export only while in the build directory; `$PWD` supplies its absolute
path, so subsequent directory changes are safe. Alternatively, prepend `./` to
the tutorial commands while working in that directory, or use an absolute
executable path. For multi-configuration generators, use the directory containing
the chosen configuration's executables instead.

Solver commands do not require the source checkout as their working directory.
Relative output paths, such as `--csv spectrum.csv`, refer to the directory
where you run the command. Only the optional repository plotting and maintenance
scripts need access to the source files.

## Install the programs

After building, install into a prefix of your choice. From the build directory:

```sh
cmake --install . --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
bethe-xxx-pbc 4
```

No administrator privileges are needed for this user-owned prefix. To select
the prefix at configure time instead, pass `-DCMAKE_INSTALL_PREFIX=/path/to/prefix`
to CMake and subsequently use `cmake --install .`. Multi-configuration builds
also need `--config Release` (or the configuration you built).

All enabled frontends are installed to `bin/`, with `COPYING` under
`share/doc/bethe/`. The directories follow CMake's `GNUInstallDirs` settings,
including `CMAKE_INSTALL_BINDIR` and `CMAKE_INSTALL_DOCDIR`. Prefer relative
directory overrides so `--prefix` and package staging continue to work; see
[CMake's install reference](https://cmake.org/cmake/help/latest/command/install.html).
`--component BetheRuntime` selects just this payload; on Unix, package builders
can stage it with `DESTDIR=/path/to/stage cmake --install . --prefix /usr`.

`BETHE_INSTALL_APPS` defaults to `ON` for standalone builds and `OFF` when
embedded. Set it to `OFF` to omit the install rules. No programs are installed
with `BETHE_BUILD_APPS=OFF`; tests, dependency development files and source
scripts are never part of the Bethe runtime component.

This is a normal dynamically linked installation, not a self-contained binary
bundle. Shared dependencies selected by the build (for example BLAS/LAPACK or
MPLAPACK) must remain available to the runtime loader. For dependencies in a
nonstandard installed location, configure an appropriate `CMAKE_INSTALL_RPATH`
or use your platform's loader configuration. Do not point an install RPATH at
temporary dependency build directories that you intend to remove. Bethe does
not preserve build-tree RPATHs. One exception is source-built oneTBB: Bethe
installs its shared runtime under `CMAKE_INSTALL_LIBDIR`, its license under
`share/doc/bethe/`, and a relative runtime search path in the parallel
frontends. A system-provided oneTBB remains a system dependency. This keeps
the default source build relocatable without installing dependency headers
or development libraries.

Installation currently covers the **programs**, not a `find_package(Bethe)`
C++ package. Uni20 does not yet export an installable package for Bethe's
transitive targets. Use [embedding](#embed-the-library) for library consumers;
the runtime installer does not publish an incomplete header-only package.

With tests enabled and Python available, `ctest -R '^Install.runtime$'` checks
the installed payload, relocation, command help and small calculations in
each enabled precision. It also tests component installation and Unix staging,
using disposable prefixes rather than your configured install prefix. This
test requires all frontend targets to have been built and relative install
directories contained within the prefix; it is not registered for absolute
destinations or parent-directory traversal overrides.

## Use a development checkout

To use a local Uni20 checkout instead of the pinned version, configure from
your build directory with absolute paths to both source checkouts:

```sh
cmake /path/to/bethe -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFETCHCONTENT_SOURCE_DIR_UNI20=/path/to/uni20
cmake --build . --parallel
ctest --output-on-failure
```

To return to the pinned revision, use a fresh build directory without the
source override. Uni20 can use installed
dependencies or fetch compatible versions of fmt, oneTBB, and mdspan; it also
configures BLAS/LAPACK. This is not yet a component-only configure, though
default builds compile only the components required by our targets.

## LTO and linker settings

Bethe and the pinned Uni20 default dependency LTO/IPO to off when embedded.
An existing cache or explicit caller setting is preserved; use
`-DUNI20_ENABLE_LTO=OFF` in an older build tree if needed. Enabling it requires
coordinating LTO and linker settings across both projects.

## Enable binary128

For binary128, use a separate build directory and add
`-DUNI20_ENABLE_FLOAT128=ON -DUNI20_BACKEND_MPLAPACK=ON`.
To use an installed MPLAPACK 3.0+ binary128
package, add `-DUNI20_USE_SYSTEM_MPLAPACK=ON` and
`-Dmplapack_DIR=/path/to/lib/cmake/mplapack`; otherwise Uni20 can fetch it.
See [Uni20's scalar and provider setup](https://github.com/Uni20-dev/uni20/blob/717c798c047c9e844535d9b21a09c152b49ea9af/docs/linalg/scalar_configuration.md).

`UNI20_USE_SYSTEM_MPLAPACK=OFF` forces a v3.0.0 source fetch. Alternatively,
build MPLAPACK v3.0.0 separately and pass its **build directory** as
`mplapack_DIR`, with `UNI20_USE_SYSTEM_MPLAPACK=ON`. Its build-tree CMake package is supported;
no system-wide installation is needed. MPLAPACK must contain components for all
enabled scalar families: with a binary128-only package, also set
`-DUNI20_ENABLE_FLOAT80=OFF`. This disables Uni20's optional alias/provider,
not Bethe's native `long-double` arithmetic or the `fp80` format check.
Use C++23 as described in Uni20's provider setup above. Keep build directories
on machine-local storage when the source checkout is shared between hosts.

Scalar availability and linear-algebra providers are separate: binary128
arithmetic can also be enabled with `UNI20_BACKEND_MPLAPACK=OFF`. Uni20 still
uses MPLAPACK's scalar headers/math, but does not register its matrix kernels.
`bethe-info` reports these capabilities separately. The scalar options accept
`ON`, `OFF`, or `AUTO`; `AUTO` detects without fetching. Use a fresh build tree
when migrating from the removed `UNI20_ENABLE_MPLAPACK`, `UNI20_ENABLE_MPFR`,
or `UNI20_ENABLE_MPC` options.

## Enable arbitrary precision

Arbitrary precision is optional. Uni20 enables its real and complex scalars
together, requiring GMP, MPFR >= 4.2, and MPC development packages. All modes of
`bethe-xxx-pbc`, `bethe-xxx-obc`, `bethe-richardson` and `bethe-central-spin`
work without MPLAPACK; Uni20 supplies CPU QR for the latter two.
In a separate build directory:

```sh
cmake /path/to/bethe -G Ninja -DCMAKE_BUILD_TYPE=Release -DUNI20_ENABLE_MPREAL=ON
cmake --build . --parallel
./bethe-xxx-pbc 32 --precision digits-100
./bethe-xxx-pbc 17 --spinons --precision bits-256
./bethe-xxx-obc 16 --sectors --precision bits-256
./bethe-richardson --levels 0,1,2,3 --pairs 2 --g 1 --precision digits-100
./bethe-central-spin --couplings 1,0.7,0.3 --field 1 --sz 0 --precision bits-256
./bethe-info
```

For finite-ring structure factors, finite-size Lieb–Liniger and Bose–Fermi, additionally
enable MPLAPACK matrix kernels. Unlike QR, Uni20's arbitrary-precision square
LU solves currently require this provider:

```sh
cmake /path/to/bethe -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DUNI20_ENABLE_MPREAL=ON -DUNI20_BACKEND_MPLAPACK=ON
cmake --build . --parallel
./bethe-xxx-structure-factor 32 --precision bits-256 --threads 4
./bethe-lieb-liniger-obc 2 --length 1 --c 1e-50 --precision digits-100
./bethe-info
```

Uni20 can fetch MPLAPACK 3.0+ or use an installed/build-tree package with
`-DUNI20_USE_SYSTEM_MPLAPACK=ON -Dmplapack_DIR=/path/to/mplapack`.
That package must include its **MPFR component**; a binary128-only build is
not sufficient. Native float80 kernels are also requested where that family is
enabled; use `-DUNI20_ENABLE_FLOAT80=OFF` with an MPFR-only package. Binary128
can be enabled alongside arbitrary precision with `-DUNI20_ENABLE_FLOAT128=ON`
on supported configurations.

This MPFR path covers periodic/open XXX real roots, sector/ground states, periodic one-spinon
branches, real-root excitation scans, prepared raising form factors and the
finite two-spinon DSF, finite-size periodic/hard-wall Lieb–Liniger states and scans,
plus Richardson and central-spin ground-state continuation.
It does **not** enable arbitrary precision throughout all models or frontends.
The XXX real-root and Gaudin solvers need just `UNI20_ENABLE_MPREAL=ON`;
Lieb–Liniger and form factors additionally need `UNI20_BACKEND_MPLAPACK=ON`. Enabling that provider
also lets Uni20 dispatch the Lieb–Liniger LU and Gaudin QR operations to MPLAPACK.
No ambient precision default determines the calculation's
working precision. At the provider boundary, Uni20 temporarily sets and restores
MPLAPACK's thread-local defaults around each synchronous call.
See [the explicit-precision C++ interface](xxx-structure-factor.md#arbitrary-precision).

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

The long group currently contains five binary128 cases: the interacting
sine-Gordon vacuum and excited level, and the Hubbard charge-continuum oracle,
native-default accuracy and convention checks. Previous local runs took roughly
one to fifteen minutes each. Their full numerical tolerances remain unchanged;
each has a 30-minute timeout. The corresponding fp64/long-double cases and
faster binary128 kernel, precision and failure-path checks stay in the routine
suite. New costly cases should be declared with `LONG_RUNNING_TESTS` alongside
their target's `bethe_discover_tests` call, using GoogleTest names.

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
Runtime-precision Lieb–Liniger and form-factor coverage belongs to
`ci-wide`; the scalar-only sanitizer configuration also checks that unavailable
provider requests fail clearly. These jobs do not claim macOS, Windows or GPU
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

## Embed the library

`BETHE_BUILD_APPS` and `BETHE_BUILD_TESTS` default to `ON` standalone and
`OFF` when embedded. No Python bindings or installed CMake package are provided
yet. Use `add_subdirectory`/FetchContent and link `bethe::bethe` in a parent.
The public numerical target does not link `uni20_cli` or CLI11. If a parent
already provides Uni20 and enables Bethe's applications, it must create
`uni20_cli` by enabling `UNI20_BUILD_CLI=ON` **before** adding Uni20.
An existing standalone cache with `UNI20_BUILD_CLI=OFF` likewise needs
`-DUNI20_BUILD_CLI=ON` when applications are enabled.

## Source layout and development

Set `BETHE_BUILD_EXAMPLES=ON` to build library composition examples in the
`examples/` subdirectory. The first example prepares one Lieb–Liniger background
asynchronously and submits both elementary branches before observing it; see
[shared backgrounds](async-execution.md#shared-backgrounds). Examples are checked
in CI and are not installed with the command-line tools.

### Where code belongs

| Location | Responsibility |
| --- | --- |
| `include/bethe/` | Scalar-templated numerical API, independent of CLI parsing and formatting |
| `include/bethe/detail/` | Shared numerical implementation helpers |
| `apps/` | Command-line frontends, model-specific arguments, formatting and report metadata |
| `tests/` | Numerical regressions, independent reference checks and frontend tests |
| `data/citations.json` | Literature metadata used to generate the C++ registry and bibliography |

Formatting stays in `apps/`, separate from the numerical API and any future
Python bindings. All executables link the private `bethe_cli` helper target;
numerical headers and data-output unit tests remain parser-independent.

The repository's `.clang-format` is copied from Uni20; use `clang-format -i`
on changed C++ files to apply the shared style.

### Linear solves and failure policy

We link `uni20_core` for scalar facilities, `uni20_common` for `half_int` and
presentation facilities, and `uni20_linalg` for square and least-squares solves.
Frontends additionally use `uni20_cli` for option parsing.

Square Newton and dressing-equation solves use Uni20's recoverable
`solve_inplace_with_info`, with LAPACK/MPLAPACK for compatible column-major
workspaces and native CPU fallback where needed. The row-major vector
adapter in `detail/newton.hpp` rearranges its owned coefficient copy in place
and exposes non-owning column-major views; the caller's coefficients are
preserved without an additional matrix allocation.

That adapter retains Bethe's explicit `64 * epsilon<Real>()` relative-pivot
threshold. Models already using Uni20 matrix workspaces retain the default
zero threshold and their own residual/condition checks. Neither policy is a
condition-number guarantee. Singular, rejected-pivot and nonfinite failures
return to the model's stalled/ill-conditioned or continuation-rejection path;
failed workspaces never become accepted iterates. Shape errors remain distinct
from numerical failures. No process-global error policy is changed.

Rectangular Newton corrections and continuation tangents use Uni20's
`least_squares_with_info` through read-only row-major views in `detail/newton.hpp`.
Normal dispatch selects LAPACK/MPLAPACK QR or Uni20's CPU Givens implementation;
native long double remains supported without MPLAPACK. Bethe retains its
`64 * epsilon` relative-diagonal rejection threshold, which is neither a
numerical-rank determination nor a condition estimate. No normal equations
are formed. A failed correction leaves the RHS unchanged; success replaces
it with the solution. An inconsistent least-squares fit can succeed: physical
residual checks and continuation acceptance still belong to the model.

The rectangular adapter also accepts an explicit finite `uni20::Precision`
for runtime-precision scalars. At this low-level boundary Uni20 converts inputs
at the requested working precision. The Richardson and central-spin model APIs
are stricter: inputs and tolerance must already have the options' finite
precision; mixed, exact and unset runtime values require explicit conversion.
Their sums, predictor/corrector state, backtracking limits and output diagnostics
retain that precision, including zero-budget and analytic cases. Numerical rank
via SVD, condition estimates and
refinement remain coordinated with
[Uni20's numerical infrastructure plan](https://github.com/Uni20-dev/uni20/blob/717c798c047c9e844535d9b21a09c152b49ea9af/docs/linalg/precision_numerics_plan.md).

### Shared iteration and continuation

`detail/newton_backtracking.hpp` shares trial construction, normalization
and damping. Models retain their pivot policy,
admissible domain, residual acceptance and iteration-counter meaning.
The string solver still requires strict decrease even below tolerance.

`detail/eigenvalue_continuation.hpp` shares the Richardson/central-spin adaptive
predictor-corrector controller and trust-box correction, including attempted-step
budgets and rejected-stage handling. Tangents, trust-radius selection, physical
energy bounds and final diagnostics remain model-specific. The tests exercise
these contracts in each enabled precision; continuation policy remains in
Bethe while square factorization and numerical solve diagnostics belong to Uni20.

### Frontend helpers

Use the shared lifecycle and output helpers; keep scientific arguments,
convergence criteria and missing-value policy in the model frontend.
Paths below are relative to `apps/`.

| Header | Responsibility |
| --- | --- |
| `cli-common.hpp`, `precision.hpp` | Common CLI utilities, checked precision requests and dispatch |
| `program-options.hpp` | Parse/information/error lifecycle, exact option adapters, Bethe identity and citations |
| `report-common.hpp`, `excitation-report.hpp` | Shared report rendering |
| `data-output-options.hpp` | Export flags for typed-table frontends |
| `run-metadata.hpp` | Uni20 run metadata projected into human overviews and legacy export keys |
| `result-output.hpp` | A frozen numerical summary shared across batch result tables |

Each frontend declares its own `TableDescriptor` list when registering output
options. The shared `table-descriptors.hpp` defines only the descriptor type and
generic validation: it has no executable-name registry or model-specific tables.
Declarations supply help, parser choices, screen flags, availability notes and
explicit primary-table eligibility. Result construction supplies the actual
available tables, which must match those declarations before exports are opened.

Keep display flags separate from production requirements. Use
`output.needs("--roots", roots)` when deciding whether root data is needed;
the unchanged `roots` flag controls display, while a named export can request
the same producer independently. Shared producers may construct sibling tables,
but only requested destinations receive them. Never manufacture an empty table
to satisfy an unavailable name. `Output.table_contracts` checks frontend help,
selection errors and auxiliary-file/screen separation.

The biquadratic frontend has private headers in `apps/biquadratic/`:

- `options.hpp`: option declarations and validation.
- `report.hpp`: common metadata, spin-content and real-root table writers.
- `real.hpp`, `analytic.hpp`, `qsystem.hpp`, `singlet.hpp`, `clusters.hpp`:
  individual calculation families.

These are frontend implementation details, not public numerical APIs.

### Citation maintenance

Literature metadata is centralized in `data/citations.json`. See
[maintaining citations](citations.md) to regenerate the C++ registry and the
bibliography in `CITATIONS.md`. Generated files are checked in, so Python is
needed only for regeneration and optional maintainer tests, not normal builds.

### Documentation equations

Use GitHub's backtick-protected inline math and fenced `math` display blocks.
In all math, use `\lt` and `\gt` instead of literal angle brackets: these can
survive Markdown processing but fail in GitHub's browser renderer. Use `\mathrm`
for named functions such as atan2; GitHub rejects `\operatorname` even though
MathJax itself supports it. End multiline display rows with `\\{}`: GitHub's Markdown
renderer adds an unwanted backslash to a bare `\\` at the end of a line.
The empty group preserves the TeX row break without relying on trailing spaces.

Check the conventions locally, or verify that every expression survives
GitHub's Markdown processing unchanged (requires authenticated `gh` and network
access; renders only, without publishing):

```sh
python3 scripts/check_markdown_math.py
python3 scripts/check_markdown_math.py --github
```

This checks the Markdown-to-TeX boundary and known browser restrictions, not
the full MathJax grammar or mathematical content. The rendering API does not
run GitHub's browser math component: confirm unfamiliar notation in a browser
too. Testing raw TeX alone does not catch Markdown escaping or GitHub's restrictions.

Next: [run a calculation and interpret its output](command-line.md).
