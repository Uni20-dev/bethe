# Building and using Uni20

[Back to the overview](../README.md)

Start with the pinned dependency for a reproducible build. Use the sibling
checkout when developing both projects, and enable binary128 only when you
need the extra precision.

CMake 3.28+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
and Uni20's numerical dependencies are required. Uni20 is pinned to a tested
commit and fetched automatically unless a parent already supplies `uni20_core`.

## Build the pinned version

From the repository root:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The examples in these guides use `build/` relative to the repository root.
Substitute your actual build directory throughout. When the checkout is shared
between hosts, build on machine-local storage: our convention is a
`build_codex` symlink to local storage, with separate build trees beneath it,
for example `-B build_codex/release`. Do not share CMake caches between machines.

## Use a development checkout

To work with the sibling Uni20 checkout instead of the pinned version:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFETCHCONTENT_SOURCE_DIR_UNI20="$PWD/../uni20"
cmake --build build --parallel
ctest --test-dir build --output-on-failure
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
`-DUNI20_ENABLE_MPLAPACK=ON`. To use an installed MPLAPACK 3.0+ binary128
package, add `-DUNI20_USE_SYSTEM_MPLAPACK=ON` and
`-Dmplapack_DIR=/path/to/lib/cmake/mplapack`; otherwise Uni20 can fetch it.
See [Uni20's provider setup](https://github.com/Uni20-dev/uni20/blob/ef7e8b38f853a109e48c6b1762cab074f076639e/docs/linalg/mplapack_binary128.md).

`UNI20_USE_SYSTEM_MPLAPACK=OFF` forces a v3.0.0 source fetch. Alternatively,
build MPLAPACK v3.0.0 separately and pass its **build directory** as
`mplapack_DIR`, with `UNI20_USE_SYSTEM_MPLAPACK=ON`. Its build-tree CMake package is supported;
no system-wide installation is needed. Enable only the binary128 backend and
use C++23 as described in Uni20's provider setup above. Keep build directories
on machine-local storage when the source checkout is shared between hosts.

## Run and extend the tests

CTest runs the GoogleTest numerical cases, CLI checks, and optional Python
citation-maintenance checks. Tests reuse a parent-provided GoogleTest target,
prefer an installed GoogleTest 1.12+ package, or fetch v1.17.0. GoogleTest is
only required when `BETHE_BUILD_TESTS=ON`; library consumers do not inherit it.

Numerical cases are instantiated separately for fp64, native `long double`,
and optional fp128. Names include the precision, so failures can be isolated:

```sh
ctest --test-dir build -N
ctest --test-dir build -R 'Hubbard.*fp128' --output-on-failure
build/tests/bethe_heisenberg_tests --gtest_filter='*/long_double.NativePrecisionSix'
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

## Embed the library

`BETHE_BUILD_APPS` and `BETHE_BUILD_TESTS` default to `ON` standalone and
`OFF` when embedded. No Python bindings or installed CMake package are provided
yet. Use `add_subdirectory`/FetchContent and link `bethe::bethe` in a parent.

## Source layout and development

`include/bethe/` contains the scalar-templated library; `apps/` contains thin
command-line front ends; `tests/` contains the regression suite. We currently
link `uni20_core` for scalar facilities, `uni20_common` for `half_int` and
the CLI presentation layer, and `uni20_linalg` for the Hubbard Newton solves.
The latter uses Uni20's native-precision dense-solve dispatch, including its
generic CPU path where needed; no separate linear algebra implementation is
vendored here. Formatting stays in `apps/`, separate from the
numerical API and any future Python bindings.
Generic CLI, precision dispatch, and report rendering live in
`apps/cli-common.hpp`, `apps/report-common.hpp`, and `apps/excitation-report.hpp`;
model-specific arguments and report metadata stay in their respective front ends.
The repository's `.clang-format` is copied from Uni20; use `clang-format -i`
on changed C++ files to apply the shared style.

Literature metadata is centralized in `data/citations.json`; see
[maintaining citations](citations.md) to regenerate the C++ registry and the
bibliography in `CITATIONS.md`. Generated files are checked in, so Python is
needed only for regeneration and optional maintainer tests, not normal builds.

Next: [run a calculation and interpret its output](command-line.md).
