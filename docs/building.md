# Building Bethe and using Uni20

[Back to the overview](../README.md)

Start with the pinned dependency for a reproducible build. Use the sibling
checkout when developing both projects, and enable binary128 only when you
need the extra precision.

CMake 3.28+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
and Uni20's numerical dependencies are required. Uni20 is pinned to a tested
commit and fetched automatically unless a parent already supplies `uni20_core`.
The exact source revision is recorded in
[`CMakeLists.txt`](../CMakeLists.txt); `bethe-info` reports the revision actually
used by an installed build, including a local source override. Use the pinned
version unless developing both libraries together. It supplies the numerical,
Async, precision and presentation APIs used throughout Bethe; an older checkout
that passes the initial configuration checks may still lack newer APIs.
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
See [Uni20's latest scalar and provider setup](https://github.com/Uni20-dev/uni20/blob/main/docs/linalg/scalar_configuration.md).
For a revision-matched reference, read `docs/linalg/scalar_configuration.md`
in the Uni20 source tree used by your build.

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

All numerical frontends support explicit MPFR precision. Uni20 enables
its arbitrary-precision real and complex scalars together, requiring GMP,
MPFR >= 4.2, and MPC development packages. Set `UNI20_ENABLE_MPREAL=ON`;
MPLAPACK matrix kernels are optional. With `UNI20_BACKEND_MPLAPACK=OFF`,
Uni20 supplies CPU real/complex LU and QR where needed; normal dispatch prefers
an enabled LAPACK provider. See the [capability inventory](app-capabilities.md)
for each model's execution path and numerical limits.
In a separate build directory:

```sh
cmake /path/to/bethe -G Ninja -DCMAKE_BUILD_TYPE=Release -DUNI20_ENABLE_MPREAL=ON
cmake --build . --parallel
./bethe-xxx-pbc 32 --precision digits-100
./bethe-xxx-pbc 17 --spinons --precision bits-256
./bethe-xxx-obc 16 --sectors --precision bits-256
./bethe-richardson --levels 0,1,2,3 --pairs 2 --g 1 --precision digits-100
./bethe-central-spin --couplings 1,0.7,0.3 --field 1 --sz 0 --precision bits-256
./bethe-su3-dispersion --branch four-soliton --folded --precision digits-100
./bethe-sine-gordon-dispersion --p 0.4 --momentum 0.7 --precision digits-100
./bethe-info
```

The same build also supports Newton solves and complex form-factor determinants:

```sh
./bethe-xxx-structure-factor 32 --precision bits-256 --threads 4
./bethe-lieb-liniger-obc 2 --length 1 --c 1e-50 --precision digits-100
./bethe-hubbard-obc 5 --u -4 --precision bits-256
./bethe-xxz-pbc 16 --delta 0.5 --precision bits-256
./bethe-xxz-obc 16 --delta 3 --precision digits-100
./bethe-su3-pbc 24 --precision digits-100
./bethe-lieb-liniger-dispersion --c 4 --momentum 1 --precision digits-100
./bethe-info
```

For optional MPLAPACK matrix kernels, add `-DUNI20_BACKEND_MPLAPACK=ON`.
Uni20 can fetch MPLAPACK 3.0+ or use an installed/build-tree package with
`-DUNI20_USE_SYSTEM_MPLAPACK=ON -Dmplapack_DIR=/path/to/mplapack`.
That package must include its **MPFR component**; a binary128-only build is
not sufficient. Native float80 kernels are also requested where that family is
enabled; use `-DUNI20_ENABLE_FLOAT80=OFF` with an MPFR-only package. Binary128
can be enabled alongside arbitrary precision with `-DUNI20_ENABLE_FLOAT128=ON`
on supported configurations.

The selected precision applies within each model's existing physical domain
and state families. More bits do not add missing sectors, remove a quadrature
cutoff or guarantee convergence. Model-specific controls remain in the
[model guides](models.md#current-coverage).

CPU and MPLAPACK paths use the same explicit-precision interfaces. Enabling
MPLAPACK lets Uni20 dispatch supported LU and QR operations to that provider;
Bethe does not select a backend manually. This CPU coverage does not imply
that every Uni20 operation has an MPFR fallback: condition estimation, SVD
and eigensystems retain their provider requirements.
No ambient precision default determines the calculation's
working precision. At the provider boundary, Uni20 temporarily sets and restores
MPLAPACK's thread-local defaults around each synchronous call.
See [the explicit-precision C++ interface](xxx-structure-factor.md#arbitrary-precision).

## Library users and contributors

For parent-project CMake setup and a small synchronous calculation, see
[the C++ library guide](library.md). For source architecture, frontend helpers,
tests and CI, see [development](development.md) and [testing](testing.md).

<a id="run-and-extend-the-tests"></a>
[Run and extend the tests](testing.md#run-and-extend-the-tests).

<a id="continuous-integration"></a>
[Continuous integration](testing.md#continuous-integration).

<a id="embed-the-library"></a>
[Embed the library](library.md#embed-the-library).

<a id="source-layout-and-development"></a>
[Source layout and development](development.md#source-layout-and-development).

<a id="where-code-belongs"></a>
[Where code belongs](development.md#where-code-belongs).

<a id="linear-solves-and-failure-policy"></a>
[Linear solves and failure policy](development.md#linear-solves-and-failure-policy).

<a id="shared-iteration-and-continuation"></a>
[Shared iteration and continuation](development.md#shared-iteration-and-continuation).

<a id="frontend-helpers"></a>
[Frontend helpers](development.md#frontend-helpers).

<a id="citation-maintenance"></a>
[Citation maintenance](development.md#citation-maintenance).

<a id="documentation-equations"></a>
[Documentation equations](development.md#documentation-equations).
