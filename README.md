# Bethe ansatz solvers

C++23 tools for finite-system Bethe ansatz calculations, complementing
[Uni20](https://github.com/Uni20-dev/uni20) and the
[Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

The first solver computes the real-root, zero-field ground state of the
even-length periodic spin-1/2 Heisenberg (XXX) antiferromagnet. It succeeds
MPToolkit's `misc/heisenberg-energy.cpp`; no legacy copy is kept here.

## Build

CMake 3.24+, a C++23 compiler supported by Uni20 (GCC 13+ or Clang 19+),
and Uni20's numerical dependencies are required. Uni20 is pinned to a tested
commit and fetched automatically unless a parent already supplies `uni20_core`.
For development with the sibling checkout:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFETCHCONTENT_SOURCE_DIR_UNI20="$PWD/../uni20"
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Omit the source override to fetch the pinned revision. Uni20 can use installed
dependencies or fetch compatible versions of fmt, oneTBB, and mdspan; it also
configures BLAS/LAPACK. This is not yet a component-only configure, though
default builds compile only the components required by our targets.

For binary128, use a separate build directory and add
`-DUNI20_ENABLE_MPLAPACK=ON`. To use an installed MPLAPACK 3.0+ binary128
package, add `-DUNI20_USE_SYSTEM_MPLAPACK=ON` and
`-Dmplapack_DIR=/path/to/lib/cmake/mplapack`; otherwise Uni20 can fetch it.
See [Uni20's provider setup](https://github.com/Uni20-dev/uni20/blob/a13cc911b58a670da48fba637ed1d13fd2d4f717/docs/linalg/mplapack_binary128.md).

`BETHE_BUILD_APPS` and `BETHE_BUILD_TESTS` default to `ON` standalone and
`OFF` when embedded. No Python bindings or installed CMake package are provided
yet. Use `add_subdirectory`/FetchContent and link `bethe::bethe` in a parent.

## Command line

```sh
build/heisenberg-energy 16
build/heisenberg-energy 64 --precision fp64 --tolerance 1e-12
build/heisenberg-energy 6 --precision long-double --roots
# In a binary128-enabled build:
build/heisenberg-energy 16 --precision fp128 --tolerance 1e-30
```

`--precision` selects `fp64`, `long-double` (the default), or optional `fp128`.
`long double` precision is platform-dependent; `fp128` uses Uni20's configured
MPLAPACK binary128 type. The fp128 CLI currently requires MPLAPACK's native
`_Float128/strfromf128` mode: other modes are rejected at compile time because
the pinned Uni20's scalar I/O narrows them to `long double`.
Calculations, tolerance parsing, and output retain the
selected precision; displayed digits are not a guarantee of energy accuracy.

The output includes energy, normalized equation residual, convergence status,
and update count. `--max-iterations` defaults to 10000; zero evaluates only the
initial zero-root guess. Exit status is 0 for convergence, 2 for budget
exhaustion, and 1 for invalid input or another error. A budget-exhausted result
is explicitly marked as an unconverged estimate. There is no silent precision
fallback. Run `--help` for the options.

## Library and conventions

```cpp
#include <bethe/heisenberg.hpp>

auto result = bethe::heisenberg::ground_state<long double>(16);
if (!result.converged) {
  // Inspect result.residual_norm and result.iterations, or increase the budget.
}
// result.energy includes N/4; result.rapidities contains the N/2 real roots.
```

The Hamiltonian is

```text
H = sum_(i=0)^(N-1) S_i . S_((i+1) mod N),     J = 1, h = 0.
```

This is the spin-1/2, not Pauli-matrix, normalization. For N=2 the periodic sum
counts the bond twice, giving E=-3/2. For N=4, E=-2; for N=6,
E=-(2+sqrt(13))/2. Odd lengths, other magnetization sectors, complex-root
states, and anisotropy are not supported by this first solver.

With zero-based i and r=N/2, the conventions are:

```text
I_i = i - (r-1)/2,       phi(z) = 2 atan(z),
F_i = N phi(z_i) - 2 pi I_i - sum_(j != i) phi((z_i-z_j)/2),
E   = N/4 - sum_i 2/(1+z_i^2).
```

The simultaneous fixed-point update starts at zero roots and follows Eq. (9)
of the [reference paper](https://arxiv.org/abs/cond-mat/9809163). It uses
compensated phase/energy sums, O(N^2) work per sweep, and O(N) storage.
Convergence means `max_i abs(F_i)/N <= residual_tolerance`; the default is
32 times the selected type's epsilon. This is a residual test, not a rigorous
energy-error bound. Unsupported inputs throw `std::invalid_argument`, nonfinite
arithmetic throws `std::runtime_error`, and budget exhaustion returns the last
roots with `converged=false`.

Tests use exact small-chain energies and roots, an independent evaluation of
the equations, a published N=16 value, reflection symmetry, iteration-budget
checks, invalid input, and precision-preserving CLI output. High-precision
tests use an irrational exact energy to detect accidental double narrowing.
The suite also compares independently converged results across precisions.

## Source and attribution

`include/bethe/` contains the scalar-templated library; `apps/` contains thin
command-line front ends; `tests/` contains the regression suite. We currently
need only `uni20_core`, for scalar traits, limits, math, and scalar I/O.

The successor code retains the original GPL-3.0-or-later licensing and
attribution; see [COPYING](COPYING) and [CITATIONS.md](CITATIONS.md).
