# Using the C++ library

[Build and precision configuration](building.md) · [Model APIs](models.md) ·
[Async composition](async-execution.md) · [Complete examples](async-examples.md)

Start with a synchronous model call and inspect its numerical status. Use Async
when you have independent solves, a parameter grid, or several curves sharing
one prepared background. It does not parallelize coupled Newton iterations.

## Embed the library

`BETHE_BUILD_APPS` and `BETHE_BUILD_TESTS` default to `ON` standalone and
`OFF` when embedded. No Python bindings or installed CMake package are provided
yet. Use `add_subdirectory`/FetchContent and link `bethe::bethe` in a parent.
The public numerical target does not link `uni20_cli` or CLI11. If a parent
already provides Uni20 and enables Bethe's applications, it must create
`uni20_cli` by enabling `UNI20_BUILD_CLI=ON` **before** adding Uni20.
An existing standalone cache with `UNI20_BUILD_CLI=OFF` likewise needs
`-DUNI20_BUILD_CLI=ON` when applications are enabled.

For a source checkout supplied by the parent project:

```cmake
cmake_minimum_required(VERSION 3.28)
project(my_calculation LANGUAGES CXX)
set(BETHE_BUILD_APPS OFF CACHE BOOL "Build Bethe frontends")
set(BETHE_BUILD_TESTS OFF CACHE BOOL "Build Bethe tests")
add_subdirectory(/path/to/bethe bethe-build)
add_executable(my_calculation main.cpp)
target_link_libraries(my_calculation PRIVATE bethe::bethe)
```

Bethe supplies its C++23 requirement transitively. Configure numerical providers
and optional scalars before adding Bethe, using the [build guide](building.md).
The installed command-line package is distinct from this source embedding route.

## A synchronous calculation

```cpp
#include <bethe/heisenberg.hpp>
#include <uni20/core/scalar_io.hpp>
#include <iostream>

int main()
{
    auto const state = bethe::heisenberg::ground_state<double>(4);
    if (!state.converged) return 2;
    std::cout << uni20::format_real(state.energy) << '\n';
}
```

This is the same four-site, periodic spin-1/2 problem as the
[first command-line calculation](getting-started.md). The total energy is -2
for J=1. The result also owns rapidities, quantum numbers, momentum and numerical
diagnostics. The [XXX guide](xxx.md) defines its coordinate and normalization.
Model result types differ; read their status and optional-value conventions
before consuming energies, especially after continuation or refinement failure.

## Explicit runtime precision

In an MPFR-enabled build, the corresponding call is:

```cpp
bethe::heisenberg::SolverOptions<uni20::mpreal> options(
    uni20::Precision::bits(256));
auto const state = bethe::heisenberg::ground_state<uni20::mpreal>(4, options);
```

Use finite precision in options and construct continuous model inputs and custom
tolerances at that precision. Do not rely on MPFR's ambient defaults or create
high-precision parameters by first rounding decimal values to binary64.
`uni20::scalar_at_precision<Real>(value, precision)` is useful in generic code;
exact integer inputs avoid an initial floating-point rounding. Check each
model's API for its accepted parameter types and precision agreement rules.

Working bits, residual tolerance, quadrature/mesh controls and iteration budgets
are different controls. Inspect convergence before interpreting more printed
digits as accuracy. Independent MPFR requests may use different precisions on
one scheduler; each request owns its controls.

## Compose larger calculations

- [Thermodynamic spinon densities](xxx-spinon-density.md): physical W2/W4/W6
  weights, rapidity conventions and numerical diagnostics.
- [Direct spinon moments](xxx-spinon-moments.md): fixed-momentum and global
  observables, physical reference denominators and separate error controls.
- [Async execution](async-execution.md): state scans, bounded grids, shared
  preparation, lifetime, cancellation and completion contracts.
- [Observable grids](parallel-grids.md): a caller-thread collection example
  suitable for streaming and ordered reductions.
- [Complete C++ examples](async-examples.md): buildable parameter sweeps and
  shared-background calculations, including native-only demonstrations.
- [Execution capabilities](app-capabilities.md): which paths each app exposes.

Wait for the final result of the operation you submitted. A scheduler-wide drain
is not the completion boundary for a library call. The
[development guide](development.md) describes reusable implementation helpers;
[testing](testing.md) covers independent numerical and lifetime checks.
