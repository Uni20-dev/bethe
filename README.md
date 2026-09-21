# Bethe ansatz solvers

C++23 tools for finite-system Bethe ansatz calculations, complementing
[Uni20](https://github.com/Uni20-dev/uni20) and the
[Matrix Product Toolkit](https://github.com/mptoolkit/mptoolkit).

The finite-size solvers cover periodic and free-end open spin-1/2 Heisenberg
(XXX) chains at zero field: ground states for even and odd lengths, lowest
energies in each magnetization sector, specified real-root states, and
energy-ordered real-root excitation scans at fixed total spin. The periodic
solver also provides the odd-chain one-spinon branch.
Analytic thermodynamic spinon dispersions are provided for XXX and
gapless XXZ. It succeeds MPToolkit's `misc/heisenberg-energy.cpp`; no legacy
copy is kept here.

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
Bethe defaults dependency LTO/IPO to off: the pinned Uni20's directory-local
LTO setting does not propagate to the application's final link step (notably
with Clang). An existing cache or explicit caller setting is preserved; use
`-DUNI20_ENABLE_LTO=OFF` in an older build tree if needed. Enabling it requires
coordinating LTO and linker settings across both projects.

For binary128, use a separate build directory and add
`-DUNI20_ENABLE_MPLAPACK=ON`. To use an installed MPLAPACK 3.0+ binary128
package, add `-DUNI20_USE_SYSTEM_MPLAPACK=ON` and
`-Dmplapack_DIR=/path/to/lib/cmake/mplapack`; otherwise Uni20 can fetch it.
See [Uni20's provider setup](https://github.com/Uni20-dev/uni20/blob/a13cc911b58a670da48fba637ed1d13fd2d4f717/docs/linalg/mplapack_binary128.md).

`UNI20_USE_SYSTEM_MPLAPACK=OFF` forces a v3.0.0 source fetch, but the pinned
Uni20 revision currently has an include-interface/export issue with that
embedded build path. A working alternative is to build MPLAPACK v3.0.0
separately and pass its **build directory** as `mplapack_DIR`, with
`UNI20_USE_SYSTEM_MPLAPACK=ON`. Its build-tree CMake package is supported;
no system-wide installation is needed. Enable only the binary128 backend and
use C++23 as described in Uni20's provider setup above. Keep build directories
on machine-local storage when the source checkout is shared between hosts.

`BETHE_BUILD_APPS` and `BETHE_BUILD_TESTS` default to `ON` standalone and
`OFF` when embedded. No Python bindings or installed CMake package are provided
yet. Use `add_subdirectory`/FetchContent and link `bethe::bethe` in a parent.

## Command line

`heisenberg-energy` is the periodic front end. For free-end OBC, use the separate
`heisenberg-open-energy` program described [below](#open-boundary-chains).

```sh
build/heisenberg-energy 16
build/heisenberg-energy 64 --precision fp64 --tolerance 1e-12
build/heisenberg-energy 6 --precision long-double --roots
build/heisenberg-energy 15 --sz 1/2
build/heisenberg-energy 16 --sectors
build/heisenberg-energy 65 --spinons
build/heisenberg-energy 32 --excitations 10 --spin 1
build/heisenberg-energy 5 --quantum-numbers -1,1 --roots
# In a binary128-enabled build:
build/heisenberg-energy 16 --precision fp128 --tolerance 1e-30
```

`--precision` selects `fp64` (the default), `long-double`, or optional `fp128`.
`long double` precision is platform-dependent; `fp128` uses Uni20's configured
MPLAPACK binary128 type. The fp128 CLI currently requires MPLAPACK's native
`_Float128/strfromf128` mode: other modes are rejected at compile time because
the pinned Uni20's scalar I/O narrows them to `long double`.
Calculations, tolerance parsing, and output retain the
selected precision; displayed digits are not a guarantee of energy accuracy.

The default is a ground state; odd lengths return one of the degenerate
`Sz=1/2` representatives. `--sz` accepts integers, fractions such as `-3/2`,
or half-integer decimals. `--sectors` reports one lowest-energy representative
in each sector, including spin-reversed partners. `--spinons` requires odd
`N>=3` and reports the one-spinon family, not the full `Sz=1/2` spectrum.
`--quantum-numbers` takes a strictly increasing comma-separated list;
an empty string selects the fully polarized state. `--excitations COUNT|all`
scans a restricted real-root family, described [below](#real-root-excitation-scans).
These five modes are mutually exclusive. `--roots` also prints the exact
Bethe quantum numbers.

Output uses Uni20's presentation layer on terminals: aligned fields and tables,
exact fractional quantum-number labels, and semantic convergence markers.
`--format auto` (the default) selects this report on a terminal and the existing
plain, whitespace-separated output when redirected. Use `--format pretty` to
save a formatted report or `--format plain` to request script-oriented output
explicitly. Both retain the selected type's full round-trip precision, including
binary128; no displayed values are narrowed to `double`.

The pretty report honors `UNI20_COLOR`, `NO_COLOR`, `UNI20_GLYPHS`, and
`UNI20_CHARSET`. For example, `UNI20_GLYPHS=ascii UNI20_COLOR=never` selects
ASCII table rules and status markers without ANSI color. Terminal width (or
`COLUMNS`) selects aligned tables or labeled records; numeric strings are never
split or truncated, so a single long value can still exceed a very narrow
terminal. Scans separate energies from convergence diagnostics, and root tables
identify their sector, hole, or excitation level.

Output includes energy, momentum, normalized equation residual, convergence
status, update count, and solver CPU time in seconds. CPU time measures process
CPU consumption during state construction and solving, not elapsed wall time;
for scans it covers the whole scan, including the ground reference and energy
ordering for excitations. Report formatting and output
are excluded. Both output formats include it, with a `# CPU time:` comment
before plain scan tables. Very short runs may report zero at the clock's
resolution; an unavailable or wrapped CPU clock is reported as `unavailable`.

`--max-iterations` defaults to 10000 **per state**;
zero evaluates only the initial zero-root guess. Exit status is 0 if all states
converged, 2 if any exhausted their budget, and 1 for invalid input or another
error. A budget-exhausted result is explicitly marked as an unconverged
estimate. There is no silent precision fallback. Run `--help` for the options.

## Library and conventions

```cpp
#include <bethe/heisenberg.hpp>

auto result = bethe::heisenberg::ground_state<long double>(16);
if (!result.converged) {
  // Inspect result.residual_norm and result.iterations, or increase the budget.
}
// result.energy includes N/4; result.rapidities contains the N/2 real roots.

using uni20::half_int;
auto sector = bethe::heisenberg::sector_ground_state<long double>(15, half_int::parse("3/2"));
auto sectors = bethe::heisenberg::sector_ground_states<long double>(16);
auto branch = bethe::heisenberg::one_spinon_branch<long double>(65);
bethe::heisenberg::QuantumNumbers numbers{half_int{-1}, half_int{1}};
auto specified = bethe::heisenberg::solve_real<long double>(5, numbers);
// Optional initial roots are the fourth solve_real argument, after options.
```

The Hamiltonian is

```text
H = sum_(i=0)^(N-1) S_i . S_((i+1) mod N),     J = 1, h = 0.
```

This is the spin-1/2, not Pauli-matrix, normalization. For N=2 the periodic sum
counts the bond twice, giving E=-3/2. For N=4, E=-2; for N=6,
E=-(2+sqrt(13))/2. Complex roots, infinite-root SU(2) descendants, finite-size
XXZ calculations, boundary fields, and other boundary conditions are not
supported yet. The API in this section is periodic; the open-chain API is separate.

For M finite real roots, the conventions are:

```text
phi(z) = 2 atan(z),
F_i = N phi(z_i) - 2 pi I_i - sum_(j != i) phi((z_i-z_j)/2),
E   = N/4 - sum_i 2/(1+z_i^2),
P   = pi M - (2 pi/N) sum_i I_i  (mod 2 pi).
```

`solve_real` uses `Sz=N/2-M` and the conventional all-1-string quantum-number
window: `M<=N/2`, `|2I_i|<=N-M-1`, and `2I_i` has the parity of `N-M-1`.
This explicitly supported window is not a classification of every real-root
solution. Quantum numbers and `Sz` use `uni20::half_int`; conversion to the
selected real type uses the doubled integer directly, never `to_double()`.
The empty quantum-number set gives the fully polarized state.

For sector minima, `M=N/2-|Sz|`. Even chains occupy the consecutive numbers
`I_i=i-(M-1)/2`, with zero-based `i`. Odd chains use `I_i=i-M/2`, one of two
reflection-related minima when `M>0`; negating and reversing that sequence
selects the other momentum. Negative `Sz` uses the all-down reference vacuum
and sets `spin_reversed=true`. Its roots count up spins rather than down spins.

`RealState<Real>` contains the roots, quantum numbers, `sz`, energy, diagnostics,
and an exact integer `momentum_index` with `P=2*pi*momentum_index/N` in `[0,2*pi)`.
The integer index is derived modulo `N`, without rounding floating-point phases.
`GroundState<Real>` remains an alias for compatibility. Results with
`converged=false` are numerical iterates, not established eigenstates; their
momentum labels specify the requested state.

The simultaneous fixed-point update starts at zero roots and follows Eq. (9)
of the [reference paper](https://arxiv.org/abs/cond-mat/9809163). It uses
compensated phase/energy sums, O(M^2) work per sweep, and O(M) storage per state.
The helpers returning all sectors or an entire branch retain every state's
roots, using O(N^2) storage. Use the single-state functions to stream large scans.
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
Independent bit-basis exact diagonalization checks every sector minimum for
N=2..9 and the odd-chain one-spinon energies and momenta through N=9. The
momentum check adds a multiple of `(T+T^-1)/2` to the Hamiltonian. This oracle
uses double precision; separate irrational analytic references test native
long-double and binary128 accuracy.

## Real-root excitation scans

Both front ends enumerate supported real-root highest-weight states at a
chosen **total spin S**, with `M=N/2-S` roots and `Sz=S`:

```sh
build/heisenberg-energy 64 --excitations 10 --spin 1
build/heisenberg-open-energy 64 --excitations 10 --spin 1
build/heisenberg-energy 64 --excitations all --spin 1
build/heisenberg-open-energy 64 --excitations all --spin 1
build/heisenberg-energy 65 --excitations 10 --spin 1/2 --precision long-double
build/heisenberg-open-energy 32 --excitations 5 --spin 2 --roots
```

`--spin` defaults to 1 for even N and 1/2 for odd N. It must be nonnegative,
no greater than N/2, and have the same integer/half-integer parity as N/2.
`--spin` and `--max-candidates` require `--excitations COUNT|all`; `--sz` remains
the separate sector-minimum mode, not a filter for this scan.

The scan returns up to COUNT lowest **converged multiplets in this family**,
including its sector minimum. Each multiplet is represented once; its `2S+1`
SU(2) partners are not listed separately. Distinct multiplets with equal
energies are retained, including periodic reflection/momentum partners. COUNT
can cut through such degeneracies. Results are sorted by computed energy;
exact ties use lexicographic Bethe quantum numbers. Near-degenerate ordering
can change with numerical precision, and the residual is not an energy-error
bound. For odd N and S=1/2 the family includes ground-state multiplets.

Use `--excitations all` to return every converged multiplet in the selected
real-root family, without needing its size in advance. This does not include
other total-spin sectors or unsupported string states, and still respects
`--max-candidates`. For N=64, S=1 it returns 528 multiplets if all converge.

Reports give S, absolute energy, `gap=E-E0` relative to the **global ground
state of the same finite chain**, quantum numbers, and convergence diagnostics.
Periodic reports additionally give lattice momentum; open reports do not.
Gaps retain the selected arithmetic precision, including fp128. Tiny negative
gaps due to roundoff are not clamped. Optional root blocks identify each level.
In plain output, S and I use exact fractions, the I column is comma-separated,
and `-` denotes an empty set of quantum numbers.

This is **not a complete low-energy spectrum**, even in the requested S
sector. It enumerates only the quantum-number windows of the existing
real-root solvers. For periodic even N, S=1 gives the conventional two-spinon
triplet family; further triplets and excited singlets require additional
families. In particular, for even N and S=0 the current window contains only
the ground-state configuration. Complex/string states and infinite-root
descendants are not solved by this mode.

There are `C(N-M,M)` candidates for either boundary: for example, N=64, S=1
has 528. The scan solves **every** candidate before returning the requested
lowest subset, not just the first COUNT configurations. The default
`--max-candidates 10000` rejects larger families before any solve; raise it
explicitly for larger jobs, including when using `all`. A numeric COUNT and
this limit must both be positive. The limit bounds the number of configurations,
not the cost of an individual solve. A bounded heap retains only O(COUNT*M+N)
data, including the global ground reference; the underlying O(M^2) work per
iteration is unchanged.
With `all`, roots for the entire converged family are retained.

Failed candidates are excluded from the energy-ordered list, but **all**
candidate failures count toward the scan status. Reports show total and
converged candidate counts, returned multiplets, and the first failed
configuration's diagnostics. An exhausted candidate makes the ordering
incomplete and returns exit status 2, even if all displayed states converged.
If only the ground reference fails, absolute-energy ordering can still be
complete within the family, but gaps are `unavailable` and the exit status
is also 2. No unconverged ground energy is used to manufacture gaps.

The independent small-chain tests compare the scans against exact
diagonalization through N=8, subtracting the Sz=S+1 spectrum from Sz=S to
check total-spin multiplicities, not just energy membership. They also check
momentum, exhaustive quantum-number coverage, bounded-prefix ordering,
omitted singlets, failure accounting, and an irrational gap beyond fp64.

Library usage:

```cpp
#include <bethe/heisenberg_excitations.hpp>

namespace xxx = bethe::heisenberg;
auto const spin = uni20::half_int{1};
auto const candidates = xxx::real_excitation_count(64, spin); // no solves
xxx::RealExcitationOptions selection{.count = 10, .max_candidates = 10000};
auto periodic = xxx::real_excitations<long double>(64, spin, selection);
auto open = xxx::open::real_excitations<long double>(64, spin, selection);
// Optional fourth argument: SolverOptions<Real>.
// Each scan.levels entry contains .state and .gap (std::optional<Real>).
// scan.family_converged() covers every candidate, not only retained levels.
// scan.converged() additionally requires the ground reference to converge.
```

`real_excitation_count` accepts an optional upper limit (default: maximum
`size_t`). Exceeding that limit, `max_candidates`, or representable count
throws `std::length_error`; invalid spin or scan options throw
`std::invalid_argument`. Nonfinite numerical failures propagate from the
underlying solver. Arithmetic, sorting, and gap subtraction never narrow
the selected real type to double.

## Single-spinon dispersion

For odd N, set `M=(N-1)/2`. The one-spinon family occupies M of the M+1 slots
`-M/2, -M/2+1, ..., M/2`, leaving a hole `I_h`. Use
`one_spinon_state<Real>(N, hole, options)` for one member, or
`one_spinon_branch<Real>(N, options)` for all `(N+1)/2` members in ascending k.
Each `SpinonState` contains its `RealState` in `.state`, the exact `.hole`, and:

```text
spinon_momentum       k = pi/2 - 2*pi*I_h/N,
bulk_subtracted_energy = E_N - N*e_inf,     e_inf = 1/4 - log(2).
```

The allowed finite-size k values run from `pi/(2N)` to `pi-pi/(2N)` in steps
of `2*pi/N`. This is a hole-based spinon convention, with
`P = pi*M + pi/2 - k (mod 2*pi)`, **not** momentum relative to the odd-chain
ground state. Both lattice P and spinon k are reported. The bulk-subtracted
energy retains finite-size corrections; it is not obtained by subtracting the
odd-chain ground-state energy, which itself contains a spinon.

For the analytic zero-field thermodynamic dispersion, include
`<bethe/spinon.hpp>` (also included by `<bethe/heisenberg.hpp>`):

```cpp
long double k = 1.0L;
auto xxx = bethe::heisenberg::spinon_energy(k);          // J=1
auto xxz = bethe::xxz::spinon_energy(k, 0.5L);           // Delta=0.5, J=1
auto scaled = bethe::xxz::spinon_energy(k, 0.5L, 2.0L); // J=2
auto e_inf = bethe::heisenberg::bulk_energy_density<long double>();
```

These implement `epsilon(k)=J*pi*sin(k)/2` for XXX and
`epsilon(k)=J*pi*sin(gamma)*sin(k)/(2*gamma)`, `gamma=acos(Delta)`, for gapless
XXZ. They require finite `k` in `[0,pi]`, positive finite J, and for XXZ
`-1<Delta<=1`. The isotropic endpoint uses the explicit XXX limit; k=0 and
k=pi return exactly zero. They are analytic functions, not a finite-size XXZ
solver. See [CITATIONS.md](CITATIONS.md) for the derivations.

## Open-boundary chains

`heisenberg-open-energy` has free ends with no boundary fields:

```sh
build/heisenberg-open-energy 16
build/heisenberg-open-energy 15 --sz 3/2
build/heisenberg-open-energy 16 --sectors
build/heisenberg-open-energy 4 --quantum-numbers 1,2 --roots --precision fp128
```

It shares precision selection, CPU timing, convergence diagnostics, and
pretty/plain reports with the periodic program, but has no `--spinons` or
boundary-selection switch. Open chains have no translation momentum, so no
momentum fields or columns are reported. The library likewise uses a separate
result type without momentum members:

```cpp
#include <bethe/heisenberg_open.hpp>

namespace obc = bethe::heisenberg::open;
auto ground = obc::ground_state<long double>(16);
auto sector = obc::sector_ground_state<long double>(15, uni20::half_int::parse("3/2"));
auto sectors = obc::sector_ground_states<long double>(16);
bethe::heisenberg::QuantumNumbers numbers{uni20::half_int{1}, uni20::half_int{2}};
auto state = obc::solve_real<long double>(4, numbers);
```

The free-end Hamiltonian and logarithmic equations use the same rapidity scale
as the periodic solver:

```text
H   = sum_(i=0)^(N-2) S_i . S_(i+1),
F_i = 2 N phi(z_i) - 2 pi I_i
      - sum_(j != i) [phi((z_i-z_j)/2) + phi((z_i+z_j)/2)],
phi(z) = 2 atan(z),
E   = (N-1)/4 - sum_i 2/(1+z_i^2).
```

Only the positive, finite real-root branch is represented: `M<=N/2`, with
distinct increasing integer labels `1<=I_i<=N-M`. Labels still use
`uni20::half_int`, but half-odd integers are rejected. The sector minimum fills
`I=1,...,M`, where `M=N/2-|Sz|`; negative sectors use spin reversal. The
ground-state wrapper selects Sz=0 for even N and Sz=1/2 for odd N.

The sum excludes **both** self-scattering terms, including the reflected root
of the same particle. The residual is `max|F_i|/(2N)`, not `max|F_i|/N`.
`SolverOptions<Real>` and update-budget semantics are shared with the periodic
solver. Initial guesses may be finite and nonnegative; the default is zero.
The physical roots are strictly positive. Zero-root Bethe vectors, complex
strings, and infinite-root descendants are outside the supported family;
this is not a complete-spectrum enumerator. N must be at least 2.

For N=2 there is just one bond, giving E=-3/4, rather than the periodic
double-bond value -3/2. N=3 gives E=-1; N=4 gives E=-3/4-sqrt(3)/2. Tests compare
every sector minimum with independent exact diagonalization through N=9,
all supported real-root configurations through N=8, and the one-magnon standing
waves. Irrational analytic references test precision beyond double separately.

## Source and attribution

`include/bethe/` contains the scalar-templated library; `apps/` contains thin
command-line front ends; `tests/` contains the regression suite. We currently
link `uni20_core` for scalar facilities and `uni20_common` for `half_int` and
the CLI presentation layer. Formatting stays in `apps/`, separate from the
numerical API and any future Python bindings.
The repository's `.clang-format` is copied from Uni20; use `clang-format -i`
on changed C++ files to apply the shared style.

The successor code retains the original GPL-3.0-or-later licensing and
attribution; see [COPYING](COPYING) and [CITATIONS.md](CITATIONS.md).
