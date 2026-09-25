# Periodic XXX chains: states and conventions

[Back to the overview](../README.md)

This guide connects a numerical result to the periodic spin-1/2 Heisenberg
model it represents. Start with the Hamiltonian normalization, then follow
the quantum numbers through to rapidities, energy, and momentum. For a first
command-line run, see [the CLI guide](command-line.md).

## Hamiltonian and normalization

The Hamiltonian is

```math
H=\sum_{i=0}^{N-1}\mathbf S_i\cdot\mathbf S_{(i+1)\bmod N},
\qquad J=1,\quad h=0.
```

This is the spin-1/2, not Pauli-matrix, normalization. For N=2 the periodic sum
counts the bond twice, giving E=-3/2. For N=4, E=-2; for N=6,
E=-(2+sqrt(13))/2. Complex roots and infinite-root SU(2) descendants are not
supported by this solver. The API here is periodic; [free-end XXX](open-chains.md) and
[finite-size XXZ](xxz.md) have separate APIs. Boundary fields and other
boundary conditions remain unsupported.

## Calculate a state from C++

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

## Rapidity and quantum-number conventions

For M finite real roots, the conventions are:

```math
\begin{aligned}
\phi(z)&=2\arctan z,\\
F_i&=N\phi(z_i)-2\pi I_i-\sum_{j\ne i}\phi\!\left(\frac{z_i-z_j}{2}\right),\\
E&=\frac N4-\sum_i\frac{2}{1+z_i^2},\\
P&=\pi M-\frac{2\pi}{N}\sum_i I_i\pmod{2\pi}.
\end{aligned}
```

`solve_real` uses $S^z =N /2-M$ and the conventional all-1-string quantum-number
window: $M \le N /2$, $\lvert 2I_{i} \rvert \le N -M -1$, and $2I_{i}$ has the parity of $N -M -1$.
This explicitly supported window is not a classification of every real-root
solution. Quantum numbers and $S^z$ use `uni20::half_int`; conversion to the
selected real type uses the doubled integer directly, never `to_double()`.
The empty quantum-number set gives the fully polarized state.

## Sectors and state labels

For sector minima, $M =N /2-\lvert S^z \rvert$. Even chains occupy the consecutive numbers
$I_{i} =i -(M -1)/2$, with zero-based $i$. Odd chains use $I_{i} =i -M /2$, one of two
reflection-related minima when $M >0$; negating and reversing that sequence
selects the other momentum. Negative $S^z$ uses the all-down reference vacuum
and sets `spin_reversed=true`. Its roots count up spins rather than down spins.

`RealState<Real>` contains the roots, quantum numbers, `sz`, energy, diagnostics,
and an exact integer `momentum_index` with `P=2*pi*momentum_index/N` in $[0,2\,\pi)$.
The integer index is derived modulo $N$, without rounding floating-point phases.
`GroundState<Real>` remains an alias for compatibility. Results with
`converged=false` are numerical iterates, not established eigenstates; their
momentum labels specify the requested state.

## Iteration and convergence

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

## Validation

Tests use exact small-chain energies and roots, an independent evaluation of
the equations, a published N=16 value, reflection symmetry, iteration-budget
checks, invalid input, and precision-preserving CLI output. High-precision
tests use an irrational exact energy to detect accidental double narrowing.
The suite also compares independently converged results across precisions.
Independent bit-basis exact diagonalization checks every sector minimum for
N=2..9 and the odd-chain one-spinon energies and momenta through N=9. The
momentum check adds a multiple of $(T +T ^{-1})/2$ to the Hamiltonian. This oracle
uses double precision; separate irrational analytic references test native
long-double and binary128 accuracy.

Next: [enumerate real-root excitations](excitations.md) or
[interpret the one-spinon branch](spinons.md). Methods and attribution are in
[References and provenance](../CITATIONS.md).
