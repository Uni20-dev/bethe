# Open XXX chains with free ends

[Back to the overview](../README.md)

Removing the periodic bond changes more than the total energy: translation
momentum is no longer a quantum number, and the Bethe equations include
reflected scattering. A separate front end and result type keep these
differences explicit. Precision and presentation controls are shared with the
[other programs](command-line.md).
For anisotropic interactions, see [free-end XXZ chains](xxz-open.md).

## Calculate a free-end state

`bethe-xxx-obc` has free ends with no boundary fields:

```sh
build/bethe-xxx-obc 16
build/bethe-xxx-obc 15 --sz 3/2
build/bethe-xxx-obc 16 --sectors
build/bethe-xxx-obc 4 --quantum-numbers 1,2 --roots --precision fp128
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

## Hamiltonian and reflected scattering

The free-end Hamiltonian and logarithmic equations use the same rapidity scale
as the periodic solver:

```text
H   = sum_(i=0)^(N-2) S_i . S_(i+1),
F_i = 2 N phi(z_i) - 2 pi I_i
      - sum_(j != i) [phi((z_i-z_j)/2) + phi((z_i+z_j)/2)],
phi(z) = 2 atan(z),
E   = (N-1)/4 - sum_i 2/(1+z_i^2).
```

## Allowed states and convergence

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

## Normalization checks and validation

For N=2 there is just one bond, giving E=-3/4, rather than the periodic
double-bond value -3/2. N=3 gives E=-1; N=4 gives E=-3/4-sqrt(3)/2. Tests compare
every sector minimum with independent exact diagonalization through N=9,
all supported real-root configurations through N=8, and the one-magnon standing
waves. Irrational analytic references test precision beyond double separately.

Next: [scan the open-chain real-root excitations](excitations.md).
Compare the [periodic XXX equations](xxx.md#rapidity-and-quantum-number-conventions)
and see [References and provenance](../CITATIONS.md) for the free-end derivation.
