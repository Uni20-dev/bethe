# Non-Hermitian quantum-group XXZ chain

**Status: native positive finite-real-root library implemented for 0<Delta<1;
frontend, complex-root branches and root-of-unity representation accounting
remain follow-ups.** This is not the existing free-end XXZ model.

## Hamiltonian and source normalization

Use spin-half operators, J=1, N>=2, and opposite imaginary end fields:

```text
H = sum_{j=1}^{N-1} (Sx_j Sx_{j+1} + Sy_j Sy_{j+1} + Delta Sz_j Sz_{j+1})
    + i*sqrt(1-Delta^2)/2 * (Sz_1-Sz_N).
Delta = cos(gamma), 0<gamma<pi/2.
```

The source equations are [Gainutdinov–Hao–Nepomechie–Sommese (2015),
Eqs. (1.1)–(1.3)](https://arxiv.org/html/1505.02104). Our Hamiltonian is
their Pauli-normalized Hamiltonian divided by four, with the boundary sign
reversed. Spatial reflection reverses that sign without changing energies.
The fully polarized reference is E_F=(N-1)*Delta/4. There is no translation
momentum for this open chain.

With eta=i*gamma, the multiplicative Bethe equations are

```text
[sinh(lambda_i+i*gamma/2)/sinh(lambda_i-i*gamma/2)]^(2N)
 = product_{j!=i} sinh(lambda_i-lambda_j+i*gamma)/sinh(lambda_i-lambda_j-i*gamma)
                 * sinh(lambda_i+lambda_j+i*gamma)/sinh(lambda_i+lambda_j-i*gamma).
E = E_F - sum_i sin(gamma)^2/(cosh(2*lambda_i)-cos(gamma)).
```

## Implemented regular branch

`bethe/xxz_quantum_group_critical.hpp` uses
`z=tanh(lambda)/tan(gamma/2)`, with ordered positive roots and
`(1-Delta)*z^2 < 1+Delta`. Writing p=1+Delta, q=1-Delta, its logarithmic residual is

```text
F_i = 4*N*atan(z_i) - 2*pi*I_i
      - 2*sum_{j!=i} [atan(Delta*(z_i-z_j)/(p-q*z_i*z_j))
                     +atan(Delta*(z_i+z_j)/(p+q*z_i*z_j))].
```

Unlike the free-end model, there is no additional boundary phase in this
equation. A simultaneous iteration solves for `atan(z_i)` from the remaining
terms. The residual norm is `max(abs(F))/(2*N)`, not an energy-error estimate.
Energy is evaluated directly as
`E_F-sum((p-q*z^2)/(1+z^2))` to avoid unnecessary hyperbolic reconstruction.

`solve_real(N, Delta, labels, options)` accepts increasing positive integer
labels, M<=floor(N/2), I<=N-M and
`I < N-M+1-(N-2*M+2)*gamma/pi`. The last condition excludes a conservative
32-epsilon relative band at the infinity threshold. These necessary bounds
do not establish completeness or guarantee convergence for every label set.
`sea_state(N, Delta, ell, options)` chooses I=1,...,(N-ell)/2; ell must have
N's parity. For an odd chain length, explicitly supply odd ell.
No unrestricted ground-state or degeneracy claim
is inferred from this label selection alone.

```cpp
#include <bethe/xxz_quantum_group_critical.hpp>
namespace qg = bethe::xxz::quantum_group::critical;
auto state = qg::sea_state<long double>(32, 0.25L, 0);
if (state.converged) {
  auto energy = *state.energy;
}
```

Existing quantum-number, compensated-summation, scalar and solver-option
facilities are reused. Native fp64, long-double and fp128 follow the same
equations; no double fallback is used. Iteration work is O(M^2), storage O(M).
Failed solves retain their last root coordinates and residual but publish no
energy. Status distinguishes iteration exhaustion from a precision/branch
limit. An unrepresentable rapidity is also a precision failure. A zero update
budget still evaluates the seed, so an analytic one-root seed or vacuum can
converge without an update.

## Independent validation

`scripts/reference_xxz_nonhermitian.py` constructs the complex spin-basis
Hamiltonian in fixed magnetization sectors and uses a general eigensolver,
not a self-adjoint solver. NumPy is only a development-oracle dependency.
`--self-test` checks the two-site spectrum, the exact one-magnon spectrum,
and the defective two-site Delta=0 endpoint. The script also compares an
exploratory real-root iteration with independent eigenvalues.

The one-down-spin spectrum is E_F together with
`E_F-Delta+cos(pi*k/N)`, k=1,...,N-1. Small-chain direct calculations through
N=8 agree with the selected sea at Delta=0.25, 0.6 and 0.9. Native tests use
two-, three- and four-site analytic energies, independent five-/six-/eight-site
matrix references and direct substitution into the complex Bethe equations
through N=32. Thus both the imaginary boundary normalization and the absence
of the free-end boundary phase are tested.

## Why root-of-unity work is separate

At Delta=0 the even-chain sea reaches infinite rapidity. Already at N=2,
the Sz=0 Hamiltonian is nonzero but H^2=0: the two algebraically repeated
zero eigenvalues have only one eigenvector. An eigenvalue list cannot stand
in for Jordan chains. Generalized eigenvectors and complete strings require
additional treatment; see [Gainutdinov–Nepomechie (2016)](https://arxiv.org/html/1603.09249).
Both papers are also recorded in [the central citation registry](../CITATIONS.md).

Do not extrapolate the massive solver's generic representation multiplicities
to these parameters. Do not replace this chain with a Hermitian matrix having
the same real eigenvalues. Even in the implemented interval, regular-root
solutions do not by themselves specify the full root-of-unity spectrum.

The next stages are a boundary-specific frontend, admissible-label scans,
the Delta=0 free-fermion/Jordan benchmark, and additional complex-root branches.
CFT fitting must identify boundary sectors and distinguish c from an effective
central charge; a numerical Casimir coefficient is not automatically c.
RSOS restrictions and periodic loop realizations are separate representations,
not a boundary-field toggle on this open spin-chain solver.
