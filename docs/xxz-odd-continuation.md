# Following odd-ring XXZ states through root collisions

[Polynomial equations](xxz-polynomial.md) | [Negative anisotropy](xxz-negative.md) | [Model catalogue](models.md)

The internal driver `bethe::xxz::detail::continue_odd_polynomial<Real>` in
[xxz_odd_continuation.hpp](../include/bethe/xxz_odd_continuation.hpp) follows
the free-fermion ground-sea polynomial towards `-1<Delta<=0` on odd periodic
rings. It uses the collision-regular [coefficient equations](xxz-polynomial.md)
derived from [Caux's coordinate Bethe ansatz](../CITATIONS.md#caux-xxz-coordinate).
All arithmetic, including linear solves, remains in the selected scalar type:
fp64, native long double, or MPLAPACK fp128.

This is an **internal numerical continuation**, not a public ground-state
solver. `equations_converged` means the numerical checks below passed; it
does not certify a nonzero physical Bethe vector or prove sector minimality.
The public XXZ input range is unchanged.

## Coordinates that retain the momentum branch

Let `M=N/2-|Sz|`, folding spin reversal so `M<=floor(N/2)`. Choose one of
the two reflection-related free seas,

```text
z_j = tan(pi*(j-M/2)/N), j=0,...,M-1,
P_total = pi*M + pi*M/N (mod 2*pi).
```

The driver represents a monic polynomial in the affine coordinate x:

```text
z = o+s*x,
o = -tan(pi/(2*N)),
s = sqrt((1+Delta)/(1-Delta) + o^2),
Q(x) = x^M + c_(M-1)*x^(M-1) + ... + c_0.
```

The original polynomial class defaults to `o=0,s=1`; its extended constructor
accepts any finite center and positive finite scale. For the minus-scattering
recurrence, the transformed factors are

```text
A(x) = 1+Delta-(1-Delta)*o^2 + (-(1-Delta)*o*s-i*Delta*s)*x,
B(x) = (1-Delta)*o*s-i*Delta*s + (1-Delta)*s^2*x,
A-B*x = 1+Delta-(1-Delta)*(o+s*x)^2.
```

The driving factor is `1+i*o+i*s*x`. Thus the same divided-Horner construction
removes self-scattering before reduction modulo Q, including at collisions.

Define `t=s/(i-o)` and the reciprocal polynomial

```text
W(t) = 1+c_(M-1)*t+...+c_0*t^M.
```

The chosen center already carries the required total momentum:
the physical root polynomial has `R(i)=(i-o)^M*W(t)`. Enforcing `Im(W)=0`
therefore fixes `conj(R(i))/R(i)` without extracting roots. Crucially, this
is a **linear** coefficient constraint. Since `Im(t)` is nonzero, eliminate

```text
c_(M-1) = -sum_(j=0)^(M-2) c_j*Im(t^(M-j))/Im(t).
```

Energy and momentum evaluation use W and `t*W'`, avoiding large powers of
`(i-o)/s`. When the coupling changes, rescale
`c_j <- c_j*(s_old/s_new)^(M-j)` to preserve the physical roots.

## Adaptive continuation and acceptance

Newton solves the first `M-1` raw coefficient equations after eliminating
the last coefficient with the momentum constraint. The analytic Jacobian
includes that elimination. **All M residuals**, including the unused Newton
row, are checked for acceptance. Backtracking requires a decrease of the
normalized full residual and rejects nonfinite trials or observable poles.

A stage is accepted only when:

- the normalized full coefficient residual and the complex momentum-phase
  error are both at most the requested `residual_tolerance`;
- the largest relative coefficient correction is at most the square root
  of that tolerance (a correction check, not an energy-error estimate);
- the reduced Jacobian reciprocal condition estimate exceeds `64*epsilon`;
- the energy change obeys `|E_new-E_old|<=N*|Delta_new-Delta_old|/4`,
  with an additional `256*epsilon*(1+|E_old|)` roundoff allowance;
- after two accepted couplings, the energy lies below the extrapolation of
  the preceding energy secant, with a rounding allowance scaled by the energies
  and the ratio of successive coupling steps.

The last bound follows from `||dH/dDelta||<=N/4` and is necessary for a
continuously tracked eigenvalue. It can reject some branch jumps but cannot
certify admissibility. The secant test uses an additional property of a
**sector minimum**: as a minimum of Rayleigh quotients affine in Delta,
its energy is concave. For successively decreasing couplings, the next energy
must lie below the extrapolation of the preceding secant. This is necessary,
not sufficient, for remaining on that minimum. In particular, the driver is
not a generic excited-state continuation algorithm.

The reciprocal condition estimate is computed from
the matrix infinity norm and the inverse obtained with the same factorization
as the Newton correction. No normal equations are formed.

The coupling step adapts to successful stages and is halved after failures.
A stage allows at most 24 accepted Newton updates before retrying a smaller
coupling step. Every accepted update, including updates in discarded stages,
counts towards **one global** `max_iterations` budget. The `M=0,1` branches
are analytic and can converge with a zero budget. No tolerance is silently
relaxed when arithmetic precision prevents convergence.

The result returns the coefficients, center, scale, requested `delta`,
attempted `root_delta`, numerical status, iteration and rejected-stage counts,
and diagnostics. On failure, energy and residual are still evaluated at the
**requested** Hamiltonian, not a nearby accepted coupling. They are not
eigenstate energies. Unavailable observables are NaN or infinite diagnostics.
The condition and correction diagnostics refer to the last Newton
linearization, potentially at `root_delta`, rather than the requested Delta.
Invalid inputs throw; exhausted budgets and numerical failures return status.

## Independent checks and remaining work

The typed tests compare every folded sector at N=3,5,7,9 against dense
spin-basis diagonalization at four negative couplings. Joint energy/translation
spectra check the momentum branch. They exercise the collision at
`Delta=-cos(pi/N)` for N=5,7,9, and compare the exact five-site two-magnon
energy at native precision, including `Delta=-1+128*epsilon`.

Larger-chain tests use independent sparse spin-basis ground energies for
N=13,17,21, `|Sz|=1/2`, at Delta=-0.5,-0.97,-0.999. The optional
[reference generator](../tests/reference_xxz_odd_ed.py) uses NumPy/SciPy,
not Bethe equations; its largest basis has 352716 states. These fp64
references complement, rather than replace, the native-precision analytic
checks. Neither NumPy nor SciPy is a build or normal test dependency.

The 17-site regression is deliberately discriminating: without the concavity
check, an adaptive step from Delta=-0.9156494140625 to -0.96014404296875
jumped to another branch. It passed the full residual, momentum, conditioning,
and Lipschitz checks, but reached approximately -4.02462082373559 at
Delta=-0.97 instead of the independent sector minimum -4.18131119490746.
The discrepancy occurred in both fp64 and long double; it was a branch-tracking
failure, not a reason to loosen the energy comparison.

This checkpoint passes the full 433-test GCC 13 Debug suite with fp128 and
299-test Clang 20 Release suite without MPLAPACK. The new driver contributes
21 and 14 typed cases respectively. The optional reference generator also
reproduces all nine sparse-ED energies, with eigenpair residuals below 6e-13.

Passing these finite-size checks does not establish general sector-minimum
tracking, particularly at singular or root-of-unity configurations. Larger
polynomial systems may still become ill-conditioned. Public integration
needs an explicit physical-state policy and global sector selection: the
smallest-|Sz| sector is not always the global ground state on a negative-Delta
odd ring. These steps remain unfinished, as does excitation classification.
