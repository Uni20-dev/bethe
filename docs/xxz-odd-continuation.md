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
Negative-Delta odd rings remain excluded from the public ground-state API;
even rings and free ends are already supported.

**Development status: deferred until needed.** Retain this implementation
and its regression tests, but further odd-ring work and public integration
are not current priorities. The [restart checklist](xxz-negative.md#deferred-odd-ring-work-restart-checklist)
records the unresolved issues and the recommended first investigation.

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
  and the ratio of successive coupling steps;
- the energy is no greater than either of two independent trial-state
  Rayleigh quotients, with the arithmetic allowance described below.

The energy-change bound follows from `||dH/dDelta||<=N/4` and is necessary for a
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

### Independent variational rejection

`OddSectorVariationalBounds<Real>` in
[xxz_odd_bounds.hpp](../include/bethe/xxz_odd_bounds.hpp) evaluates two
explicit trial states without using Bethe equations or continuation history.
For any nonzero state in the same magnetization sector, its Rayleigh quotient
is an **upper** bound on the true sector minimum. A candidate above that
bound cannot be the minimum, even if its equations and momentum are correct.
Being below the bound is necessary, not sufficient.

Write `R=sin(pi*M/N)/sin(pi/N)` and `c=-cos(pi/N)`. The free-fermion ground
sea at Delta=0, held fixed as the coupling changes, has expectation

```text
E_sea(Delta) = c*R + Delta*[N/4 - M + (M*M-R*R)/N].
```

The first term is the sum of occupied one-particle cosine energies. For the
second, the density is `rho=M/N`, the nearest-neighbor one-body correlation
has magnitude `R/N`, and Wick's theorem gives
`<n_j*n_(j+1)>=rho*rho-(R/N)^2`. This includes the appropriate odd-ring
free-sea shift, not an even-ring occupation formula.

The [projected spin helix](xxz-spin-helix.md) of pitch `pi+pi/N` has equal
amplitude magnitude in every configuration in the sector. Set
`A=M*(N-M)/(N-1)`. Counting unlike neighbors and their exchange phases gives

```text
E_helix(Delta) = c*A + Delta*(N/4-A).
```

The trial vector remains valid at **every** coupling, even though it is an
eigenvector only at its commensurate coupling `Delta=c` (apart from trivial
sectors). This comparison only requires matching magnetization, not momentum:
it bounds a minimum in the whole magnetization sector, not a momentum block.
At Delta=0 the free sea saturates the sector bound; at the all-phantom
collision the helix expectation is `N*Delta/4`. Neither statement alone
classifies the spectrum at other couplings.

Each converged continuation stage must satisfy

```text
E <= min(E_sea,E_helix) + 256*epsilon*(1+N+abs(E)+abs(min(E_sea,E_helix))).
```

The allowance is a native-precision numerical safeguard, not an
outward-rounded error certificate. Failure rejects and shortens the step,
without changing the requested coupling or loosening the equation tolerance.
`variational_rejections` counts stages rejected by this check;
`variational_upper_bound` always reports the bound at the **requested** Delta,
including on failed returns. The state and global-minimum caveats still apply.

Tests independently enumerate the free Slater determinant and helix vector,
apply the spin Hamiltonian directly, and compare both expectations through
N=9 in all supported precisions. Small-sector ED checks the bound direction.
At N=17, M=8, Delta=-0.97, the helix bound excludes the previously observed
wrong branch by more than 0.15, without relying on the preceding secant or ED.

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

The 17-site regression is deliberately discriminating: before the concavity
and variational guards, an adaptive step from Delta=-0.9156494140625 to -0.96014404296875
jumped to another branch. It passed the full residual, momentum, conditioning,
and Lipschitz checks, but reached approximately -4.02462082373559 at
Delta=-0.97 instead of the independent sector minimum -4.18131119490746.
The discrepancy occurred in both fp64 and long double; it was a branch-tracking
failure, not a reason to loosen the energy comparison.

The driver tests run in fp64, long double, and optional fp128, under GCC and
Clang. The optional reference generator reproduces all nine minimal-|Sz|
sparse-ED energies, with eigenpair residuals below 6e-13.

An independent [quantum-Wronskian diagnostic](xxz-wronskian.md) now checks
the generic-q coefficient identity, and projected spin-helix tests establish
nonzero eigenvectors at the phantom collision. The diagnostic remains
separate from this driver: neither result is a blanket admissibility policy.
The complementary [regular-state criterion](xxz-regularity.md) now tests
the nonzero factors and full Jacobian required by the Gaudin norm formula
without extracting roots. Infinite/colliding/string limits remain separate.

Passing these finite-size checks does not establish general sector-minimum
tracking, particularly at singular or root-of-unity configurations. Larger
polynomial systems may still become ill-conditioned. Public integration
needs an explicit physical-state policy: the
smallest-|Sz| sector is not always the global ground state on a negative-Delta
odd ring. The all-sector candidate comparison below addresses that selection
mechanism, but does not certify its inputs. Public integration and excitation
classification remain unfinished.

## Comparing sectors without assuming the answer

The internal `scan_odd_polynomial_sectors<Real>` in
[xxz_odd_sectors.hpp](../include/bethe/xxz_odd_sectors.hpp) follows every
spin-reversal-distinct branch, `M=0,...,floor(N/2)`, at the requested coupling.
It deliberately returns an `OddSectorScan`, not a public `GroundState`.
The entries are ordered by increasing M, with positive
`Sz=N/2-M`; the negative-Sz partners have the same energy.

This matters even on five sites. At Delta=-0.9, the polarized sector has
energy -1.125, below the smallest-|Sz| value of approximately -0.9898034892.
Selecting `Sz=1/2` unconditionally would return the wrong global candidate.
The scan compares all sectors instead of assuming a universal phase boundary
from this or any other finite-size example.

The result keeps three distinct questions visible:

1. `equations_complete`: did every continuation converge, with a finite
   energy? If not, `lowest_index` is absent. The successful sectors cannot
   determine a minimum while a failed sector might lie below them.
2. Each sector's optional `wronskian`: what does the independent
   [Wronskian test](xxz-wronskian.md) say? It is absent for failed
   continuations. `wronskians_consistent` summarizes these diagnostics but
   is not a physical-state certificate; root-of-unity and conditioning
   limitations still apply.
3. `lowest_index`: which of the **complete set of numerically followed
   branches** has the lowest reported energy? This does not prove that the
   branches are physical or that each is its sector's minimum.

Each converged entry also carries a `regularity` result from the
[Gaudin-based regular-state test](xxz-regularity.md).
`regular_states_complete` summarizes those results independently of
`wronskians_consistent`. In particular, a regular vacuum at a root of unity
can pass the former while the generic Wronskian is inconclusive. Exceptional
states require a limiting construction, not automatic rejection.

The [explicit helix diagnostic](xxz-spin-helix.md) now recognizes the
all-phantom collision when the regular check is unresolved.
`state_checks_complete` reports whether every converged sector passes the
regular check, matches that helix, or has a resolved numerical phantom
witness; none constitutes a sector-minimum
test. Coupling mismatches are recorded, not rounded away, and failed
continuations are still never filled in or omitted.

The [numerical phantom-witness diagnostic](xxz-phantom-check.md) combines
mixed-root reduction, root recovery, and the coordinate dressing map. It
requires a dressed amplitude resolved against propagated root uncertainty
and an arithmetic allowance. Failed attempts and explicit work limits
remain visible in `phantom_lifts` and `phantom_work_limited`; neither means
the state vanishes. The map has a kernel, so reduction alone is never
enough. These numerical witnesses are not rigorous interval certificates.

The minimum is chosen by native-precision comparison, without rounding
energies or discarding a sector because its Wronskian is inconclusive.
An exactly equal numerical value keeps the earlier entry. `nearby_indices`
lists values within the reporting band

```text
128*epsilon*max(N,abs(minimum_energy)).
```

This band neither changes the selected minimum nor establishes degeneracy
or an energy-error bound. At the phantom collision, all folded sectors in
the small-ring tests fall within it; the scan retains their original energies.

The Newton budget applies **per sector**, as in the existing sector APIs.
The scan's `iterations` sums all accepted updates, including those in rejected
continuation stages. It finishes the other sectors after a numerical failure
so the caller receives their diagnostics, but never manufactures a minimum
from the incomplete set. Changing only the Wronskian tolerance cannot change
the continuation energies, iteration counts, or selected index. The same
is true of changing only the phantom-witness options (the fifth argument).

Tests compare every folded sector and the selected minimum with spin-basis
ED through N=9, at Delta=0 and four negative couplings. Analytic N=5
comparisons use the selected native precision, including near Delta=-1.
Budget tests cover incomplete scans, and root-of-unity tests ensure that
Wronskian failure is not mislabeled as failed continuation or physical
inadmissibility.

All folded sectors at N=13 and N=17 are also checked against independent
sparse-ED references at Delta=-0.7 and -0.99. These cases select smallest
|Sz| and fully polarized sectors respectively. Regenerate the references with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 tests/reference_xxz_odd_ed.py \
    --all-sectors --sites 13 17 --deltas -0.7 -0.99
```

The largest basis in this scan audit is 24310 states; the observed eigenpair
residuals are below 8e-13. The existing default reference-generator command
still checks minimal-|Sz| sectors through N=21. NumPy/SciPy remain optional
maintainer tools, not solver or normal-test dependencies.
