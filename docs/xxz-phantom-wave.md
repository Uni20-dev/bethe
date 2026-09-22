# Constructing phantom-dressed wavefunctions

[Mixed-root reduction](xxz-phantom.md) | [Spin helices](xxz-spin-helix.md) | [Odd-ring continuation](xxz-odd-continuation.md)

`PhantomDressing<Real>` in
[xxz_phantom_wave.hpp](../include/bethe/xxz_phantom_wave.hpp) implements the
coordinate-space map that adds p same-chirality phantom roots to a state
with r finite particles. It consumes a callback for the finite state's
amplitudes, rather than extracting roots or allocating a full spin basis.
This is currently an internal construction/validation tool, not a new
ground-state frontend or an automatic admissibility flag.
The [finite-root coordinate evaluator](xxz-coordinate-wave.md) can supply
that callback from momentum factors, without enumerating all permutations.

## The dressing map

Let `q=exp(i*k0)`, with `Delta=cos(k0)`, and let `M=r+p`. For ordered
occupied sites `J=(j_0,...,j_(M-1))`, select p indices
`b_0<...<b_(p-1)` to be phantom particles. The other sites form the
finite-particle configuration A. Our coordinate-ansatz deduction is

```text
(D_p Phi)(J) = sum_B q^[sum_h (j_(b_h)+2*(r+h-b_h))] * Phi(A).
```

Each phantom contributes its plane wave `q^j`. Moving it past a finite
particle contributes `q^2`; there are `r+h-b_h` finite particles to its
right. This generalizes the one-finite-particle formula in the
[reduction guide](xxz-phantom.md#an-independent-mixed-state-wavefunction).
It also gives the identity map for p=0 and the projected helix for r=0.

The physical starting point is the phantom/reduced-root construction of
[Popkov, Zhang and Klümper](../CITATIONS.md#popkov-2021), equation (15) and
its reduced twisted equations. The explicit subset map above and its
implementation are our deduction, not a quoted formula from that paper.

## What the Hamiltonian identity does—and does not—establish

At `q^(N-2*r)=1`, the map obeys

```text
H_periodic^(r+p) * D_p = D_p * H_twisted^r,
exp(i*phi) = q^(-2*p).
```

Both Hamiltonians include the same polarized energy `N*Delta/4`. Our
boundary convention assigns coefficient `exp(-i*phi)/2` to the hop from
site N-1 to site 0, and its conjugate to the reverse hop. It gives
`exp(i*N*k)=exp(i*phi)` for a single particle. Local terms in the
intertwining identity cancel; periodic closure requires the separate
commensurability condition.

Thus, **if the dressed vector is nonzero**, an eigenvector of the reduced
twisted Hamiltonian yields a periodic eigenvector at the same energy.
The nonzero qualification matters. At N=7, r=2, p=1 and
`q=exp(2*pi*i/3)`, take the nonzero reduced eigenvector
`Phi(a,b)=q^(a+b)`. Then

```text
(D_1 Phi)(j_0,j_1,j_2) = q^(j_0+j_1+j_2)*(1+q^2+q^4) = 0.
```

The map therefore has a genuine kernel, not merely a poorly conditioned
matrix. This example itself contains phantom roots and is outside the
finite regular-root family. It does **not** show that a regular reduced
Bethe state is annihilated; it shows why injectivity of the whole map cannot
be silently assumed. A suitable nonzero-lift criterion for the followed
regular branches remains a separate requirement.

## API, work limits, and normalization

```cpp
using Real = long double;
using Complex = std::complex<Real>;
using bethe::xxz::detail::PhantomDressing;

Complex q{-Real{1}/Real{2}, -std::sqrt(Real{3})/Real{2}};
PhantomDressing<Real> dressing(7, 2, 1, q); // N=7, r=2, p=1
std::array<std::size_t, 3> occupied{0, 2, 5};
auto psi = dressing.amplitude(occupied, finite_amplitude);
// finite_amplitude(span<const size_t>) supplies Phi on two ordered sites.
```

`for_each_term` exposes each finite subset and its complex coefficient
instead of summing amplitudes. Its spans refer to scratch storage and must
not be retained by the visitor. Coordinate validation happens before the
first callback. The map is defined even away from commensurability;
`commensurability_error()` records `|q^abs(N-2*r)-1|`, and no coupling or
phase is snapped to make it vanish.

Each amplitude has exactly `binomial(M,p)` terms. Construction checks this
integer count against `max_terms` (default 100000), including overflow,
before invoking any callback. A failed budget raises `length_error`; it
never returns a partial sum. Per-amplitude scratch memory is O(M), while
work is O(`binomial(M,p)*(M+p*log(N))`) plus the callback cost. This is not
a polynomial-cost method for large r and p simultaneously.

The amplitude sum uses native-precision compensated real and imaginary
accumulation. Nonfinite callback values, terms, or sums raise an error.
No wavefunction normalization is imposed, and a small output is not
automatically called zero or nonzero. Applications must compare with the
input normalization and actual numerical resolution.

## Independent tests

The tests run in fp64, long double, and fp128:

- Every column of the full Hamiltonian identity at `(N,r,p)` equal to
  `(5,1,1)`, `(7,2,1)`, `(9,3,1)`, `(9,2,2)`, `(8,1,2)`, and `(9,6,1)`,
  including even N and negative-Sz sectors, for both chiralities.
  Wrong-twist and noncommensurate controls fail that identity.
- The explicit nonzero-input/zero-output kernel example above, including
  direct verification that the source is an eigenstate.
- Continued minimal-|Sz| states at N=7,9,11,13 with two finite roots and
  respectively one, two, three, and four phantom roots. Only in the test,
  the finite quadratic is solved and the subset coordinate evaluator is
  checked against the explicit two-magnon formula. Both its twisted
  Hamiltonian residual and the nonzero dressed vector's periodic residual
  are checked directly at native precision.
- Vacuum/identity limits, exact callback counts, integer-count overflow,
  and rejection before callbacks for invalid coordinates or work budgets.
- General [root recovery](polynomial-roots.md) and dressing at
  `(N,r,p)=(9,3,1),(11,4,1),(11,3,2),(13,4,2),(13,3,3)`, with nonzero
  vectors and independent twisted/periodic Hamiltonian residuals.

The construction is available for arbitrary finite-state amplitude
callbacks, but these tests do not establish a nonzero lift for every
continued branch or prove sector minimality. The sector scan now uses a
separate [numerical witness check](xxz-phantom-check.md), with root-uncertainty
propagation and explicit budgets. Public solver domains are unchanged.
