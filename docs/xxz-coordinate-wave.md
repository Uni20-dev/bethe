# Finite-root coordinate wavefunctions

[Phantom dressing](xxz-phantom-wave.md) | [Polynomial equations](xxz-polynomial.md) | [Model catalogue](models.md)

`bethe::xxz::detail::CoordinateBetheWave<Real>` in
[xxz_coordinate_wave.hpp](../include/bethe/xxz_coordinate_wave.hpp) evaluates
an ordered spin configuration from finite, nonzero momentum factors
`v_j=exp(i*k_j)`. Complex factors are supported; they need not lie on the
unit circle. It is an internal wavefunction building block, not an
eigensolver or a new frontend.

## Coordinate convention

Starting from [Caux's coordinate ansatz](../CITATIONS.md#caux-xxz-coordinate),
we use its denominator-free numerator with

```text
F(i,j) = 1 + v_i*v_j - 2*Delta*v_i,
Psi(x) = sum_P sign(P) * product_(a<b) F(P_a,P_b)
                       * product_a v_(P_a)^x_a.
```

For two particles, the exchanged coefficient divided by the direct
coefficient is `-F(2,1)/F(1,2)`, as in the two-magnon scattering relation.
Sites are zero-based and strictly increasing. The public XXZ solvers'
scaled rapidities z convert to `v=-(1-i*z)/(1+i*z)` in all their supported
anisotropy regimes. This evaluator takes v explicitly and does not perform
root extraction or decide whether a root is singular.

Each unordered pair is divided by the positive scale
`max(1,abs(F(i,j)),abs(F(j,i)))`. This changes all configurations by the
same global normalization, not their relative amplitudes. It does not
divide by a possibly vanishing scattering factor. The returned vector is
not unit-normalized. At exceptional inputs the numerator itself can
vanish; a regularized limiting Bethe vector is a separate construction.

## Subsets instead of permutations

Let D(S) sum the permutations of the labels in S on the first `|S|`
occupied sites, with `D(empty)=1`. Select the last label j and set
`T=S\{j}`. The pair factors involving j no longer depend on the ordering
of T, so

```text
D(S) = sum_(j in S) D(T) * v_j^x_(|S|-1)
          * product_(i in T) F(i,j) * (-1)^(number of i in T with i>j).
```

The implementation uses the normalized pair factors above. This recurrence
is our regrouping of the coordinate ansatz, not a separate physical
assumption or a formula attributed to Caux. It needs O(`r^2*2^r`)
arithmetic and O(`2^r+r^2`) memory per amplitude, plus
O(`r^2*log(N)`) work to form the integer plane-wave powers. It avoids r!
enumeration but remains exponential in the number of particles r.

Construction checks `2^r` against `max_subsets` (default 65536) before
allocating root or pair storage. Excess work or integer-shift overflow
raises `length_error`, never a partial answer. The budget counts subsets,
not bytes or operations; r=0 still requires one subset.

## Using the result

```cpp
using Real = long double;
using Complex = std::complex<Real>;
using bethe::xxz::detail::CoordinateBetheWave;

std::array<Complex, 2> v{Complex{-1, 0}, Complex{0, 1}};
CoordinateBetheWave<Real> wave(8, Real{1}/Real{2}, v);
std::array<std::size_t, 2> occupied{0, 3};
auto result = wave.evaluate(occupied);
// result.value: the coordinate numerator with global pair normalization.
// result.absolute_term_sum: sum of absolute permutation terms.
```

The example illustrates the API, not an on-shell state. To supply the
[phantom map](xxz-phantom-wave.md), use a callback returning
`wave.evaluate(selected).value`. Neither component silently certifies a
nonzero eigenstate.

The absolute-term sum is accumulated through an unsigned version of the
same subset recurrence. Comparing `abs(value)` with it reveals cancellation
at that configuration. It is **not** a wavefunction norm, a rigorous error
bound, or a bound on root uncertainty. Real and imaginary sums use native
compensated accumulation; fp64, long double, and fp128 follow the same code.
Nonfinite factors, powers, terms, or sums raise `overflow_error`. Underflow
and severe cancellation remain possible: a floating-point zero is not a
proof that the exact state is zero. Pair scaling does not guarantee that
all intermediate plane-wave powers are representable.

## Validation and remaining work

Tests compare the recurrence with all permutations through six particles,
including complex off-shell inputs, and check the sign under root exchange.
A twelve-particle free-fermion configuration has an independently known
Fourier-determinant magnitude and absolute-term sum, testing 4096 subsets
in place of 479001600 permutations. Direct periodic spin-Hamiltonian
checks cover three and four roots at N=8,10 and
`Delta=-1/2,0,1/2,1,2`; perturbing a momentum fails the eigenstate check.
The two-finite-root phantom-dressing tests also use this evaluator, with
independent two-magnon amplitudes and twisted/periodic Hamiltonian checks.

Recovering general finite roots from the continued polynomial with adequate
accuracy remains a separate step. These tests do not establish a nonzero
lift of every mixed-phantom branch, prove sector minimality, or broaden the
public odd-ring ground-state domain.
