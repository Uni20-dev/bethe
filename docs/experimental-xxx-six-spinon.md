# Selected finite six-spinon states

The experimental `bethe::experimental::xxx::six_spinon` API solves selected
triplets of the periodic spin-1/2 XXX chain with two two-strings or one
three-string. It retains their actual finite-size deviations and computes
normalized raising form factors from a prepared singlet. This is a selected-state
solver, not a complete six-spinon spectrum or a thermodynamic extrapolation.
See the [experimental build contract](experimental.md).

We use spin operators sigma/2, exchange J=1, and standard ABA rapidities:

```math
H=\sum_{j=1}^{N}\mathbf S_j\cdot\mathbf S_{j+1},\qquad
E=\frac{N}{4}-\frac12\sum_a\frac{1}{\lambda_a^2+1/4}.
```

The singlet input to `PreparedSinglet` uses the supported real solver's
z=2 lambda convention; the preparation converts once. The returned excited
roots are lambda. Fourier raising weights are

```math
w=N\left|\langle e|S_0^+|0\rangle\right|^2.
```

There is no 2 pi multiplier or thermodynamic sector coefficient in this weight.
The [complex-root guide](experimental-xxx-complex-roots.md) describes its Gaudin
norm matrices, determinant evaluation and independent normalization controls.

At validated zero momentum transfer, total spin raising annihilates the
singlet. The form-factor result is `converged` with exactly zero `weight` and
absent `log_weight`. Family labels, original equations, precision and norms
are still checked before applying this selection rule; unresolved states
remain failures. The family passes its revalidated momentum label to the
shared complex-root helper, rather than identifying a zero from a small
computed weight.

## Requests and supported geometry

`Configuration` contains `sites`, `family`, six hole labels and the string
branch labels. Labels use `uni20::half_int`; `uni20::from_twice(j)` constructs
one from its doubled integer value. Invalid sizes, parity, repeated holes or
out-of-window labels throw before numerical work. N must be divisible by four
and at least 12. These requests select highest-weight S=1 states with N/2-1 roots.

| Family | Hole window, with unit spacing | String branches |
| --- | --- | --- |
| `Family::two_pairs` | Integers from -N/4 to N/4 | Two distinct values from -3/2, -1/2, 1/2, 3/2 |
| `Family::triple` | Half-integers from -N/4-1/2 to N/4+1/2 | One of -1, 0, 1 |

`labels(request)` returns the remaining Takahashi numbers, their string lengths
and the momentum index, with momentum 2 pi times that index divided by N.
Hole/branch order does not change the requested physical state; the result
retains the supplied request for accounting.

Ordinary geometry uses real roots and conjugate pairs. The two mapped N=12
three-string exceptions instead allow the complete root set to be self-conjugate
without each nominal string being self-conjugate. Their doubled Takahashi labels
are (3,5,2) and (-5,-3,-2). The mapping follows Deguchi and Giri's
[non-self-conjugate solutions](https://arxiv.org/abs/1408.7030), Table I.
It is a bounded mapping, not a general exceptional-string algorithm. Other
collapsed, singular or unresolved configurations can fail explicitly.

The Takahashi solution supplies a seed. Damped Newton iteration then solves the
full finite-root equations with an analytic Jacobian. Success additionally
requires the original multiplicative equations, momentum and recovered quantum
numbers to agree with the request. The latter check distinguishes different
states at the same momentum; its cut conventions follow the treatment of
[deformed strings](https://arxiv.org/abs/0707.2803) by Hagemans and Caux.

## A fresh calculation

From a build directory outside the source tree:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_EXAMPLES=ON
cmake --build . --target bethe_example_finite_six_spinon
./experimental/bethe_example_finite_six_spinon
```

The [complete example](../experimental/examples/xxx_six_spinon.cpp) prepares an
N=12 singlet and solves this two-pair request:

```cpp
namespace six = bethe::experimental::xxx::six_spinon;
six::Configuration request{
    12, six::Family::two_pairs, {-3, -2, -1, 0, 1, 3},
    {uni20::from_twice(-1), uni20::from_twice(3)}};
auto result = six::solve<double>(request);
if (result.converged()) {
    auto factor = six::raising_weight(prepared_singlet, result);
    if (factor.converged()) {
        auto weight = *factor.weight;
    }
}
```

Its momentum index is 3, energy is approximately 0.16684618680572643 and
Fourier raising weight is approximately 4.93041741948e-10. The example computes
them afresh; no catalogue, stored ground factors or private archive is required.

`SolveResult` carries the request, status, optional state, iteration count,
diagnostic text and attempted coordinates on failure. A successful state carries
the roots, parameterization, energy, momentum and measured residuals. The
raising evaluator checks the roots and family labels again and uses the prepared
ground's precision. Altering a result's cached status cannot make invalid roots
produce a weight. A failed state or form factor is an absent value, never a
zero-intensity spectral line.

## Precision and rejection

Native fp64, long double and enabled binary128 use `solve<Real>(request)`.
For MPFR/MPC, pass an explicit finite exemplar:

```cpp
using R = uni20::mpreal;
R zero(0, uni20::Precision::bits(256));
auto result = six::solve<R>(request, zero);
```

Prepare the singlet at that same precision. Exact, unset or mixed runtime
precision is rejected; no ambient MPFR default is used. Raising the requested
precision does not widen an already rounded ground state automatically.

The current unreduced representation rejects roots whose separation from poles
or coincidences is at most 64 sqrt(epsilon). It does not replace unresolved
roots by ideal strings. This gate and the determinant pivot diagnostic are
arithmetic rejection criteria, not relative-error estimates. In particular, a
known N=20 two-pair deviation near 2.435e-43 cannot be retained next to a height
of one half in native binary128; the corresponding request remains a failure.
The stabilized single-pair representation is a different family-specific method.

## Reproducing the independent controls

The compact C++ fixtures include both ordinary families at N=12, 16 and 20,
plus the two mapped N=12 exceptional states. The tests compare individual roots,
energies and form factors; check analytic Jacobians; and exercise wrong labels,
modified roots, unresolved native deviations and runtime precision validation.
MPFR checks use 128 and 256 bits. Selected success does not establish completeness
at any of those sizes.

The independent mpmath reference solves multiplicative equations and evaluates
direct products/determinants, separate from the C++ logarithmic form-factor
implementation. With mpmath >=1.2 installed, regenerate the fixtures:

```sh
python3 experimental/reference/xxx_six_spinon.py --digits 110 --verify
python3 experimental/reference/xxx_six_spinon.py --digits 130 --verify
```

Both compare the committed root/energy/weight values at at least 85 decimal
digits; the 130-digit run checks 90 digits. These are reference calculations,
not the achieved accuracy of a native-precision solve. Source provenance is the
finite-six research implementation at commit `1222b7ed`, with explicit precision,
failure records and shared preparation added for this experimental API.
