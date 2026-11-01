# Experimental finite XXX four-spinon states

[Experimental components](experimental.md) · [Complex-root controls](experimental-xxx-complex-roots.md)

The finite-chain four-spinon contribution to the zero-field transverse structure
factor comes from triplet intermediate states. This experimental API solves a
selected nonsingular one-two-string triplet and evaluates its normalized raising
weight. It preserves the finite string deviation, including values too small to
change the rounded imaginary parts of the complex roots.

The current family requires N at least 8 and divisible by four. Four distinct
integer holes lie in [-N/4,N/4], and the two-string branch is -1, 0 or 1.
Centred/singular pairs and seeds outside the narrow-pair domain are explicitly
unsupported. This is neither a complete finite-chain spectrum nor the
[thermodynamic four-spinon integral](xxx-four-spinon-thermodynamics.md).

## A selected state

Enable `BETHE_BUILD_EXPERIMENTAL=ON`, link `bethe::experimental`, and include
`<bethe/experimental/xxx_four_spinon.hpp>`:

```cpp
namespace xxx = bethe::experimental::xxx;
auto gs = bethe::heisenberg::sector_ground_state<double>(12, uni20::half_int{0});
xxx::PreparedSinglet<double> ground(12, gs);
auto selected = xxx::solve_four_spinon<double>({12, {-3,-1,0,2}, 1});
if (ground.ready() && selected.converged()) {
    auto ff = xxx::four_spinon_raising_weight(ground, *selected.state);
    if (ff.converged()) {
        double gap = selected.state->energy() - *ground.energy();
        double weight = *ff.weight;
        // Add this discrete line at momentum index 9, gap and weight.
    }
}
```

The weight is the squared matrix element of the normalized Fourier raising
operator, without the spectral delta function's 2*pi factor. Its value for this
state is about 0.000012014419465322. At zero field the corresponding longitudinal
weight is half the raising weight. The singlet preparation can be reused for
other selected states; all per-state matrix work is separate.

The supplied example performs a fresh ground/state solve and prints the gap,
weight, deviation and original-equation residual. From your chosen build directory:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_EXAMPLES=ON
cmake --build . --target bethe-experimental-four-spinon-example
./experimental/bethe-experimental-four-spinon-example
```

The example uses 128-bit MPFR when configured, otherwise binary64. These headers
and the example are build-tree components and are not installed.

## Precision, states and failures

`solve_four_spinon<Real>(request, exemplar)` derives working precision from a
real exemplar. It is optional for native types; MPFR requires an explicit
finite-precision value, for example `mpreal(0, Precision::bits(256))`. The ground
solve and selected state must use the same precision. No ambient MPFR default
or backend is selected by this code.

Invalid requests or precision inputs throw `std::invalid_argument`. A valid
request returns a `FourSpinonSolve<Real>` with its sorted, canonical request and
one of `converged`, `unsupported_geometry`, `unconverged` or `precision_limit`.
Only a converged result contains a state. Numerical failures carry a diagnostic;
they never stand for a zero spectral weight. The form-factor call separately
returns the supported `heisenberg::FormFactor<Real>` status and optional weight.

`FourSpinonState<Real>` is an immutable shared snapshot. It retains the request,
real rapidities, centre, signed deviation, log absolute deviation, momentum
index and J=1 energy. Its rational residual is the unscaled maximum of the
original multiplicative equations; its quantum residual is the maximum of the
resolved logarithmic equations divided by N. Copying keeps the validated
snapshot alive. A moved-from state is not ready and its value accessors throw.

`roots()` returns standard ABA rapidities, with pair members at
centre +/- i(1/2+delta). This array alone is insufficient when delta is tiny:
`deviation()` and `log_abs_deviation()` remain part of the numerical state.
The public API does not accept edited root arrays as validated states.

## Finite deviations and stabilized matrices

Ideal Takahashi roots provide only the initial guess. Newton iteration then
solves for real roots, centre and log absolute deviation, retaining the sign.
The original pair self factor delta/(1+delta) is evaluated explicitly, avoiding
subtraction of rounded imaginary parts. Original Bethe equations and the exact
momentum label must pass before a state is returned.

The form factor uses algebraic row and column transformations of the finite-root
Gaudin and raising matrices. Both final matrices are assembled directly as real
nonsymmetric matrices. This removes inverse-deviation cancellation while
retaining the finite deviation; it is not an ideal-string approximation.
Logarithmic column scales avoid overflowing individual products. Uni20's normal
log-determinant dispatch chooses the factorization backend.

Small residuals and pivot diagnostics are rejection checks, not certified
relative error bounds. Increase arithmetic precision to test a difficult result.
General collapsed/singular-string regularization is outside this family.

## Numerical controls

- N=12 roots, energy and raising weight agree with the independent direct
  complex-equation/determinant fixture and its small-chain ED control described
  in [the complex-root guide](experimental-xxx-complex-roots.md).
- A negative N=16 deviation is checked against an independently evaluated
  arbitrary-precision weight, approximately 8.285230871529277e-8.
- An N=64 case retains delta approximately 1.41905199058e-94 when its displayed
  roots have imaginary parts exactly +/-1/2. Its stabilized raising weight is
  approximately 1.339859187815072e-11; the unreduced evaluator rejects those
  rounded roots. The compact N=16/N=64 references carry about 40 decimal digits.
- Tests also check the analytic Jacobian by finite differences, reflection,
  canonical labels, snapshot lifetime, explicit failures and MPFR at 128/256 bits.

The form-factor tests load 100-digit independent ground roots so they exercise
preparation and normalization without repeating the real-root ground iteration
for every scalar. Regenerate or verify those backgrounds using
`python experimental/reference/xxx_ground_fixture.py --digits 130 --verify`.
The example still solves the ground state from scratch.

The full complex-root weight controls can also be repeated independently:

```sh
python experimental/reference/xxx_finite_four.py --digits 150 --verify
python experimental/reference/xxx_finite_four.py --digits 180 --verify
```

This mpmath evaluator solves the original multiplicative equations and uses
unreduced complex determinants. It agrees with the N=16/N=64 deviation and
weight fixtures to 40 relative decimal places. It needs enough precision to
resolve the tiny pair in its complex roots; the stabilized C++ representation
avoids that particular precision requirement.

These selected controls validate this family and normalization, not completeness
of a spinon sum or the accuracy of a thermodynamic extrapolation.
