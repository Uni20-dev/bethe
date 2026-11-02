# Experimental finite XXX complex-root helpers

[Experimental components](experimental.md) · [Finite XXX structure factor](xxx-structure-factor.md)

`<bethe/experimental/xxx_complex_roots.hpp>` supplies regular, unreduced
finite-root equations and form factors in
`bethe::experimental::xxx::complex_roots`. They are mathematical building
blocks for family-specific solvers, not a complete finite-chain state scanner.
All roots are standard ABA rapidities lambda; the public real-root solver's
coordinate is z=2 lambda.

## Equations and precision

The helpers evaluate the multiplicative Bethe equations through their logarithms:

```math
L_i=N\left[\log(\lambda_i+\mathrm{i}/2)-\log(\lambda_i-\mathrm{i}/2)\right]
+\sum_{j\ne i}\left[\log(\lambda_i-\lambda_j-\mathrm{i})
-\log(\lambda_i-\lambda_j+\mathrm{i})\right].
```

`rational_residual` measures the largest absolute value of exp(L_i)-1,
using an expm1/trigonometric expression near zero. `gaudin` returns the
unscaled Gaudin matrix; differentiation gives dL_i/dlambda_j = -i G_ij.
`energy` evaluates J=1 energy and rejects a nonreal or nonfinite result.

The residual allowance accounts for proximity to the equation poles at
lambda = +/- i/2 and root differences = +/- i. A root coincidence is still
rejected by the separate geometry gate, but does not enlarge this allowance:
the scattering factors are regular there. This distinction prevents an
iterate approaching duplicate real roots from passing a progressively looser
equation check. The allowance remains a numerical screening rule, not a
rigorous bound on root or form-factor errors.

| Helper | Inputs after the real template argument |
| --- | --- |
| `log_equation` | sites, complex-root span, root index, real exemplar |
| `rational_residual`, `gaudin`, `residual_tolerance`, `energy` | sites, complex-root span, real exemplar |
| `separation` | complex-root span, real exemplar |
| `phase` | real angle |
| `bethe_atan` | complex value |

Span helpers require an explicit finite-precision exemplar for MPFR, including
empty spans. Native callers may omit it. `phase` and `bethe_atan` derive precision
from their scalar argument. Exact, unset or mixed-precision roots must be
converted explicitly before a numerical call. Arithmetic stays in the selected
type: fp64, native long double, optional binary128 or MPFR/MPC.

`bethe_atan` fixes the value on the imaginary-axis cuts to the odd
quantum-number convention, independently of signed real zero. Away from those
cuts it uses Uni20's principal complex inverse tangent. This convention belongs
to the Bethe equations; it does not change the underlying scalar math library.

## Unreduced form factors

`raising_weight(prepared_singlet, roots)` reuses a validated
`PreparedSinglet<Real>` and constructs the excited Gaudin and raising matrices
through the current Uni20 determinant dispatch. It returns the same
`heisenberg::FormFactor<Real>` status/optional-weight record as the supported
real-root API. The weight is N times the normalized local raising matrix element
squared, with no 2*pi factor or thermodynamic spinon coefficient.

Family callers can supply a third argument, the validated **absolute excited
momentum index** in `[0,N)`. It is checked against the roots using the existing
precision-scaled residual tolerance, which must also distinguish neighbouring
lattice momenta. When that index equals the singlet ground momentum, total
spin raising annihilates the singlet: the result is `converged` with exactly
zero `weight` and absent `log_weight`. The ground momentum index is zero for
N divisible by four and N/2 otherwise. Root, precision, original-equation and
Gaudin-norm checks still precede this selection zero; a bad state cannot acquire
a successful zero by supplying a momentum label. A mismatched label gives
`roots_unconverged`; an index outside `[0,N)` throws `std::invalid_argument`.

The two-argument call evaluates the unreduced determinant without inferring a
discrete label from a rounded phase. At a symmetry-forbidden state it may give
a precision failure or a small numerical residual weight. Family solvers
should supply the momentum obtained from their validated quantum numbers.

The two products in a complex-root raising column are evaluated independently.
They are not conjugates, unlike the real-root specialization. Ground-only
normalization is reused from the prepared context. The rational determinant
conventions follow the Caux–Hagemans–Maillet formulas described in the
[form-factor guide](xxx-structure-factor.md).

Malformed counts or scalar precision throw `std::invalid_argument`. Numerical
failure returns an absent weight: off-shell roots give `roots_unconverged`,
while unresolved poles/determinants give `precision_limit`. Family labels,
conjugate geometry, momentum and the intended branch must still be checked by
the family-specific state solver. A small residual alone does not identify a
state family or certify its weight's relative error.

The unreduced formula deliberately rejects nearly singular strings. Its
separation gate and pivot diagnostic are rejection indicators, not error
bounds. A rounded ideal-string pair cannot substitute for a separately retained
finite deviation. Cancellation-aware one-pair and higher-family representations
are separate work; no ideal-string fallback is used here.

## Independent controls

An analytic N=4 control gives excited energy -1 and Fourier raising weight 4/3.
An N=6 triplet checks that zero absolute momentum is not confused with zero
momentum transfer from its momentum-pi ground state. A nonsingular N=12
zero-momentum triplet checks exact zero output, a resolved norm, and rejection
of off-shell, singular, incorrectly labelled and mixed-precision inputs.
The tests also differentiate the logarithmic equations numerically to check
Gaudin entries and check both signed-zero sides of the atan convention. A resolved
N=12 state with holes {-3,-1,0,2}, two-string branch 1 and momentum index 9
checks genuinely complex products against an independent mpmath evaluator:

- Energy: approximately -1.43478610347334284.
- Fourier raising weight: approximately 0.0000120144194653219899.
- Spin-basis exact diagonalization agrees in energy within 5e-16 and in weight
  within 7e-15 relative error; this small native check is not a claim of that
  accuracy for arbitrary inputs.

The committed fixture carries 100 decimal digits. Regenerate the mathematical
control from coarse seeds, without any private data or compiled Bethe extension:

```sh
python3 -m venv /path/to/reference-env
/path/to/reference-env/bin/pip install 'mpmath>=1.2'
/path/to/reference-env/bin/python experimental/reference/xxx_regular_roots.py --verify
/path/to/reference-env/bin/python experimental/reference/xxx_regular_roots.py --digits 130 --verify
/path/to/reference-env/bin/python experimental/reference/xxx_zero_momentum.py --verify
```

The independent program solves all five complex multiplicative equations and
uses direct determinant/product formulas. It does not share the C++ logarithmic
form-factor code or a family solver's real/pair Jacobian. The routine C++ tests
consume the compact fixture and do not run Python or diagonalize a large matrix.
