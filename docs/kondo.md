# Kondo: universal zero-temperature impurity energy

**Status: equation-level target and independent reference oracle established;
native library and frontend not yet implemented.**

This is the next implementation target, not a claim that the repository
already solves finite-band or finite-size Kondo Hamiltonians.

## Physical scope and conventions

Take the isotropic, antiferromagnetic, single-channel spin-1/2 Kondo model
in its universal scaling limit, at zero temperature, with a **uniform**
magnetic field on impurity and host with equal g factors. This does not
cover ferromagnetic exchange, multiple channels, other impurity spins,
anisotropy, finite bandwidth, or a Kondo lattice.

The proposed input `b` is the full Zeeman splitting in energy units:
`H_field=-b*(S_imp^z+S_host^z)`. The positive scale `T_B` is defined below,
not an unspecified "Kondo temperature". Both are energies; set k_B=1.

Define `E_imp(b)=E_with_impurity(b)-E_clean_host(b)`, using the same host
and field in the subtraction. Return `Delta E_imp(b)=E_imp(b)-E_imp(0)`,
not the cutoff-dependent absolute impurity energy. The conjugate response
is `M_imp=-d Delta E_imp/db`: impurity-induced **total** magnetization,
including the host change, not an arbitrary finite-band local spin expectation.

[Barcza et al.](../CITATIONS.md#barcza-2020), Sec. VII, discuss why absolute
Bethe-ansatz ground energies cannot simply be compared with a tight-binding
host. Their Eqs. (4)-(5) use `B=b/2`. In Eqs. (124)-(126), set
`T_B=2*T1`, hence `x=|b|/T_B=B/T1`. Missing this factor of two would give
the wrong energy derivative and susceptibility normalization.

The exact zero-field scale calibration is

```text
a0 = 1/sqrt(2*pi*e)
chi0 = dM_imp/db at b=0 = a0/T_B
T_chi = 1/(4*chi0), so T_B = 4*a0*T_chi.
```

No automatic conversion from a bare exchange J and bandwidth D is planned:
that would require a specified regularization and scale-matching convention.

## Source equations

For x>=0, write `m(x)=M_imp(|b|)`. Eqs. (124)-(125) give

```text
a_n = (n+1/2)^(n-1/2) * exp(-n-1/2) / (2*sqrt(pi)*n!)
m(x) = sum_(n>=0) (-1)^n a_n x^(2n+1),                  0<=x<=1

A(w) = Gamma(w+1/2) * exp(w) * w^(-w),   A(0)=sqrt(pi)
C = 1/(2*pi^(3/2))
m(x) = 1/2 - C * integral_0^infinity sin(pi*w)/w * A(w) * x^(-2w) dw,
                                                              x>=1.
```

The magnetization is odd in b. The x=1 integral is only conditionally
convergent: `A(w)->sqrt(2*pi)`, so its tail is proportional to
`sin(pi*w)/w`. Ordinary finite-cutoff quadrature is not a sufficient
crossover algorithm. It is not acceptable to leave an unreported hole
around x=1 or publish a slowly truncated integral as converged.

## Energy derived from the response

Set `e(x)=Delta E_imp/T_B`. Integrating `e'(x)=-m(x)` gives

```text
e(x) = -sum_(n>=0) (-1)^n a_n x^(2n+2)/(2n+2),          0<=x<=1

e(x) = e(1) - (x-1)/2
       + C * integral_0^infinity sin(pi*w)/w * A(w)
             * (x^(1-2w)-1)/(1-2w) dw,                 x>=1.
```

The apparent singularity at w=1/2 is removable, with value `log(x)`.
Use `expm1((1-2*w)*log(x))/(1-2*w)` nearby. Energy is even in b, with
`e(0)=0`, `e(x)=-a0*x^2/2+O(x^4)`, and `e(x)/x -> -1/2` at high field.
These formulas and the energy normalization above are our derivation from
the response, not a finite-band ground-energy formula quoted from the paper.

## Native numerical implementation plan

Reuse `detail/quadrature.hpp`, compensated sums, Uni20 scalar precision,
and the common CLI/run-context/table/citation infrastructure. Keep this
model's crossover and subtraction logic separate from finite-root solvers.

1. **Low-field coefficients:** generate them by the stable recurrence
   `a_n/a_(n-1)=(n-1/2)/n * exp((n-1/2)*log1p(1/(n-1/2))-1)`.
   This avoids overflowing factorials/powers and does not require Gamma.
2. **Crossover:** investigate Euler acceleration of alternating series and
   of unit-interval sine lobes. Compare independent refinement levels and
   propagate quadrature/roundoff estimates through the transformation.
   A small last transformed term alone is not sufficient evidence.
3. **Scaled Gamma amplitude:** evaluate log(A) directly in native precision,
   not separate huge Gamma and tiny powers. A shifted Stirling expansion
   with a checked remainder is a candidate. For large w, the leading part
   is `w*log1p(1/(2*w))-1/2+log(2*pi)/2`, avoiding cancellation of
   `w*log(w)` terms. Audit Uni20's available scalar support before choosing
   a shared special-function facility; no hidden double fallback.
4. **Large ratios:** avoid requiring `b/T_B` to be representable when b and
   T_B individually are. Logarithmic ratios and the scaled energy `e(x)/x`
   give a path to high-field evaluation. Set explicit representability and
   tolerance contracts rather than silently overflowing.
5. **Failure semantics:** separate series, quadrature, oscillatory-tail and
   precision failures. All physical observables remain absent on failure.
   Analytic b=0 may bypass numerical work but not input validation.

The first public interface should accept signed field and positive T_B and
return energy change, induced magnetization and the zero-field susceptibility
calibration. Finite-temperature TBA and excited/finite-size states are later
capabilities, not implied by this zero-temperature implementation.

## Independent oracle and acceptance evidence

`scripts/reference_kondo.py` uses mpmath's arbitrary-precision alternating
summation and oscillatory quadrature. It is a development oracle, not a
production dependency or an implementation in the requested native scalar.

```sh
python3 scripts/reference_kondo.py --digits 40
python3 scripts/reference_kondo.py --digits 65
```

It separately evaluates both magnetization representations at x=1 and
rejects disagreement at the requested precision. Reference values at
x=0.1, 0.5, 1, 2 and 10 must stabilize with working precision before being
used in native regression tests. The energy reference uses the independently
summed low-field anchor e(1), not a fitted integration constant.

The 40- and 65-digit runs agree to at least 39 decimal places at these
points; values rounded here are development references, not production output:

| x | m(x) | e(x) |
| --- | --- | --- |
| 0.1 | 0.0241204367140108666055443531468 | -0.00120793395564716881820549509988 |
| 0.5 | 0.112566048096924547924356687485 | -0.0291472098839518339995212394642 |
| 1 | 0.192111646210937727816038464640 | -0.106792242303035303923031867865 |
| 2 | 0.275163637910673046871367248582 | -0.346291700274578677809068078930 |
| 10 | 0.387179325881091209470177660865 | -3.17987143725043217796651309659 |

Acceptance tests for the native implementation must cover:

- both representations at crossover and neighboring fields;
- magnetization oddness, energy evenness, and scale covariance;
- `-d Delta E_imp/db=M_imp` with native-precision finite differences;
- the susceptibility calibration and quadratic low-field energy;
- logarithmic approach to saturation and the high-field energy slope;
- fp64, long-double and fp128 references that detect narrowing;
- separate budgets, near-zero fields, extreme ratios, invalid inputs, and
  missing failed observables in every export format.

The original [Andrei solution](../CITATIONS.md#andrei-1980) remains the
catalogue's foundational reference. This first implementation deliberately
uses an observable whose subtraction and scale can be stated without
pretending that a continuum BA cutoff is a generic lattice band.
