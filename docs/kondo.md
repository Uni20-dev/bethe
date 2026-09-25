# Kondo: universal zero-temperature impurity energy

**Status: native-precision zero-temperature response library and frontend implemented.**

This is a scaling-limit impurity response, not a finite-band or finite-size
Kondo Hamiltonian solver.

## Command-line use

```sh
bethe-kondo-response --field 2 --scale 1
bethe-kondo-response --field -0.1 --scale 1 --precision long-double
bethe-kondo-response --field 2 --scale 1 --precision fp128 --json response.json --csv response.csv
bethe-kondo-response --references
```

Both `--field` and `--scale` are required. Precision defaults to fp64;
fp128 requires an enabled MPLAPACK build. The `response` table contains the
impurity energy change, induced magnetization, zero-field susceptibility,
error estimates, work counters and convergence status. The preamble and
export metadata record the conventions below, numerical controls and CPU time.
JSON, CSV and TSV files use the shared output options; references appear only
with `--references`.

`--tolerance`, `--max-series-terms`, `--max-lobes`, `--max-evaluations` and
`--max-quadrature-levels` expose the library controls below. An incomplete
calculation exits with status 2 and publishes missing observables (JSON null,
empty CSV/TSV cells), not a partial physical answer. Invalid input exits with
status 1 before opening output files, including existing `--force` targets.

## Physical scope and conventions

Take the isotropic, antiferromagnetic, single-channel spin-1/2 Kondo model
in its universal scaling limit, at zero temperature, with a **uniform**
magnetic field on impurity and host with equal g factors. This does not
cover ferromagnetic exchange, multiple channels, other impurity spins,
anisotropy, finite bandwidth, or a Kondo lattice.

The input $`b`$ is the full Zeeman splitting in energy units:
$`H_{\mathrm{field}} =-b \,(S_{\mathrm{imp}} ^z +S_{\mathrm{host}} ^z)`$. The positive scale $`T_{B}`$ is defined below,
not an unspecified "Kondo temperature". Both are energies; set k_B=1.

Define `E_imp(b)=E_with_impurity(b)-E_clean_host(b)`, using the same host
and field in the subtraction. Return `Delta E_imp(b)=E_imp(b)-E_imp(0)`,
not the cutoff-dependent absolute impurity energy. The conjugate response
is `M_imp=-d Delta E_imp/db`: impurity-induced **total** magnetization,
including the host change, not an arbitrary finite-band local spin expectation.

[Barcza et al.](../CITATIONS.md#barcza-2020), Sec. VII, discuss why absolute
Bethe-ansatz ground energies cannot simply be compared with a tight-binding
host. Their Eqs. (4)-(5) use $`B =b /2`$. In Eqs. (124)-(126), set
`T_B=2*T1`, hence `x=|b|/T_B=B/T1`. Missing this factor of two would give
the wrong energy derivative and susceptibility normalization.

The exact zero-field scale calibration is

```math
a_0=\frac1{\sqrt{2\pi e}},\qquad
\chi_0=\left.\frac{dM_{\mathrm{imp}}}{db}\right|_{b=0}=\frac{a_0}{T_B},
\qquad T_\chi=\frac1{4\chi_0},\quad T_B=4a_0T_\chi.
```

No automatic conversion from a bare exchange J and bandwidth D is planned:
that would require a specified regularization and scale-matching convention.

## Source equations

For x>=0, write $`m (x)=M_{\mathrm{imp}} (\lvert b \rvert)`$. Eqs. (124)-(125) give

```math
\begin{aligned}
a_n&=\frac{(n+1/2)^{n-1/2}e^{-n-1/2}}{2\sqrt\pi\,n!},\\{}
m(x)&=\sum_{n\ge0}(-1)^na_nx^{2n+1},\qquad 0\le x\le1,\\{}
A(w)&=\Gamma(w+1/2)e^ww^{-w},
\qquad A(0)=\sqrt\pi,\qquad C=\frac1{2\pi^{3/2}},\\{}
m(x)&=\frac12-C\int_0^\infty\frac{\sin(\pi w)}{w}
A(w)x^{-2w}\,dw,\qquad x\ge1.
\end{aligned}
```

The magnetization is odd in b. The x=1 integral is only conditionally
convergent: $`A (w)\to \sqrt{2\,\pi }`$, so its tail is proportional to
$`\sin (\pi \,w)/w`$. Ordinary finite-cutoff quadrature is not a sufficient
crossover algorithm. It is not acceptable to leave an unreported hole
around x=1 or publish a slowly truncated integral as converged.

## Energy derived from the response

Set $`e (x)=\Delta E_{\mathrm{imp}} /T_{B}`$. Integrating `e'(x)=-m(x)` gives

```math
\begin{aligned}
e(x)&=-\sum_{n\ge0}\frac{(-1)^na_nx^{2n+2}}{2n+2},\qquad 0\le x\le1,\\{}
e(x)&=e(1)-\frac{x-1}{2}
+C\int_0^\infty\frac{\sin(\pi w)}{w}A(w)
\frac{x^{1-2w}-1}{1-2w}\,dw,\qquad x\ge1.
\end{aligned}
```

The apparent singularity at w=1/2 is removable, with value `log(x)`.
Use `expm1((1-2*w)*log(x))/(1-2*w)` nearby. Energy is even in b, with
`e(0)=0`, `e(x)=-a0*x^2/2+O(x^4)`, and `e(x)/x -> -1/2` at high field.
These formulas and the energy normalization above are our derivation from
the response, not a finite-band ground-energy formula quoted from the paper.

## Library use and numerical contract

```cpp
#include <bethe/kondo.hpp>
auto state = bethe::kondo::ground_response<long double>(2.0L, 1.0L); // b, T_B
if (state.converged) {
  auto energy_change = *state.energy_change;
  auto induced_magnetization = *state.magnetization;
  auto chi0 = *state.zero_field_susceptibility;
}
```

`Options<Real>` defaults to tolerance `1048576*epsilon<Real>`, at most 256
low-series terms, 256 unit-interval Fourier lobes, 200000 total quadrature
evaluations and 12 quadrature refinement levels per lobe. Series and lobe
counts have a hard work ceiling of 4096. Real parameters must be finite,
with strictly positive T_B and tolerance. Negative b is supported.

The tolerance is absolute in **magnetization and Delta E/|b|**, not in
physical energy or Delta E/T_B. Returned `magnetization_error` and
`scaled_energy_error` use those same conventions. At b=0, both observables
and error estimates are zero and no numerical work is required; chi0 still
must be representable. Input validation is never bypassed.

Errors combine successive refinement differences, propagated quadrature
estimates and conservative floating-point estimates. They are **not rigorous
interval enclosures**. Two consecutive agreements on 16-term/lobe refinement
blocks are required; a single small last term is not the stopping criterion.
For small x, the low-field calculation resolves M/x and e(x)/x^2 before
rescaling, avoiding a spurious loss of relative low-field accuracy.

All three physical observables are optional and remain absent on failure:
`series_limit`, `quadrature_limit`, `tail_limit` or `precision_limit`.
The work counters include attempted calculations. Error estimates that were
not completed remain infinite. A nonzero field producing an underflowed
magnetization or energy, or an unrepresentable chi0, is reported as
`precision_limit`, not a successful zero response.

## Native numerical implementation

The library reuses `detail/quadrature.hpp`, compensated sums and Uni20 scalar
precision. Its crossover and subtraction logic is separate from finite-root
solvers. No mpmath, Gamma-library or double-precision fallback is used at run time.

1. **Low-field coefficients:** generated by the stable recurrence
   `a_n/a_(n-1)=(n-1/2)/n * exp((n-1/2)*log1p(1/(n-1/2))-1)`.
   This avoids overflowing factorials/powers and does not require Gamma.
2. **Crossover:** Euler acceleration is used for both the alternating
   low-field series and the unit-interval sine lobes. Half-differences keep
   intermediate values bounded while propagating input-error estimates.
   Numerical integration uses the shared native tanh-sinh quadrature.
3. **Scaled Gamma amplitude:** log(A) is evaluated directly in native precision,
   not as separate huge Gamma and tiny powers. Recurrence shifts the Gamma
   argument to at least 32, followed by 16 Stirling corrections. The next
   omitted term bounds the positive-real asymptotic remainder;
   floating-point error is estimated separately. See
   [DLMF 5.11.1 and 5.11(ii)](https://dlmf.nist.gov/5.11).
   For large w, the leading part
   is `w*log1p(1/(2*w))-1/2+log(2*pi)/2`, avoiding cancellation of
   $`w \,\log (w)`$ terms. This is a model-specific scaled amplitude, not a new
   general-purpose Gamma function.
4. **Large ratios:** high-field evaluation uses logarithmic ratios and
   $`e (x)/x`$; b/T_B itself need not fit the scalar. The energy integrand uses
   $`(x ^{-2w }-x ^{-1})/(1-2w)`$, evaluated with `expm1` near its removable
   singularity. Physical energy is recovered as $`\lvert b \rvert \,e (x)/x`$.

The frontend uses the common CLI, run context, tables and citations.
Finite-temperature TBA and excited/finite-size states
are later capabilities, not implied by this zero-temperature implementation.

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

Additional 45-digit oracle values at x=1.0001, 1.01 and 1000000 test the
high-field representation close to crossover and far into saturation.

Native regression tests cover:

- both representations at crossover and neighboring fields;
- magnetization oddness, energy evenness, and scale covariance;
- `-d Delta E_imp/db=M_imp` with native-precision finite differences;
- the susceptibility calibration and quadratic low-field energy;
- logarithmic approach to saturation and the high-field energy slope;
- fp64, long-double and fp128 references that detect narrowing;
- separate budgets, near-zero fields, overflowing ratios, invalid inputs,
  and missing failed observables;
- frontend precision dispatch, metadata, JSON/CSV/TSV export, budget forwarding,
  streaming output and invalid-input protection of existing files.

The original [Andrei solution](../CITATIONS.md#andrei-1980) remains the
catalogue's foundational reference. This first implementation deliberately
uses an observable whose subtraction and scale can be stated without
pretending that a continuum BA cutoff is a generic lattice band.
