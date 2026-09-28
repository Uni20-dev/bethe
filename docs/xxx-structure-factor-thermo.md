# Thermodynamic XXX two-spinon structure factor

`bethe-xxx-structure-factor-thermo` evaluates the **spectral density** of the
infinite spin-1/2 antiferromagnetic XXX chain at zero field and temperature.
It implements the exact two-spinon contribution, not the full structure factor.
Unlike the [finite-ring tool](xxx-structure-factor.md), there are no Bethe roots
to solve, intermediate states to enumerate, or delta functions to broaden.

See the [heat-map tutorial](tutorials/xxx-structure-factor-thermo.md) for an
iMPS-oriented introduction.

## Run a frequency cut or a heat map

```sh
bethe-xxx-structure-factor-thermo --momentum 1.5 --points 501 --csv cut.csv
bethe-xxx-structure-factor-thermo --momentum 1.5 --omega 1.8 --precision long-double
bethe-xxx-structure-factor-thermo --momentum-points 129 --points 401 \
  --threads 4 --csv spectrum.csv --csv-table continuum=edges.csv
```

Commands assume an installed executable on `PATH`; from its build directory,
use `./bethe-xxx-structure-factor-thermo`. The common
[precision and export controls](command-line.md) apply, including optional fp128.
`--threads` uses the Uni20 scheduler for independent spectral points; output order
and numerical results do not depend on the worker count.

| Option | Meaning |
|---|---|
| `--momentum Q` | One momentum in radians/site, in $`[0,2\pi]`$ |
| `--momentum-points N` | Uniform momentum grid including both endpoints; default 65 |
| `--omega W` | One frequency, instead of a frequency grid |
| `--points N` | Uniform frequency grid including both endpoints; default 201 |
| `--omega-min W`, `--omega-max W` | Frequency range; default $`[0,\pi J]`$ |
| `--exchange J` | Positive antiferromagnetic coupling; default 1 |
| `--channel zz\|raising` | Longitudinal density or $`S^{-+}=2S^{zz}`$; default zz |
| `--tolerance T` | Absolute error target for the transition-rate logarithm; default 4096 machine epsilon |
| `--max-evaluations N` | Kernel integrand evaluations per point; default 200000 |

Single-point options exclude their corresponding grid options. Grids are limited
to one million total points. Inputs are validated before opening output files.

The frontend declares two tables:

- **`spectrum`** (primary): `q`, `omega`, `density`, `estimated_error`,
  `evaluations`, `status`.
- **`continuum`**: `q`, `lower`, `upper` at each momentum.

For example, `--csv-table continuum=edges.csv` exports the named auxiliary table.
An unknown table name is an error. See [output and exports](output.md).

## Conventions and exact expression

```math
H=J\sum_j \mathbf S_j\cdot\mathbf S_{j+1},\qquad J\gt0,
\qquad S_q^z=\frac1{\sqrt N}\sum_j e^{-iqj}S_j^z.
```

The returned density has the same convention as the finite tool:
$`S(q,\omega)=2\pi\sum_n w_n\delta(\omega-E_n+E_0)`$ before taking
the thermodynamic limit. It is **not** $`S/(2\pi)`$ and is not a discrete weight.

For $`J=1`$, write

```math
L(q)=\frac\pi2|\sin q|,\qquad U(q)=\pi\sin(q/2),\qquad
\cosh(\pi\rho)=\sqrt{\frac{U^2-L^2}{\omega^2-L^2}}.
```

Inside $`L\lt\omega\lt U`$,

```math
S^{zz}_2(q,\omega)=\frac{e^{-I(\rho)}}{2\sqrt{U^2-\omega^2}},\qquad
I(\rho)=\int_0^\infty\frac{e^x}{x}
\frac{\cosh(2x)\cos(4\rho x)-1}{\cosh x\sinh(2x)}\,dx.
```

These are Eqs. (9)–(12) of
[Caux and Hagemans](https://arxiv.org/abs/cond-mat/0611319).
Coupling scales energies and the density as
$`S_{2,J}(q,\omega)=S_{2,1}(q,\omega/J)/J`$.

## Thresholds are not numerical failures

The density vanishes outside the continuum and approaches zero at its upper
edge. At the lower edge it diverges, with a logarithmic enhancement of the
inverse-square-root singularity. At $`q=\pi`$ the low-frequency behavior is
$`\sqrt{\log(1/\omega)}/\omega`$.
[Karbach, Müller and Bougourzi](https://arxiv.org/abs/cond-mat/9606068)
analyze these limits.

| Status | Density/error fields | Exit status |
|---|---|---|
| `converged` | Finite value and estimated numerical error; upper-edge limit is zero | 0 |
| `outside_continuum` | Exactly zero | 0 |
| `lower_threshold` | Empty (JSON null): no finite point value exists | 0 |
| `numerical_failure` | Empty (JSON null): budget, precision or representable-range limitation | 2 |

At $`q=0,2\pi`$, the density is zero by total-spin conservation, including at
zero frequency. Negative frequencies are zero at zero temperature.
The continuum can become unresolvably narrow near zero momentum; this is not
silently treated as a computed peak.

The numerical error is an **estimate for the evaluated point**, not an interval
certificate. It does not account for uncertainty in the supplied frequency or
its distance from a threshold. Raising the precision is useful near an edge.
A plotted point grid is not an integration rule for these singular functions;
do not infer sum-rule saturation from a coarse rectangular sum. In particular,
the fixed-$`q=\pi`$ zeroth moment diverges even though the momentum-integrated
intensity is finite. Frequency-bin integrals and broadening are not yet provided
by this frontend.

## Stable evaluation and validation

The original integral has a conditionally convergent oscillatory tail.
With $`a=4|\rho|`$, the implementation subtracts
$`2(1-e^{-2x})\cos(ax)/x`$ and integrates that part analytically:
$`\log(1+4/a^2)`$. With $`u=e^{-2x}`$, the remainder is

```math
R(x,a)=\frac{2u^2(3-u^2)\cos(ax)-4u}{x(1-u)(1+u)^2}.
```

Near zero, an equivalent expression using `sin`, `tanh` and `expm1` avoids
cancellation. Beyond $`X\ge1`$, its absolute tail is bounded by
$`8e^{-2X}/X`$. Composite native-precision Gauss–Legendre rules use panels
no wider than half an oscillation. Rule differences, a roundoff allowance and
the tail bound supply the error estimate. Energies are compared through
factored differences and the transition rate is assembled logarithmically.

Tests cover independent 65-digit quadrature references, all enabled scalar
precisions, reflection symmetry, coupling scaling, threshold behavior, failure
statuses, and the published two-spinon fractions: approximately **72.89%** of
the full integrated intensity and **71.30%** of the first frequency moment.
These are physical partial-sum fractions, not tolerances or fitted scale factors.
The optional [reference script](../scripts/reference_xxx_thermodynamic_structure_factor.py)
uses mpmath and a different analytic tail subtraction.

The C++ entry point is
`bethe::heisenberg::ThermodynamicTwoSpinonStructureFactor<Real>` in
[`xxx_thermodynamic_structure_factor.hpp`](../include/bethe/xxx_thermodynamic_structure_factor.hpp).
It is immutable and can be shared between scheduler jobs. The rapidity-difference
kernel is separate from continuum kinematics for reuse by four-spinon weights.

## Next checkpoint: four spinons

The two-spinon checkpoint is implemented; four spinons are not yet exposed.
The next steps are the complex gamma-ratio matrix-element series, its truncation
diagnostics, and two-dimensional integration over the allowed pair energy and
momentum. Both momentum sectors modulo $`2\pi`$ must be included. Numerical
integration errors must remain distinct from missing six-and-higher-spinon weight.

Use Caux–Hagemans together with the
[author's correction to Eq. (31)](https://scipost.org/commentary/10.1088/1742-5468/2006/12/P12013/):
the upper frequency conditions for the excluded $`K_{2c}`$ and $`K_{2d}`$
intervals are interchanged in the original paper. The corrected conditions use
$`\pi\cos(k/2)`$ for $`K_{2c}`$ and $`\pi\sin(k/2)`$ for $`K_{2d}`$.
Verify the geometry independently by intersecting the two pair continua before
implementing a spectral grid. Then test permutation symmetry, coincident-rapidity
limits, precision/refinement stability and first-moment coverage before publishing
two-plus-four-spinon plots.
