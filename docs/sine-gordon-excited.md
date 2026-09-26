# Sine-Gordon: selected exact finite-volume excitations

[Overview](../README.md) · [Vacuum and conventions](sine-gordon.md) ·
[Particle lines and Bethe–Yang](sine-gordon-excitations.md) · [Exports](output.md)

`bethe-sine-gordon-excited` solves an excited-state nonlinear integral equation
(NLIE), including the finite-volume sea contribution. Unlike the Bethe–Yang
tool, it does not discard wrapping corrections. “Exact” describes the continuum
equations; their numerical solution still has a tolerance and can fail.

## Which states?

The implemented family has **p>=1**, two real holes at opposite rapidities,
half-odd Bethe numbers **I=1/2 or 3/2**, total momentum zero and topological
winding charge **+2 or -2**. These are two solitons or two antisolitons, not a
neutral soliton–antisoliton pair. Charge conjugation makes the two signs
degenerate, so `--charge` changes the sector label, not the equation.

This is the delta=0, no-complex-root, no-special-root family in
[Feverati–Ravanini–Takács](../CITATIONS.md#feverati-ravanini-takacs-1999),
Eqs. (3.13), (3.15), (3.18), and Sec. 5.2.2. The half-odd quantization is
part of the physical state selection: switching to integer numbers is not
another level of this same local sine-Gordon sector. I=1/2 approaches its
lowest conformal primary; I=3/2 selects a particular left/right descendant.
There is no completeness or degeneracy-counting claim.

Neutral pairs, moving pairs, other Bethe labels, attractive excited states,
complex roots and continuation through special-root rearrangements are not
implemented. In particular, this finite-volume solver **does not cover p=1/3**:
the Schwinger comparison in the particle guide remains a separate light-sector
scaling-limit benchmark, not a finite-volume gauge-theory solution.

## Running and comparing

```sh
bethe-sine-gordon-excited --p 2 --length 1
bethe-sine-gordon-excited --p 2 --length 0.1 --number 1.5 --json descendant.json
bethe-sine-gordon-excited --p 2 --length 1 --precision long-double --csv-table gap=gap.csv
bethe-sine-gordon-excited --p 1 --length 1 --precision fp128 --json free.json
```

The physical inputs are soliton mass M (default one), circumference L, coupling
p and velocity=hbar=1, with the [vacuum solver's canonical normalization](sine-gordon.md#physical-convention).
The default precision is fp64; long-double and optional MPLAPACK fp128 remain
native throughout the source, quadrature, nonlinear solve and hole quantization.
Interacting fp128 calculations can take many minutes; the last example is free.

Three tables keep the energy reference and the numerical diagnostics explicit:

- `levels`: separately computed vacuum and two-soliton **bulk-subtracted**
  energies, their scaling functions, errors and work counters.
- `source`: positive Bethe number, charge, positive rapidity, quantization
  residual, hole-position error estimate, source quadrature estimate and Fourier
  tail bound. The other hole is its negative.
- `gap`: vacuum-relative energy and `scaled_gap=L*gap/(2*pi)`. This is the table
  to compare with a finite-volume excitation calculation.

```math
E_C=E-L e_{\mathrm{bulk}},\qquad Y=L E_C,\qquad
\Delta E=E_C^{\mathrm{exc}}-E_C^{\mathrm{vac}}.
```

Neither level is an absolute total energy or energy per lattice site. A common
bulk subtraction cancels in the gap. We do not evaluate the convention-dependent
bulk term, including at its resonant couplings. To compare a lattice MPS model,
first match its velocity, mass scale, continuum coupling and charge/boundary
sector; finite-spacing effects are not numerical errors of this solver.

Exit codes are 0 for success, 1 for invalid input/output errors, and 2 for an
incomplete solve. A failed state's energy and rapidity are missing, not an
unconverged last iterate. A gap is published only if **both** states converge.
Unavailable estimates are null in JSON and empty in CSV/TSV. Exports use the
shared metadata, CPU-time, streaming and overwrite-protection conventions;
`--references` prints the literature separately from ordinary output.

## Equation and branch audit

Let G be the existing native Fourier kernel and define the odd scattering phase
by its zero at the origin:

```math
\chi(z)=2\pi\int_0^zG(w)\,dw,\qquad
g(z)=\chi(z-H)+\chi(z+H),\qquad u=ML.
```

For a positive contour shift eta, the real-axis counting equation is

```math
Z(\theta)=u\sinh\theta+g(\theta)
+2\,\mathrm{Im}\int_{-\infty}^{\infty}
G(\theta-x-i\eta)\log(1+e^{iZ(x+i\eta)})\,dx,
\qquad Z(H)=2\pi I.
```

The shared vacuum iteration uses `epsilon(x)=-i*Z(x+i*eta)` and
`A(x)=log(1+exp(-epsilon(x)))`. The excited solver adds `-i*g(x+i*eta)` to
its driving term and couples that solution to a safeguarded real hole solve.
The same hole position is used in the contour source and in the real-axis
counting kernel. Each hole trial resolves the nonlinear sea; it is not a
Bethe–Yang root followed by an energy correction.

```math
E_C^{\mathrm{exc}}=2M\cosh H-\frac{M}{\pi}\,\mathrm{Im}
\int_{-\infty}^{\infty}\sinh(x+i\eta)\,A(x)\,dx.
```

The sea energy uses the **bare** driving term, not the added hole source.
Parity permits the same half-grid convolution as the vacuum. The regular
principal-log implementation requires `Re(1+exp(-epsilon))>0` on its contour
at every iteration. A violation is reported as a numerical/branch failure;
the solver does not silently unwrap logarithms or insert extra roots.
Two independently resolved contours must agree in energy and hole coordinate.
This check is numerical evidence for the selected regular branch, not a general
classification of sine-Gordon states.

## Error controls and validation

`--tolerance` is absolute in Y per state, default **16777216 native epsilons**.
This is an estimated numerical accuracy, not an interval certificate or a
relative bound for exponentially small finite-volume corrections. A gap can
accumulate errors from both level solves; no stronger gap bound is claimed.

The [vacuum numerical options](sine-gordon.md#running-the-vacuum-solver) also
apply here. `--max-iterations` (32768) is a total nonlinear-update budget **per state**,
including every hole trial and both contours. `--max-root-iterations` (1024)
is the total number of hole trials across meshes/cutoffs/contours.
Both exhausted iteration budgets use `iteration_limit`; their separate counters
identify the exhausted resource. Fourier work is budgeted per kernel/source
table, not globally; total evaluation counts are reported.

Source quadrature and its analytic Fourier tail bound are separate from hole
quantization. `hole_error` estimates the effect of root displacement and inner
noise on Y and `2*u*sinh(H)`. The inner NLIE is solved more tightly than the
outer request. Two successive rapidity mesh agreements, adjacent cutoff
agreements and an independent contour comparison are required; these compare
both Y and the hole coordinate. `kernel_error` estimates the sea-kernel table
error, while `source_quadrature_error` and `source_tail_bound` describe the hole
tables. Small residuals alone do not certify a resolved level.

Within a fixed grid, the previous converged sea can seed the next hole trial.
A converged coarse-grid hole position likewise seeds the next finer grid, but
the sea is rebuilt and the full new counting equation must converge. Seeds are
not carried between the two independently solved contours. These initial guesses
reduce redundant work; they do not replace any refinement or residual check.

At the free Dirac point the kernel vanishes, and the sea is exactly the vacuum
sea. This gives an especially direct normalization and gap check:

```math
\Delta E=2\sqrt{M^2+(2\pi I/L)^2},\qquad p=1.
```

In the ultraviolet, FRT Sec. 5.2.2 gives

```math
R^2=\frac{p+1}{2p},\qquad
\Delta_+=\Delta_-=\frac{R^2}{2}+I-\frac12,\qquad
\frac{L\Delta E}{2\pi}\longrightarrow R^2+2I-1.
```

Regression checks approach these limits for both labels at p=1, 1.5, 2 and 3,
and recover Bethe–Yang at large volume while resolving its omitted sea term.
The [independent reference script](../scripts/reference_sine_gordon_excited.py)
uses the closed p=2 rapidity kernel, a nonuniform Gauss grid and a separate root
solver. Refining 256 to 512 nodes with two contours/cutoffs gives, at M=L=1,
bulk-subtracted energies 5.11055662981005 and 17.14609634255262 for the two
labels. These are fp64 reference values, not high-precision fixtures.
Native closed-kernel source tests and free-point tests check for precision
narrowing; budget and input tests check that failed solves publish no level.

## Library entry point

```cpp
#include <bethe/sine_gordon_excited.hpp>
namespace sg = bethe::sine_gordon;
sg::TwoSolitonOptions<long double> options;
auto state = sg::two_soliton_level(1.0L, 1.0L, 2.0L, uni20::from_twice(1), options);
if (state.converged) {
  // *state.casimir_energy is E-L*e_bulk, NOT the vacuum-relative gap.
  // *state.rapidity is H; the other hole is -H.
}
```

`TwoSolitonState` shares the vacuum's diagnostic fields, but its inherited
`effective_central_charge` is intentionally absent. A charged excited level is
not a vacuum central-charge measurement. Subtract a separately converged
`vacuum_energy` with the same mass, circumference, coupling and conventions to
obtain the physical excitation gap, as the frontend does.
