# Sine-Gordon: finite-volume vacuum roadmap and kernel foundation

**Status: native-precision complex scattering kernel implemented.** The vacuum
integral-equation solver and frontend are not yet implemented. This is a field
theory calculation, not a classical sine-Gordon PDE evolution or a finite-site
spin-chain diagonalization.

## Physical convention for the planned vacuum solver

Use the canonically normalized kinetic term `(1/2)(partial phi)^2`, coupling
`0 < beta^2 < 8*pi`, and the parameter

```text
p = beta^2 / (8*pi-beta^2).
```

Then p>1 is repulsive, 0<p<1 attractive, and p=1 the free massive Dirac point.
In Rutkevich's normalization, his beta squared is our `beta^2/(8*pi)`,
his xi is p, and his gamma is `pi/(1+p)`. State these conventions explicitly;
the different meanings of beta must not be mixed.

Take the soliton mass M as the physical input scale, set velocity and hbar to
one, and use `u=M*L`. The intended observable is the **bulk-subtracted** vacuum
energy `E_C=E_0-L*e_bulk`, or the dimensionless scaling function `Y=L*E_C`.
No nonuniversal absolute bulk energy is inferred. In the attractive regime M
still denotes the soliton mass, not necessarily the lightest breather mass.

The starting point is [Destri–de Vega](../CITATIONS.md#destri-de-vega-1992).
For explicit shifted-contour formulas and checks on normalization we use
[Rutkevich (2020)](../CITATIONS.md#rutkevich-2020), equations (22), (68), (69).

## Implemented kernel

```cpp
#include <bethe/sine_gordon.hpp>
auto g = bethe::sine_gordon::scattering_kernel(std::complex<double>(0.8, 0.4), 2.7);
if (g.converged) {
  // *g.value is the complex kernel, not an energy.
}
```

The kernel is

```text
G_p(z) = integral_0^infinity dk cos(k*z)/(2*pi)
         * sinh((p-1)*pi*k/2) / [sinh(p*pi*k/2)*cosh(pi*k/2)].
```

Its Fourier representation is used only inside
`|Im(z)| < pi*min(1,p)`. Invalid couplings, nonfinite arguments, points outside
this strip and nonpositive tolerances throw. We enforce the same strip at p=1
even though the identically zero kernel can itself be continued further.

`KernelOptions<Real>` provides an absolute component tolerance (1024 native
epsilons), up to 100000 integrand evaluations, 16 quadrature levels, and 64
Fourier-cutoff attempts. All arithmetic stays in the chosen fp64, long-double
or fp128 type. The implementation reuses the common compensated tanh-sinh
quadrature; it does not introduce a double-only special-functions dependency.

The hyperbolic ratio is evaluated with `expm1` and decaying exponentials.
For complex arguments, growth from `cos(k*z)` is combined with the Fourier
decay *before* either factor is evaluated. This avoids an overflowing cosh
multiplied by an underflowing ratio on otherwise legal contours. The same
formula retains the small multiplier for p only a few native epsilons from one.

Let `d=pi*min(1,p)-|Im(z)|>0`. Beyond a cutoff K, an absolute tail bound is

```text
exp(-d*K) / [pi*d*(1-exp(-p*pi*K))].
```

The returned `tail_bound` and `quadrature_error` are separate. The latter is a
successive-mesh estimate, not a rigorous enclosure. A kernel value is optional
and is published only when their sum meets the requested tolerance. Failure
statuses distinguish cutoff, quadrature and finite-range limitations. Unknown
error diagnostics remain infinite rather than appearing to be zero.

Independent checks include the exact transforms

```text
G_(1/2)(z) = -1/(2*pi*cosh(z)),
G_2(z) = z/(2*pi^2*sinh(z)),  G_2(0)=1/(2*pi^2),
G_1(z) = 0,                  dG_p(0)/dp at p=1 = 1/8.
```

The first two are test references, not implementation shortcuts. Tests also
cover evenness, complex conjugation, a separate 90-digit quadrature reference,
overflow-prone Fourier factors, native near-free couplings, and failure budgets.

## Next implementation: shifted-contour vacuum NLIE

For `0 < eta < pi*min(1,p)/2`, use

```text
A(theta) = log(1+exp(-epsilon(theta)))
epsilon(theta) = -i*u*sinh(theta+i*eta)
                - integral G_p(theta-t)*A(t) dt
                + integral G_p(theta-t+2*i*eta)*conj(A(t)) dt
Y = -u/pi * Im integral sinh(theta+i*eta)*A(theta) dtheta.
```

The finite contour shift regularizes the counting-function logarithms without
crossing kernel poles. Do not move directly to eta=pi/2 in the attractive
regime. The chosen eta is numerical, not a physical model parameter: agreement
at two independently resolved shifts is an important check.

Implementation requirements:

1. Tabulate the real and shifted kernels on a difference grid so repeated
   convolutions reuse them. Keep their error budgets below the desired final
   energy accuracy. Reuse the stable complex `log(1+z)` utility.
2. Solve for the interaction correction to the known driving term, with a
   reported nonlinear residual and an explicit iteration budget.
3. Refine rapidity spacing and increase the rapidity cutoff independently.
   A small nonlinear residual alone must not publish a converged vacuum energy.
4. Check contour independence, the free Dirac integral, the ultraviolet
   `Y -> -pi/6` limit, and infrared wrapping contributions. Attractive-coupling
   tests must include the breather contribution, not just a soliton gas.
5. Publish M, L, p, the subtraction convention, and independent numerical
   diagnostics through the common CLI/run-context/table infrastructure.

This first target is the untwisted zero-topological-charge vacuum. Excited
states require additional source terms and branch/quantization bookkeeping;
large-volume Bethe–Yang equations alone are not an exact finite-volume solver.
