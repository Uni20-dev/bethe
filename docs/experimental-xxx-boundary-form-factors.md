# Free-end XXX boundary raising weights

[Experimental components](experimental.md)

`bethe/experimental/xxx_boundary_form_factors.hpp` evaluates normalized local
raising weights between an even, free-end XXX singlet and a supported
positive-real-root highest-weight triplet. It reuses the open-chain solver
and caches an immutable ground-state normalization for repeated evaluations.

## Conventions and use

```math
H=\sum_{j=1}^{N-1}\mathbf S_j\cdot\mathbf S_{j+1},\qquad
J=1,\qquad S^\alpha=\sigma^\alpha/2,\qquad u=z/2.
```

For normalized states the result is

```math
w_e^+=|\langle e;1,1|S_1^+|g\rangle|^2,\qquad w_e^{zz}=w_e^+/2.
```

There is no Fourier factor of N or 2 pi. Reflection gives equal left/right
endpoint squared weights. The complete singlet raising norm is 1/2, and the
longitudinal norm is 1/4; a restricted real-root family retains its actual
weight. It is never rescaled to exhaust either norm. Squared weights alone
do not determine signs of end-to-end or same-sector matrix elements.

```cpp
#include <bethe/experimental/xxx_boundary_form_factors.hpp>

namespace obc = bethe::heisenberg::open;
namespace boundary = bethe::experimental::xxx::open;
auto ground = obc::ground_state<double>(8);
auto excited = obc::sector_ground_state<double>(8, uni20::half_int{1});
boundary::BoundaryFormFactors<double> prepared(8, ground);
auto factor = prepared(excited);
if (!factor.converged())
    throw std::runtime_error("boundary form factor rejected");
auto raising_weight = *factor.weight;
auto longitudinal_weight = raising_weight / 2;
```

With `BETHE_BUILD_EXPERIMENTAL=ON` and `BETHE_BUILD_EXAMPLES=ON`, build and
run the fresh-calculation example from your build directory:

```sh
cmake --build . --target bethe_example_boundary_form_factors
./experimental/bethe_example_boundary_form_factors
```

It prints the N=8 ground energy, lowest-triplet gap and raising/longitudinal
weights. It needs neither archived data nor a cluster account. Experimental
examples and headers are build-tree components, not installed programs.

For MPFR, give the solver an explicit precision; every supplied root must
have that same finite working precision:

```cpp
using R = uni20::mpreal;
bethe::SolverOptions<R> options(uni20::Precision::bits(160));
auto ground_mp = obc::ground_state<R>(8, options);
auto excited_mp = obc::sector_ground_state<R>(8, uni20::half_int{1}, options);
boundary::BoundaryFormFactors<R> prepared_mp(8, ground_mp);
auto factor_mp = prepared_mp(excited_mp);
```

The prepared object owns its snapshot. Copies share immutable storage; later
input mutation or destruction is harmless. Moves leave the source unavailable,
with absent residual, energy and working precision. Concurrent const evaluations
use their own scratch. `ground_energy()` is recomputed from the validated roots;
`ground_residual()` is measured from the original equations rather than copied
from the input. The dimer's empty excited-root set uses the prepared precision.

## Direct boundary inverse problem

Start from the homogeneous, zero-field rational limit of the boundary scalar
product, Eq. (3.1)–(3.2) in Grijalva, De Nardis and Terras (2019):

- https://arxiv.org/html/1901.10932v4#S3
- https://doi.org/10.21468/SciPostPhys.7.2.023

Use the free-end double-row monodromy built from
$`R(v)=vI+iP`$ with $`v=u-i/2`$. In this convention its creation entry obeys

```math
\mathcal B(-i/2)=2(-1)^{N+1}S_N^-.
```

This is an operator identity, not just its action on Bethe eigenstates.
At $`v=-i`$, $`R(-i)=-2iP_{\rm singlet}`$. Sandwiching the inner double row
between the auxiliary/end-site singlet projectors reduces it to their partial
trace. Repeating the contraction eliminates the interior sites, leaving
the coefficient above. The small-chain ABA oracle tests it on every bit-basis
state at N2 through N5.

Because the creation entries commute, the needed matrix element is therefore
an on-shell/off-shell scalar product:

```math
\langle g|S_N^-|e\rangle
=\frac{(-1)^{N+1}}2
\langle g|\mathcal B(-i/2)\prod_b\mathcal B(v_b)|\uparrow\cdots\uparrow\rangle.
```

Only the singlet bra needs to be on shell. The ket has the triplet's finite
roots plus **one fixed imaginary root**, not an approximate infinite root.
The paper's Hamiltonian and Pauli-spin conventions differ from ours:
$`H_{\rm paper}=4H-(N-1)`$ at the XXX point. The formula below includes the
operator conversion and uses our energy gaps and time units.

## Scaled real determinant

Let $`g=(g_1,\ldots,g_M)`$ be the ground-state rapidities, with $`M=N/2`$,
and $`v=(v_1,\ldots,v_{M-1})`$ the highest-triplet rapidities. Define

```math
\begin{aligned}
q(x,y)&=x^2-y^2,\\{}
U_b&=\prod_{a=1}^M\frac{q(v_b+i,g_a)}{|q(v_b+i,g_a)|},\\{}
\eta_b&=\left(\frac{v_b+i/2}{v_b-i/2}\right)^N,\\{}
K_{ab}&=\frac{2\,\Im\left[\eta_b^* U_b/q(v_b+i,g_a)\right]}{q(v_b,g_a)},
\quad b\lt M,\\{}
K_{aM}&=\frac1{(g_a^2+1/4)^2}.
\end{aligned}
```

In the scalar-product determinant, each finite-root column has a common
unit-modulus phase and a product of $`|q(v_b+i,g_a)|`$. Remove these before
factorization. At the fixed imaginary root, one of the two ABA terms vanishes
exactly and its column reduces to $`K_{aM}`$, with its remaining common product
restored in logarithmic form. The final matrix is real. The factors with
powers proportional to N cancel between the scalar product and the norms,
so none of those large products are evaluated numerically.

For a finite on-shell root set $`r`$, the Gaudin matrix is

```math
\begin{aligned}
G_{aa}(r)&=\frac{2N}{r_a^2+1/4}
-\sum_{c\ne a}\left[\frac2{1+(r_a-r_c)^2}+\frac2{1+(r_a+r_c)^2}\right],\\{}
G_{ab}(r)&=\frac2{1+(r_a-r_b)^2}-\frac2{1+(r_a+r_b)^2},\quad a\ne b.
\end{aligned}
```

Writing $`A(x,y)=[1+(x-y)^2][1+(x+y)^2]`$, the normalized raising weight is

```math
w_e^+=\frac{|\det K|^2}{|\det G(g)|\,|\det G(v)|}
\frac{\prod_a4g_a^2(g_a^2+1/4)^2}{\prod_b(v_b^2+1/4)^2}
\frac{\prod_{a,b}A(g_a,v_b)}
{\prod_{a\lt c}A(g_a,g_c)\prod_{b\lt d}A(v_b,v_d)}.
```

Empty products and the empty determinant are one. In particular, the N2
singlet has $`g_1=1/2`$, $`K_{11}=4`$, $`G_{11}=8`$, giving
$`w^+=1/2`$ and $`C^{zz}(t)=e^{-it}/4`$ exactly.

## Numerical acceptance and validation

Ground and excited states must have the stated sectors and valid quantum
numbers. Each root must be finite, positive, strictly increasing and at the
prepared working precision. Malformed geometry or mixed/exact/unset MPFR
precision throws. Unconverged roots give `roots_unconverged`; rejected measured
residuals, unresolved pivots, near-coincident ground/excited roots and weight
underflow give `precision_limit`. Failed results have no numerical weight.

The recomputed residual uses the solver's max|F|/(2N) units and must be at most
256 N epsilon. Coincident-root rejection uses
256 N epsilon times (1+|u|+|v|). These are working-precision acceptance thresholds,
not mathematical error estimates. No confluent determinant or complex/string
triplet family is implemented here.

Products are accumulated in logarithms. The shared determinant helper applies
column equilibration and ordinary Uni20 LU dispatch, retaining both log magnitude
and sign/phase. Its scaled pivot margin must exceed 128 M epsilon. This is a
rejection diagnostic, not a condition number or certified error bar. A weight
above 1/2 by more than 4096 M epsilon is rejected; accepted values are not clamped.
The stored determinant sign is not, by itself, a physical matrix-element phase.

Independent tests establish the fixed-root identity on every bit-basis state
for N=2 through 5, and compare individual normalized weights to finite-space ABA
vectors and bit-basis exact diagonalization for N=2,4,6,8. The dimer and four-site
chain exhaust the raising norm; larger real-root families omit physical weight.
Native binary64, native extended precision and optional binary128 use these
controls. Explicit MPFR 128/256-bit checks include the dimer, four-site chain,
ownership, stale roots, determinant phases and precision rejection. The small
ABA oracle is a test implementation, not a scalable state solver.

This checkpoint supplies matrix elements. The real-root dynamical scan and
signed endpoint observables are separate consumers; arbitrary bulk-site
matrix elements and complete OBC spectra are outside its coverage.
