# Magnetization after a boundary spin flip

This experimental API calculates the time-dependent magnetization at **both endpoints**
after applying $`S_1^+`$ to the singlet ground state of an even, free-end XXX
chain. The initial implementation retains all states in the existing finite
real-root triplet family. It is not the full excitation spectrum, and it does
not yet give the spatial profile at interior sites.

This is a different observable from the
[boundary autocorrelation](experimental-xxx-boundary-form-factors.md). It needs signed,
consistently phased matrix elements **between excited states**, rather than
just positive ground-to-excited-state weights.

## State and observable

Use our usual spin and energy conventions:

```math
H=\sum_{j=1}^{N-1}\mathbf S_j\cdot\mathbf S_{j+1},\qquad
S^\alpha=\sigma^\alpha/2,\qquad J=1.
```

The exact, normalized perturbed state and its local magnetization are

```math
|\psi(t)\rangle=\sqrt2 e^{-i(H-E_g)t}S_1^+|g\rangle,\qquad
m_j(t)=\langle\psi(t)|S_j^z|\psi(t)\rangle.
```

The factor $`\sqrt2`$ follows from
$`\langle g|S_1^-S_1^+|g\rangle=1/2`$. The full state has total magnetization
one and $`m_1(0)=1/2`$. For a different readout site,
$`m_j(0)=-2\langle g|S_1^zS_j^z|g\rangle`$. The opposite endpoint is therefore
not initially unmagnetized: the local action conditions an entangled ground
state. The singlet itself has zero local magnetization everywhere.

Applying $`S_1^z`$ instead produces a spin-flip-symmetric state with zero local
magnetization at every time. That state still has nontrivial correlations; it
is not the spin-raising calculation implemented here.

For retained eigenstates A, choose a gauge in which
$`a_n=\langle n;1,1|S_1^+|g\rangle`$ is real and positive. Then

```math
\begin{aligned}
|\psi_A(t)\rangle&=\sqrt2\sum_{n\in A}a_n e^{-i(E_n-E_g)t}|n\rangle,\\{}
s:=\langle\psi_A|\psi_A\rangle&=2\sum_A a_n^2=4\sum_A w_n^{zz},\\{}
O_{mn}&=\langle m;1,1|S_1^z|n;1,1\rangle,\\{}
m_{1,A}(t)&=2\sum_n a_n^2O_{nn}
+4\sum_{m\lt n}a_ma_nO_{mn}\cos[(E_n-E_m)t].
\end{aligned}
```

We leave this projection **unnormalized**: we do not divide by s or force its
initial endpoint value to 1/2. Its total magnetization is s. The cosine
coefficients are signed interference terms, not probabilities or DSF weights.
The expectation is real and even in time; there is no imaginary component to
plot, although the underlying state has complex time-dependent amplitudes.

Reflection gives $`\langle m|S_N^z|n\rangle=r_mr_nO_{mn}`$. Thus the same
matrix calculation supplies the opposite endpoint: multiply each cosine
coefficient by the product of the two
[reflection eigenvalues](experimental-xxx-boundary-dynamics.md).

## Phased determinant in standard rapidities

Our derivation starts from the finite-size boundary scalar product and local
spin reconstruction in Grijalva, De Nardis and Terras, *Open XXZ chain and
boundary modes at zero temperature*, SciPost Phys. **7**, 023 (2019), §3:

- https://arxiv.org/html/1901.10932v4#S3
- https://doi.org/10.21468/SciPostPhys.7.2.023

The following is our homogeneous, zero-field **rational specialization**, in
spin rather than Pauli normalization. The paper's Hamiltonian is
$`4H-(N-1)`$ at the XXX point; we use our own energies and time units throughout.
Let X and Y be sorted positive on-shell root sets with
$`M=N/2-1`$ roots, in standard rapidities $`u=z/2`$. Define

```math
\begin{aligned}
p(u)&=\frac1{u^2+1/4},\qquad k(u)=\frac2{1+u^2},\qquad q(u,v)=u^2-v^2,\\{}
Q&=\frac{\prod_{x\in X}p(x)}{\prod_{y\in Y}p(y)},\qquad
d_p=\sum_{y\in Y}p(y)-\sum_{x\in X}p(x),\\{}
A_j&=\frac{\prod_{\ell\ne j}q(y_j,y_\ell)}{\prod_\ell q(y_j,x_\ell)},\qquad
B_j=\prod_\ell\frac{q(y_j-i,x_\ell)}{q(y_j-i,y_\ell)}.
\end{aligned}
```

Using the on-shell equations for Y removes the complex phases from the
transformed matrix:

```math
\begin{aligned}
\overline M_{jk}
&=4y_j A_j\Im B_j\,\delta_{jk}-k(y_j-y_k)+k(y_j+y_k),\\{}
\overline P_{jk}&=4y_jy_k p(y_j)Q[d_p-p(y_j)],\qquad
D(X,Y)=\overline M-2\overline P.
\end{aligned}
```

Let $`G(X)`$ be the ordinary open-chain Gaudin matrix, defined in the
[raising-form-factor derivation](experimental-xxx-boundary-form-factors.md),
and define the positive root-only factor

```math
S_X=\prod_{x\in X}4x^2
\prod_{a\lt b}\frac{q(x_a,x_b)^2}
{[1+(x_a-x_b)^2][1+(x_a+x_b)^2]}.
```

Write $`t_X=\mathrm{sgn}(\det K_X)`$ for the **real raising determinant** of
that state and $`r_X`$ for its reflection eigenvalue. In the positive-coupling
gauge the signed, normalized spin matrix element is

```math
O_{XY}=(-1)^M r_Xr_Yt_Xt_Y
\frac{\det D(X,Y)}{2\sqrt{\det G(X)\det G(Y)}}
\sqrt{\frac{S_X}{S_Y}}
\frac{\prod_{x\in X}(x^2+1/4)}{\prod_{y\in Y}(y^2+1/4)}.
```

The last, **directed product factor is $`1/Q`$**. Retain it from the original
operator formula, Eq. (3.8), through the determinant transformation. It cancels
in the product of forward and reverse matrix elements, so tests of squared
weights alone cannot detect its omission. Dropping it in our rational formula
breaks Hermiticity and gives wrong interference magnitudes even at N4. This is
a normalization safeguard in our derivation, not a claim about the paper's
subsequent thermodynamic results.

The signs $`t_X`$ alone are not physical amplitudes. They must be combined with
the phase conversion and reflection factors above. We test the invariant
cyclic products
$`\langle g|S_1^-|m\rangle O_{mn}\langle n|S_1^+|g\rangle`$
against independently constructed finite-Hilbert-space Bethe vectors.
Taking arbitrary square roots of the old positive weights is insufficient.

For identical states, the off-diagonal expression has removable poles. Its
analytic coincident-state limit is

```math
O_{XX}=\frac12-2u^{\mathsf T}G^{-1}b,\qquad
b_i=2u_i p(u_i)^2.
```

The code solves $`Gx=b`$ and contracts $`u^{\mathsf T}x`$; it does not form an
inverse, displace the roots or evaluate a singular off-diagonal determinant.
The empty-root fully polarized state has $`O_{XX}=1/2`$.

## Truncation bound

Let $`\epsilon=1-s`$ be the exact omitted norm of the normalized raised state.
Its full and projected density operators obey

```math
\left\|\,|\psi\rangle\langle\psi|
-|\psi_A\rangle\langle\psi_A|\,\right\|_1
=\sqrt{\epsilon^2+4s\epsilon}
=\sqrt{\epsilon(4-3\epsilon)}.
```

This follows by diagonalizing their difference in the two-dimensional span
of the retained and omitted components. Since $`\|S_j^z\|=1/2`$,

```math
|m_j(t)-m_{j,A}(t)|\le\frac12\sqrt{\epsilon(4-3\epsilon)}.
```

This uniform absolute bound includes interference with omitted states and
scales as $`\sqrt\epsilon`$, not $`\epsilon`$. The stronger missing-weight
bound for a two-point autocorrelation cannot be reused here. A high retained
norm consequently does not guarantee that every small magnetization feature
is accurate. Arithmetic/root errors are separate from this physical bound.
Negative omitted norm is reported as unresolved, never silently clamped.

## Build and run

This is an opt-in build-tree component, with the same
[experimental build requirements](experimental.md) as the boundary correlator.
From a new build directory, configure the source directory and build the target:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_APPS=ON
cmake --build . --target bethe-experimental-xxx-boundary-magnetization -j2
./experimental/bethe-experimental-xxx-boundary-magnetization 16 --threads 2 --csv time.csv
./experimental/bethe-experimental-xxx-boundary-magnetization 8 --precision bits-160 --json result.json
```

The default time grid has 101 points from zero to 100. Set `--time-max` and
`--time-points` to change it. `--precision fp64`, `long-double`, `fp128` and
`bits-N` use the configured scalar/provider capabilities; runtime precision
requires an MPFR-enabled build. The preamble records the working precision,
root controls, family, pair budget, CPU time and elapsed time.

Tables are declared in `--help`:

| Table | Content | Screen option |
|---|---|---|
| `time` | Both real endpoint magnetizations | Default, primary CSV/TSV table |
| `summary` | Projected norm, omitted norm and uniform bound | Default |
| `states` | Every canonical state, amplitude and preparation diagnostics | `--states` |
| `pairs` | Every upper-triangle pair, signed element and cosine term | `--pairs` |
| `roots` | Ground and attempted excited roots in solver coordinates | `--roots` |

For example, `--csv-table pairs=pairs.csv` exports all pair records without
showing them on screen. `--json result.json` exports all declared tables.
Unknown table names are errors. A missing scalar is absent/null, never a
successful zero. Failed preparations retain every requested state and pair;
an unresolved pair produces **no time curve**, and the executable returns 2.

State IDs are one-based indices in the lexicographic quantum-number order,
with explicit candidate strings. They are not momenta or energy ranks.
Pair frequencies are the absolute energy differences: cosine is even, so
this preserves the phase convention without sorting the states by energy.

## Callable interface and execution

```cpp
#include <bethe/experimental/xxx_boundary_magnetization.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace boundary = bethe::experimental::xxx::open;
uni20::async::TbbScheduler scheduler(2);
auto pending = boundary::boundary_real_triplet_magnetization_async<double>(
    8, {}, 250000, {&scheduler});
auto const& result = pending.get_wait(scheduler);
if (!result.converged()) throw std::runtime_error("unresolved boundary pairs");
auto left = result.at(10.0);
auto right = result.at(10.0, true);
```

The synchronous `boundary_real_triplet_magnetization` wrapper awaits the same
final result. The Async task owns its spectrum dependency, state preparations
and pair records; the caller keeps the scheduler alive until completion.
Pair workers draw continuously from a bounded queue, without fixed batch
barriers. Active work drains before a cancellation or exception is delivered.
Dropping a result handle does not request cancellation.

`BoundarySpinFormFactors<Real>` also provides individual same-sector elements.
Its `prepare(state)` rechecks the original equations and computes the raising
gauge, cached products and analytic diagonal. Prepared states are immutable;
copies share owned data and moves invalidate the source. Only their original
calculator, or a copy of that calculator, can combine them. This prevents
mixing ground preparations, phases or working precisions. Runtime roots and
time arguments must use the same explicit precision.

Every solve and phased log determinant uses normal Uni20 dispatch. The
Gaudin diagonal response uses the shared square-solve helper and an independent
backward-residual check. Positive row scales are removed before factorization
and restored logarithmically. Unresolved root coincidences, pivots and
underflow are failures, not exact zeros. The empty-root dimer is handled by
the analytic diagonal, without a zero-dimensional linear solve.

The upper-triangle pair budget is checked before any roots or quadratic
storage. N64 has 528 states and 139656 pairs; the default 250000-pair budget
rejects N256. The spectrum and pair stages expose separate parallelism.
State preparation currently re-evaluates the raising normalization to retain
its independent validation contract; future reuse must preserve that check.

## Validation and limits

Focused tests compare the gauge-invariant signed cyclic products and both
endpoint curves with independently constructed finite-Hilbert-space Bethe
vectors through N8. They check both directions of every matrix element,
reflection signs, total projected spin, and serial/concurrent identity.
Small-chain ED checks the complete evolution and the omitted-state bound.
The dimer has constant magnetization 1/2 at both ends; N4 recovers the full
raised state. Native binary64/binary80/binary128 and MPFR 128/256-bit tests
use the original precision throughout the determinant and wavefunction checks.

This is a restricted real-root-family calculation, not arbitrary bulk-site
OBC dynamics or a complete finite-chain spectrum. The physical omitted-norm
bound is independent of any arithmetic-error estimate. A rounded norm slightly
above one leaves that bound absent; no clipping or amplitude fitting is done.
The norm acceptance allowance is a working-precision diagnostic, not a
certified error bar. Exact-zero off-diagonal determinants and cross-state root
coincidences need separate limiting representations before they can be accepted.
