# A bound triple scattering with a single defect

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Isolated triples](biquadratic-bound-triples.md) · [Two pairs](biquadratic-two-pairs.md)

A three-string plus one real root is another four-defect family, in the
TL module `ell=N-8`. The selected-state library API accepts odd/even N>=8:

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;
auto state = ferro::triple_defect<long double>(128, 125, 121);
if (state.reference.converged) {
    auto gap = state.tl_energy; // E-(N-1), evaluated directly
    auto alpha = state.reference.rapidity;
    auto a = state.reference.center;
    auto L = state.reference.log_deviation;
    auto phi = state.reference.deviation_phase;
}
```

The labels are $`1\le I \le N -3`$ for the real root and $`1\le J \le N -7`$ for the
triple. They are neither momenta nor energy ranks. The corner
$`(I,J)=(N -3,N -7)`$ approaches gap **3 = 2+1**, the separated triple-plus-single
threshold, not the minimum of the full four-defect module. These are TL
singlet insertions, not physical spin flips; no SU(2) decomposition or
spectral weights are implied.

## Command-line selections and scans

```sh
bethe-biquadratic-obc 128 --ferromagnetic --triple-defect 125,121 --roots
bethe-biquadratic-obc 128 --ferromagnetic --triple-defects all --mixed-window 3
bethe-biquadratic-obc 128 --ferromagnetic --triple-defects 4 --mixed-window 3 \
  --precision fp128 --json triples.json --csv triples.csv
```

`--triple-defects COUNT|all` scans $`(N -3)\,(N -7)`$ label pairs unless
`--mixed-window WIDTH` selects only the highest WIDTH values of I and J.
For this family `1<=WIDTH<=N-7`, giving WIDTH^2 candidates. The default
`--max-candidates 10000` bounds all solves, not just retained output rows;
overflow-safe counts are checked before allocation. `--real-defects` is
not applicable: this family has exactly one real root alongside its triple.

The shared cluster reporter sorts by direct gap, with ties ordered by
`(I,J)`, and retains the lowest COUNT converged candidates. Tables are
`states`, `reference`, `labels` (I, J, alpha), and `string` (a, L, phi),
plus optional `roots`. The four root indices are 0 for the central triple
root, 1 and 2 for its conjugate pair, and 3 for the additional real root.
The usual JSON/CSV/TSV exports, streaming, and no-retain output are available.

If any candidate fails, exit status is 2 and the summary remains partial.
One failed diagnostic row follows the retained converged levels; its
verified `gap` is null, while energy and `tl_energy` are only estimates.
Metadata records all scanned, converged, failed, and retained counts.
Neither `all` nor a successful scan certifies a complete module spectrum.

## Reusing the isolated triple

The reference is still quantum-group XXZ at Delta=3/2 with opposite end
fields, not the zero-end-field chain. Let eta=acosh(Delta). Its roots are

```math
\begin{aligned}
u_0&=\frac{ia}{2},\qquad u_+=\eta+\frac{ia}{2}+z,\qquad u_-=\overline{u_+},\\{}
v&=\frac{i\alpha}{2},\qquad z=e^{-L+i\phi}.
\end{aligned}
```

We retain the isolated triple's three regularized equations and analytic
Jacobian. Only the external scattering is new. Define the logarithm of
the direct and reflected factors in
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020):

```math
\begin{aligned}
T(u,v)={}&\log\sinh(u-v+\eta)-\log\sinh(u-v-\eta)\\{}
&+\log\sinh(u+v+\eta)-\log\sinh(u+v-\eta).
\end{aligned}
```

The shared `detail/open_string_scattering.hpp` evaluates T and its two
analytic complex derivatives. It is only for nonsingular external
factors; the triple's singular internal factor is still evaluated using
L and phi, never by subtracting rounded roots.

To choose the continuous product-phase branch, put

```math
\begin{aligned}
F(a,\alpha)&=\sum_{s=\pm1}[\Theta(a+s\alpha;\eta)+\Theta(a+s\alpha;2\eta)],\\{}
\mathrm{lift}(\mathrm{raw},\mathrm{target})&=\mathrm{target}+\mathrm{wrap}(\mathrm{raw}-\mathrm{target}).
\end{aligned}
```

Theta and wrap have the same definitions as in the isolated-triple guide.
Our coupled residuals subtract

```math
\begin{aligned}
&\frac{\mathrm{lift}(\mathrm{Im}[T(u_0,v)+2T(u_+,v)],F(a,\alpha))}{2N}
&&\text{from the product phase},\\{}
&\frac{\mathrm{Re}T(u_+,v)}{2N}
&&\text{from the outer-root modulus},\\{}
&\mathrm{Im}T(u_+,v)
&&\text{inside the wrapped outer-root phase}.
\end{aligned}
```

The fourth equation is

```math
\Theta(\alpha;\eta/2)
-\frac{2\pi I+\mathrm{lift}(\mathrm{Im}[T(v,u_0)+T(v,u_+)+T(v,u_-)],F(\alpha,a))}{2N}=0.
```

The ideal fused phase fixes a multiple of 2pi only; it does not replace
finite-deviation scattering. Its derivative cancels inside/outside the
local lift, so the Jacobian uses derivatives of the original factors.
Both sides of $`\alpha =a`$ are tested. The actual coincident-root point is
excluded, and the inherited triple domain requires $`\lvert z \rvert \lt \eta /4`$.

A linear initial triple center and a bare-phase real-root seed start the
solve. During initialization, the two center equations use the analytically
fused phase F, while L and phi use the isolated triple equations. This avoids
the removable external-factor poles at $`\alpha =a`$ in the ideal limit.
The finite stage then restores **all** coupled equations above; convergence
is never accepted using the initializer alone. The common Newton driver
uses one iteration budget and the requested native precision throughout.
It solves four real unknowns, independent of N. Failed iterates retain
consistent diagnostics and estimates, not verified energies.

Once exp(-L) underflows, the two internal deviation equations are linear
in L and phi. After a finite-stage Newton trial rounds the angles, those
two rows are solved at the represented angles. This prevents sub-ulp
angular steps near coincident centers from inducing a spurious modulus
stall. It is part of the trial normalization and line search, not extra
iterations, a tolerance relaxation, or a change of precision. It is never
used during ideal initialization.

The lower API is
`bethe::xxz::quantum_group::triple_defect::solve(N,Delta,I,J,options)`.
Finite Delta>1 is accepted, but neither convergence nor branch existence
is guaranteed for every parameter or near unresolved root collisions.

## Coverage and validation

At Delta=3/2, combining the four-real-root, pair-plus-two-real-root,
two-pair, and triple-plus-real-root families matches independent module ED
as a multiset at N=8,9,10. The last family contributes $`(N -3)\,(N -7)`$ levels:

| N | Previously covered | Triple + real root | Module dimension | Remaining |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 8 | 5 | 14 | 1 |
| 9 | 28 | 12 | 42 | 2 |
| 10 | 66 | 21 | 90 | 3 |

The remaining counts are filled by the separately implemented
[four-string droplet family](biquadratic-bound-quartets.md). The combined
multiset now matches the complete small-chain module spectra, but this
counting is not a general completeness proof.

Additional tests check the original outer-, central-, and real-root
equations, every Jacobian column, iteration-budget failures, and the
threshold through N=100000 in fp64, long double, and enabled fp128. L and
phi remain meaningful even when exp(-L) underflows.
All selected labels at N=8,...,10 also match distinct ED levels at
Delta=1.25, 2 and 3. The initializer's Jacobian is checked even at coincident
ideal centers, while the final solver rejects coincident physical roots.
