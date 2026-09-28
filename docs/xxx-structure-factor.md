# XXX two-spinon structure factor on a finite ring

[Overview](../README.md) · [Plotted tutorial](tutorials/xxx-structure-factor.md) ·
[XXX conventions](xxx.md)

`bethe-xxx-structure-factor` adds **intensities**, not just energies, to the
existing real-root two-spinon triplet family. It supports even periodic chains
at zero field, spin-1/2 operators, and $`H=\sum_j\mathbf S_j\cdot\mathbf S_{j+1}`$
with $`J=1`$. It does **not** compute the full dynamical structure factor.

```sh
bethe-xxx-structure-factor 64 --csv spectrum.csv
bethe-xxx-structure-factor 64 --channel raising --precision long-double \
  --json spectrum.json --csv-table moments=moments.csv
bethe-xxx-structure-factor 8 --diagnostics --roots
```

The complete supported family is always requested: no energy truncation or
`--excitations all` is needed. There are $`N(N+2)/8`$ candidates: 528 at N=64.
`--max-candidates` defaults to 10000 and rejects a larger family before solving.
`--max-iterations` and `--tolerance` control the existing root solver.
fp64, native long-double, and optional fp128 arithmetic are supported throughout.

## Operator, momentum and spectral normalization

Define the raising operator and its spectral measure by

```math
S_q^+=\frac1{\sqrt N}\sum_{j=0}^{N-1}e^{-iqj}S_j^+,
\qquad
S^{-+}_2(q,\omega)=2\pi\sum_{n\in 2\mathrm{sp}}
 w_n\,\delta(\omega-E_n+E_0),
\qquad
w_n=|\langle n|S_q^+|0\rangle|^2.
```

The raising channel is named `raising`, corresponding to $`S^{-+}`$, not
$`S^{+-}`$. At zero field the two transverse channels coincide. The default
`--channel zz` uses SU(2) invariance of the singlet:

```math
S^{zz}_2=\frac12 S^{-+}_2.
```

There is no need to construct an infinite-rapidity descendant in the
$`S^z=0`$ sector. `state_id` still identifies the computed highest-weight
$`S=S^z=1`$ representative, even for `zz` output. Do not multiply its intensity
by the triplet degeneracy.

Our state momentum obeys $`P=\pi M-(2\pi/N)\sum I_j`$. Translation moving
site $`j`$ to $`j+1`$ has eigenvalue $`e^{-iP}`$, hence the above Fourier sign
selects $`q=P_0-P_n\pmod{2\pi}`$. This subtraction matters on N=4m+2 rings,
where $`P_0=\pi`$. The `q` column is **momentum transfer**, not the excited
state's absolute momentum. Selection uses exact integer momentum indices.
At the selected momentum, $`w_n=N|\langle n|S_0^+|0\rangle|^2`$.

`gap` is $`E_n-E_0`$ for the same finite ring. `weight` is the coefficient
before the delta function **without** the $`2\pi`$; it is not a sampled
spectral density. No broadening, bin widths or sum-rule rescaling are applied.

## Tables and failure semantics

| Table | Contents | Screen selection |
| --- | --- | --- |
| `spectrum` (primary) | `state_id`, `momentum_index`, `q`, `gap`, `weight`; accepted lines only | Default |
| `moments` | Partial weight and first moment per momentum, full first-moment sum rule and fraction | Default |
| `diagnostics` | Root residuals, iterations, convergence, form-factor status and pivot indicator | `--diagnostics` |
| `roots` | Ground state (`state_id=0`) and retained excited-state Bethe numbers and $`z=2\lambda`$ roots | `--roots` |

Named CSV/TSV exports request these tables independently of screen flags:

```sh
bethe-xxx-structure-factor 16 --csv-table roots=roots.csv \
  --tsv-table diagnostics=diagnostics.tsv --csv-table spectrum=lines.csv
```

IDs are local to one run: near-degenerate energy ordering can change with
precision. Quantum numbers identify states across runs. Root failures are
counted in metadata; only the first failed root state is retained as an
example in `diagnostics`, with no `state_id`. Ground-root iterates carry their
convergence flag. `iterate_energy` is diagnostic even when roots failed.

Exit 0 means the **supported two-spinon family** succeeded, not that the full
DSF was computed. Exit 2 means root or form-factor failures; no failed weight
is replaced by zero. The moments then sum only accepted lines. `precision_limit`
is distinct from `roots_unconverged`; try higher precision and inspect the
diagnostics. Invalid inputs exit 1 before export files are opened.

The exact zero-momentum moment is zero; its first-moment fraction is missing
rather than the indeterminate ratio 0/0.

## Sum rules: measure missing weight, do not hide it

For the full raising-channel spectrum, with **unshifted total** ground energy
$`E_0`$,

```math
\frac1N\sum_q\sum_n w_n(q)=\frac12,
\qquad
\sum_n(E_n-E_0)w_n(q)=-\frac{4E_0}{3N}(1-\cos q).
```

Both right-hand sides are halved for `zz`. We report the computed partial
moments and their fractions of these full sum rules. The integrated first
moment divides $`\sum_{q,n}(E_n-E_0)w_n`$ by $`-4E_0/3`$ in the raising channel.
At N=64 the supported family captures about 97.71% of the integrated weight
and 96.22% of the first moment. These are **finite-size results**. In particular,
the familiar thermodynamic two-spinon fractions 72.89% and 71.30% are not
finite-ring normalization targets. Never rescale the exported weights to
either those fractions or 100%.

The sum rules and thermodynamic benchmarks are discussed in
[Caux, Mossel and Pérez Castillo, Eqs. (62)–(63)](https://arxiv.org/abs/0806.3069).

## Determinant formula and implementation

The implementation takes the rational limit of the normalized transverse
formula in [Caux, Hagemans and Maillet, Eqs. (11)–(13)](https://arxiv.org/abs/cond-mat/0506698),
based on [Kitanine, Maillet and Terras](https://arxiv.org/abs/math-ph/9807020).
Convert the solver's roots once: $`\lambda_{\mathrm{ABA}}=z/2`$.
Below $`\mu`$ denotes the M=N/2 ground roots and $`\lambda`$ the M−1 excited
roots, all in this standard ABA convention. Hermitian conjugation lets us
evaluate the lowering matrix element in the opposite direction.

For any root set $`x`$, define the Gaudin matrix

```math
\Phi_{ab}(x)=\delta_{ab}\left[
\frac{N}{x_a^2+1/4}-\sum_{k\ne a}\frac{2}{1+(x_a-x_k)^2}\right]
+(1-\delta_{ab})\frac{2}{1+(x_a-x_b)^2}.
```

The M by M transverse matrix has entries

```math
\begin{aligned}
H_{ab}&=\frac{\prod_{j\ne a}(\mu_j-\lambda_b-i)
-r_b\prod_{j\ne a}(\mu_j-\lambda_b+i)}{\mu_a-\lambda_b},
\quad b=1,\ldots,M-1,\\{}
r_b&=\left(\frac{\lambda_b+i/2}{\lambda_b-i/2}\right)^N,
\qquad H_{aM}=\frac1{\mu_a^2+1/4}.
\end{aligned}
```

The Fourier weight is

```math
w=N\,
\frac{\prod_{a=1}^M(\mu_a^2+1/4)}{\prod_{b=1}^{M-1}(\lambda_b^2+1/4)}
\frac{|\det H|^2}{|\det\Phi(\mu)\det\Phi(\lambda)|}
\frac{1}{\prod_{a\gt b}[1+(\mu_a-\mu_b)^2]
\prod_{a\gt b}[1+(\lambda_a-\lambda_b)^2]}.
```

Thus a raw determinant of the solver's normalized z-coordinate Jacobian is
**not** the norm needed here. Coordinate scale and ABA product factors matter.
We divide each regular H column by
$`\prod_a\sqrt{1+(\mu_a-\lambda_b)^2}`$, restore its logarithm in the final
weight, and evaluate equilibrated determinants with complete-pivot LU.
Products and determinants are combined in log space.

The smallest scaled pivot is a rejection diagnostic, **not a condition number
or certified error estimate**. Unresolved pivots, near-coincident inter-set
roots, nonfinite results and weight underflow are rejected. The equations are
re-evaluated independently of the input convergence flag, with a gate of
$`256N\epsilon`$ in normalized residual units. This gate also is not a bound
on form-factor error. Very large rings can require more precision or a more
stable determinant representation; an arbitrarily loose `--tolerance` cannot
bypass these checks.

Validation includes analytic N=2 and N=4 weights, independent bit-basis ED
through N=10 (summed over degenerate subspaces), native-precision coordinate
Bethe-wave overlaps, exact momentum conventions, both sum rules, failure
paths and comparisons between fp64, long-double and fp128.

## C++ entry points and scope

```cpp
#include <bethe/xxx_structure_factor.hpp>
auto result = bethe::heisenberg::two_spinon_structure_factor<long double>(64);
// result.lines are raising-channel LehmannLine<long double> entries.
// Divide weights and moments by two for Szz; gaps/momenta are unchanged.
if (!result.converged()) {
  // Inspect result.scan, result.form_factors and missing candidate counts.
}
```

`xxx_form_factors.hpp` owns the model-specific normalization;
`xxx_structure_factor.hpp` combines it with the existing scanner.
`dynamical_structure_factor.hpp` contains only model-independent accepted
Lehmann entries and moment sums. The small internal log-determinant helper
does not add a general solver API or a new Uni20 kernel.

Strings, four-spinon states, finite fields, XXZ weights, finite temperature,
OBC form factors and a thermodynamic spectral-density evaluator remain
separate future milestones. This is not an ABACUS-style adaptive state search.
