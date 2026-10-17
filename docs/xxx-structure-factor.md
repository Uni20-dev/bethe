# XXX two-spinon structure factor on a finite ring

[Overview](../README.md) · [Plotted tutorial](tutorials/xxx-structure-factor.md) ·
[XXX conventions](xxx.md)

For an infinite-chain **spectral density** without broadening, see the separate
[thermodynamic tool](xxx-structure-factor-thermo.md) and its
[heat-map tutorial](tutorials/xxx-structure-factor-thermo.md).

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
fp64, native long-double (also selectable as fp80 when it is extended80), and
optional fp128 arithmetic are supported throughout.

## Parallel calculations

`--threads N` selects a Uni20 oneTBB scheduler with at most N participants;
the default is 1. Root solves and form factors run in bounded batches of
independent states. The ground reference, sorting, moment sums, and output
remain serial. Worker count does not change state IDs, weights, or summation
order at the same precision and with the same executable.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 bethe-xxx-structure-factor 512 \
  --max-candidates 40000 --threads 36 --quiet \
  --csv spectrum.csv --csv-table moments=moments.csv
```

This requests all 32,896 two-spinon states; choose the thread limit for your
machine. CPU time is summed across workers, not elapsed time. Extra workers
also need extra determinant workspaces. Keep numerical-library threading at
one when parallelizing over states to avoid nested oversubscription.

Library callers can pass `bethe::ExecutionOptions{&scheduler, batch_size}` as
the fourth argument of `two_spinon_structure_factor`, or through
`RealExcitationOptions::execution` for shared real-root scans. The scheduler
is borrowed and must outlive the synchronous call; exceptions join active
work before propagating. A null pointer selects Uni20's active scheduler
(serial by default). Bethe does not replace the global scheduler. The default
batch size is 128; changing it does not change the canonical result order.

```cpp
#include <bethe/xxx_structure_factor.hpp>
#include <uni20/async/tbb_scheduler.hpp>

uni20::async::TbbScheduler scheduler(8);
auto spectrum = bethe::heisenberg::two_spinon_structure_factor<double>(
    64, {}, 10000, {&scheduler});
```

Link the scheduler consumer to `TBB::tbb`; Bethe's library target supplies
`uni20_async` but does not impose a particular concurrent backend.

### Example scaling

On Polaron (two Xeon Gold 6254 CPUs, 36 physical cores), the Clang 20 Release
fp64 executable evaluated all 2,080 N=128 lines with these elapsed times:

| Scheduler limit | Elapsed seconds | Speedup |
| --- | ---: | ---: |
| 1 | 14.05 | 1.0× |
| 4 | 3.66 | 3.8× |
| 18 | 0.97 | 14.5× |
| 36 | 0.58 | 24.2× |

These are single-run measurements, including CSV export, with BLAS/OpenMP
thread limits set to one; they are not a portable performance guarantee.
All four runs exported byte-identical spectral rows (metadata and timing
naturally differ). Small systems may not benefit from extra workers.

With the same executable, N=512 completed all 32,896 lines in 332.91 seconds
using 36 workers (10,781.6 process CPU seconds, about 180 MiB peak resident
memory). Its partial spectrum contains 93.41% of the full integrated weight
and 91.12% of the first moment; these fractions are not renormalized.

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
For the same **unscaled** logarithmic equation vector, the chain rule gives
$`G_\lambda=2G_z`$ and $`\det G_\lambda=2^M\det G_z`$.
We divide each regular H column by
$`\prod_a\sqrt{1+(\mu_a-\lambda_b)^2}`$, restore its logarithm in the final
weight, and evaluate column-equilibrated determinants using Uni20's
partial-pivot LU and log-determinant interfaces, with normal backend dispatch
(LAPACK/MPLAPACK when available, otherwise a supported CPU implementation).
Products and determinants are combined in log space.

The pivot diagnostic is `min(abs(diag(U)))/max(1,max(abs(U)))` after column
equilibration; values at or below `128*M*epsilon` are rejected.
This is **not a condition number
or certified error estimate**. Unresolved pivots, near-coincident inter-set
roots, nonfinite results and weight underflow are rejected. The equations are
re-evaluated independently of the stored residual, with a gate of
$`256N\epsilon`$ in normalized residual units. This gate also is not a bound
on form-factor error. Very large rings can require more precision or a more
stable determinant representation; an arbitrarily loose `--tolerance` cannot
bypass these checks. Ground-state validation happens once during preparation;
excited-state validation happens on every evaluation. A false input convergence
flag is rejected before either check.

Validation includes analytic N=2 and N=4 weights, independent bit-basis ED
through N=10 (summed over degenerate subspaces), native-precision coordinate
Bethe-wave overlaps, exact momentum conventions, both sum rules, failure
paths and comparisons between fp64, long-double and fp128.
Prepared/unprepared regression tests cover complete two-spinon families through
N=32, alongside tests of the coordinate Jacobian and snapshot ownership.

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

### Arbitrary precision

With [MPFR/MPC and MPLAPACK MPFR enabled](building.md#enable-arbitrary-precision-for-xxx):

```sh
bethe-xxx-structure-factor 32 --precision bits-256 --threads 4 \
  --json spectrum-256.json
```

Alternatively, `--precision digits-100` requests 100 decimal digits of working
precision (333 significand bits). `bethe-info` reports build capabilities.
The C++ solver takes an explicit finite precision, rather than a global default:

```cpp
using Real = uni20::mpreal;
auto p = uni20::Precision::bits(256);
bethe::SolverOptions<Real> controls(p); // default tolerance = 32*epsilon(p)
// Optional: controls.residual_tolerance = uni20::parse_real<Real>("1e-70", p);
auto ground = bethe::heisenberg::ground_state<Real>(32, controls);
bethe::heisenberg::PreparedRaisingFormFactors<Real> factors(32, ground);
auto excited = bethe::heisenberg::sector_ground_state<Real>(32, uni20::half_int{1}, controls);
auto weight = factors(excited);
auto spectrum = bethe::heisenberg::two_spinon_structure_factor<Real>(32, controls);
```

All roots in a prepared context and its excited states must carry the same
finite precision. Exact, unset or differently rounded roots are rejected,
as are mismatched solver tolerances and warm starts. To change precision,
explicitly convert the initial roots and re-solve before preparing a new context;
padding already-rounded roots with more bits does not restore lost accuracy.
Independent contexts can use different precisions concurrently.

Spectral moment reduction also accepts an explicit finite-precision zero:
`bethe::spectral_moments<Real>(N, lines, Real(0, p))`. This keeps even empty
momentum bins at the selected precision. Tests check analytic small systems
at 128/256 bits, convergence against a 512-bit reference, and serial/parallel
agreement. These tests do not certify arbitrary large-chain form factors.

### Reusing a ground state

For individual form factors, prepare once and reuse the context:

```cpp
#include <bethe/xxx_form_factors.hpp>
namespace xxx = bethe::heisenberg;
auto ground = xxx::ground_state<long double>(64);
xxx::PreparedRaisingFormFactors<long double> const factors(64, ground);
if (!factors.ready()) {
  // Inspect factors.status(), ground_residual() and root_tolerance().
  // Improve the ground solve/precision before calculating weights.
  return;
}
auto excited = xxx::sector_ground_state<long double>(64, uni20::half_int{1});
auto factor = factors(excited);
if (factor.converged()) {
  // *factor.weight = N |<excited|S_0^+|ground>|^2, with no 2*pi.
}
```

The context owns a validated snapshot; it does not borrow `ground`. It stores
standard ABA rapidities (exposed read-only as `aba_rapidities()`), the ground
Gaudin contribution and the ground-only product/column factors. Both input
`RealState` objects still use the solver coordinate **z=2λ**; do not halve their
roots before passing them in. Each excited state gets its own root validation,
Gaudin matrix and transverse determinant.

Const evaluations can run concurrently: scratch buffers belong to each call,
and no lazy mutable cache or global precision setting is involved. The
two-spinon scanner shares one context across its scheduled evaluations.
Copying a context copies its owned data; moving transfers them together with
the source's preparation status. The source of a move is explicitly
non-ready (`precision_limit`), and returns no weights for a valid excited-state
input; assignment can restore it. Self-move assignment is a no-op.

Invalid state shape, sector or quantum-number labels throw
`std::invalid_argument`. Numerical failures are reported as
`roots_unconverged` or `precision_limit`, with absent weights, never a fabricated
zero. `ground_residual()` is the measured normalized residual, not the caller's
stored value; it is absent if validation did not reach a finite residual.
Preparation validates the supplied roots; it does not solve or refine them.

This replaces the former `raising_form_factor(N, ground, excited)` function.
For a one-shot call use `PreparedRaisingFormFactors<Real>(N, ground)(excited)`;
for a scan, keep the context outside the loop. No compatibility wrapper is kept.

### Reproducible timing control

With tests enabled, an optional benchmark compares the prepared evaluator with
the original unprepared formula, retained **only under `tests/`**. From a
configured build directory:

```sh
cmake --build . --target bethe_xxx_form_factor_benchmark
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ./tests/bethe_xxx_form_factor_benchmark 64 fp128 3
```

It solves the states once, verifies every weight in the requested precision,
then times both paths on identical inputs with alternating order. Root-scan and
one-time preparation costs are reported separately; the speedup includes that
preparation cost. This measures the form-factor stage, not a whole-scan speedup
or a determinant-backend change. The benchmark is not built or run by default.

An illustrative N=64 run (528 states, 16 October 2026) on an Intel i7-12700KF,
GCC 13.3 `Release` (`-O3`, no fast-math), with one worker and one BLAS thread:

| Arithmetic | Preparation (ms) | Unprepared (s) | Prepared (s) | Speedup including preparation |
| --- | ---: | ---: | ---: | ---: |
| Binary64 (`fp64`) | 0.080 | 0.06765 | 0.05535 | 1.22× |
| 80-bit extended (`long-double`) | 0.136 | 0.2814 | 0.2243 | 1.25× |
| Binary128 (`fp128`) | 1.442 | 6.955 | 6.236 | 1.12× |

Times are medians of five paired evaluations for binary64/extended precision,
three for binary128, using Uni20 `a25159cc` and the unchanged complete-pivot
determinant helper. Maximum relative weight changes from regrouping the
logarithmic sums were respectively $`7.12\times10^{-14}`$,
$`3.04\times10^{-17}`$ and $`6.88\times10^{-32}`$.
These are local wall-time measurements, not performance guarantees. Root
solving remains dominant: the separately measured binary128 scan took 61.8 s,
so the whole-scan improvement is much smaller than the form-factor-stage gain.

### Module boundaries

`xxx_form_factors.hpp` owns the model-specific normalization;
`xxx_structure_factor.hpp` combines it with the existing scanner.
`dynamical_structure_factor.hpp` contains only model-independent accepted
Lehmann entries and moment sums. The small internal log-determinant helper
does not add a general solver API or a new Uni20 kernel.

Strings, four-spinon states, finite fields, XXZ weights, finite temperature,
OBC form factors and a thermodynamic spectral-density evaluator remain
separate future milestones.
