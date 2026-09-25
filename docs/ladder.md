# An integrable two-leg spin ladder

[Overview](../README.md) · [SU(3) chain](su3.md) · [Model catalogue](models.md)

`bethe-ladder-pbc` implements Wang's periodic ladder, with **a particular
four-spin interaction required for integrability**. It is not a solver for
an ordinary two-leg Heisenberg ladder at arbitrary couplings. There are L
rungs, two spin-1/2 operators S and T per rung, and longitudinal field h:

```math
\begin{aligned}
H={}&\sum_j\left[\mathbf S_j\cdot\mathbf S_{j+1}+\mathbf T_j\cdot\mathbf T_{j+1}
+4(\mathbf S_j\cdot\mathbf S_{j+1})(\mathbf T_j\cdot\mathbf T_{j+1})\right]\\{}
&+J_r\sum_j\mathbf S_j\cdot\mathbf T_j-h\sum_j(S_j^z+T_j^z).
\end{aligned}
```

The leg coefficient is 1, the four-spin coefficient is 4, and $`J_{r}`$ can
have either sign. The field h has energy units (g*mu_B absorbed) and defaults
to zero. Indices are periodic; L>=2, including odd lengths. At
L=2 the closing bond is counted separately, as in our other periodic chains.

## First calculations

```sh
build/bethe-ladder-pbc 12 --rung 1
build/bethe-ladder-pbc 12 --rung -1 --sectors
build/bethe-ladder-pbc 6 --rung 0 --singlets 4 --roots
build/bethe-ladder-pbc 24 --rung 2 --precision long-double
build/bethe-ladder-pbc 12 --rung 5 --field 3 --sectors
build/bethe-ladder-pbc 6 --rung 0 --singlets 4 --field -0.125 --roots
build/bethe-ladder-pbc 12 --rung 1 --sz 2 --sectors
build/bethe-ladder-pbc 6 --rung 0 --sz 0 --singlets 4 --roots
# Requires a binary128-enabled build:
build/bethe-ladder-pbc 8 --rung 1 --precision fp128 --roots
```

The default compares all singlet-count sectors. `--singlets NS` selects a
single sector, while `--sectors` lists all L+1 sector minima and reports the
lowest one in the overview. These two options are mutually exclusive.
Within each singlet sector, the triplet populations are minimized too.
`--sz M` additionally fixes total physical Sz to the integer M in $`[-L,L]`$;
without it, Sz is minimized as well. The implementation returns one
minimizing representative, not every degenerate state. Reported momenta
refer to translation by one rung, `P=2*pi*momentum_index/L`; the reflected
index identifies its parity partner (possibly the same momentum).

The `magnetization` column gives total `M=N_t+ - N_t-`, not M per rung.
At a level crossing only one minimizer is returned, not all degenerate states
or a thermal average. Without `--sz`, exactly at zero field the original balanced-triplet
representative is retained; its M is not a unique zero-field response.

Excited-state enumeration, twists, open ends, arbitrary
four-spin couplings, correlation functions and thermodynamics are not
implemented.

## Fixed magnetization sectors

`--sz M` can be combined with `--singlets NS` or `--sectors`, and with
either sign of `--field`. At fixed M the field contributes only `-h*M`:
it does not change the minimizing roots or singlet count. This differs
from choosing a field and letting M minimize freely; some finite-size
sectors need not become global minima for any field.

The rung contains two spin-1/2 sites, so total Sz is always an integer,
even for odd L. There is no extra parity constraint on $`L -N_{s} -M`$, because
the triplet t0 has zero projection. Feasibility is $`\lvert M \rvert \le L -N_{s}`$.
With `--sectors`, only `N_s=0,...,L-|M|` are listed. Incompatible explicit
constraints are input errors, not empty successful calculations.
At $`\lvert M \rvert =L`$ the unique polarized state is returned analytically, even when
it is *not* the unconstrained ground state and even with zero solve budgets.

These results are minima of **Sz sectors**, not states of prescribed total
spin S, and do not enumerate excited levels within a sector. At h=0 a
spin-S multiplet contributes to all $`\lvert M \rvert \le S`$; equal sector energies can
therefore be ordinary SU(2) degeneracies.

## Why a ladder becomes a four-color chain

Each rung has one singlet s and three triplets t+, t0, t-. Their rung
energies are $`-3J_{r} /4`$, $`J_{r} /4-h`$, $`J_{r} /4`$, $`J_{r} /4+h`$, respectively. The identity

```math
P_{\mathrm{rung}}(j,k)=\left(2\mathbf S_j\cdot\mathbf S_k+\frac12\right)
\left(2\mathbf T_j\cdot\mathbf T_k+\frac12\right).
```

turns the leg/four-spin part into $`P_{\mathrm{rung}} -1/4`$. Thus

```math
\begin{aligned}
H&=\sum_jP_{\mathrm{rung}}(j,j+1)-\frac L4+J_r(L/4-N_s)-hM,\\{}
E&=E_{\mathrm{perm}}-\frac L4+J_r(L/4-N_s)-hM.
\end{aligned}
```

The singlet count N_s is conserved. Only a chemical-potential term changes
as J_r and h vary; a particular multiplet branch has roots independent
of both, although the minimizing multiplet can change. The field is a linear
combination of conserved color populations, so the same Bethe equations apply.
We use the permutation form in
[Wang's Eqs. (2)–(5)](https://arxiv.org/html/cond-mat/9901168), with
$`J_{r} =2J`$ in those equations and the energy constants restored as above.

Two useful exact checks follow. The all-singlet product has
$`E =3L (1-J_{r})/4`$. It is a global ground state for $`J_{r} -\lvert h \rvert \ge 4`$, so the default
calculation returns it analytically, even with zero iteration/branch
budgets. This sufficient condition follows from the operator inequality
`sum(P-1)>=-4*N_triplet`; it does not assert uniqueness on the boundary or locate
every finite-ring crossing. The lowest one-triplet energy above the product is
`J_r-|h|-4*sin²(pi*floor(L/2)/L)`. The fully polarized triplet product is
also analytic when $`\lvert h \rvert \ge 4+\max (J_{r},0)`$, with
$`E =L \,(3/4+J_{r} /4-\lvert h \rvert)`$ and `M=sign(h)*L`. The same bound applies to defects
of the favored triplet color: every other rung color then costs at least 4.
These are sufficient conditions, not all finite-ring phase boundaries.
In the N_s=0, h=0 sector the model reduces to an SU(3) triplet permutation
chain plus the constant $`L \,(J_{r} -1)/4`$.

## Finite rings: highest weights are not physical populations

An SU(4) highest weight is a Young diagram with row lengths
$`\mu_{1} \ge \mu_{2} \ge \mu_{3} \ge \mu_{4} \ge 0`$, summing to L. A multiplet contains a population
vector w exactly when mu dominates its sorted entries: every prefix sum
of mu is at least the corresponding prefix sum of sorted(w).

This matters on a ring. It is tempting to solve just the packed real sea
for the physical populations and call it the sector minimum. That is
wrong even for six rungs:

| Physical populations (s,t+,t0,t-) | Own highest-weight sea E_perm | Actual minimizing highest weight | Minimum E_perm |
| --- | --- | --- | --- |
| 4,1,1,0 | -1 | 4,2,0,0 | $`1-\sqrt{5}`$ |
| 2,2,1,1 | $`-(5+\sqrt{17})/2`$ | 2,2,2,0 | $`-1-\sqrt{13}`$ |

In both cases the lower state is an SU(4) **descendant**. Lowering operators
change its populations without changing the permutation energy or momentum;
the physical rung contribution still uses the requested N_s. The report
therefore keeps physical populations, highest-weight rows, and the
descendant flag separate. `--roots` prints the finite roots of that
highest-weight representative, **not fabricated finite roots for its
descendant**. With `--sectors --roots`, only the selected best state's
roots are printed.

At h=0, for N_t=L-N_s triplets, it suffices to use the population representative
$`(N_{s},q +[r \gt 0],q +[r \gt 1],q)`$, where $`N_{t} =3q +r`$. Every SU(3) triplet multiplet
contains this balanced weight. Scanning all compatible SU(4) highest
weights therefore retains every triplet-multiplet choice at fixed N_s.
This weight-space statement should not be confused with assuming a
generic ordering of periodic multiplet energies:
[Hakobyan's periodic ordering theorem](https://arxiv.org/abs/cond-mat/0403587)
has an additional row-parity condition.

At nonzero h we instead minimize the Zeeman term within each multiplet.
Writing its rows as `lambda_0,...,lambda_3`, restriction to the three triplet
colors gives interlacing diagrams $`\lambda_{i} \ge \mu_{i} \ge \lambda _{i +1}`$ with
`sum(mu)=L-N_s`. A singlet count occurs precisely when
$`\lambda_{3} \le N_{s} \le \lambda_{0}`$. For positive h the extremal triplet weight
maximizes $`\mu_{0} -\mu_{2}`$. With $`T =L -N_{s}`$, choose

```math
\begin{aligned}
\mu_0&=\min(\lambda_0,T-\lambda_2-\lambda_3),\\{}
\mu_2&=\max(\lambda_3,T-\mu_0-\lambda_1),\\{}
\mu_1&=T-\mu_0-\mu_2.
\end{aligned}
```

Assign $`(\mu_{0},\mu_{1},\mu_{2})`$ to `(t+,t0,t-)`, swapping t+ and t- for negative
h. This finite-dimensional weight optimization is checked against exhaustive
dominance enumeration through 16 rungs, and the resulting energies against
independent color-word Hamiltonians through six rungs. For L=6, N_s=4 and
small positive h, the minimizing population is `(4,2,0,0)`, not `(4,1,1,0)`:
simply shifting the old balanced representative loses a unit of Zeeman energy.

For a fixed M, there is no need to solve every individual triplet population.
Dominance of a four-color weight by lambda is equivalent to bounding all
single entries between $`\lambda_{3}`$ and $`\lambda_{0}`$, and all pair sums above
by $`P =\lambda_{0} +\lambda_{1}`$. At fixed N_s, with $`T =L -N_{s}`$, every triplet
population consequently lies in the same interval

```math
\mathrm{lower}=\max(\lambda_3,T-P),\qquad
\mathrm{upper}=\min(\lambda_0,P-N_s).
```

For nonnegative M write $`(N _+,N_{0},N _-)=(x +M,T -M -2x,x)`$ and intersect these
three bounds over integer x. Negative M swaps t+ and t-. An empty interval
means the multiplet is absent from this sector. Otherwise one representative
suffices: its permutation energy is independent of x. This constant-time
membership calculation is exhaustively checked against all population
vectors through 16 rungs. The shared scan still solves each relevant
highest-weight sea only once and reuses it across singlet-count sectors.

## Three nested real seas

Remove trailing zero rows and let k<=4 be the remaining number of colors.
There are k-1 root families, with
$`M_{0} =L`$, `M_a=sum_(b>a) mu_b`, and $`M_{k} =0`$. Define
$`\theta_{w} (x)=2 \arctan (x /w)`$. The equations solved are

```math
\begin{aligned}
L\theta_{1/2}(\lambda_j^1)-\sum_{\mathrm{same}}\theta_1(\lambda_j^1-\lambda_m^1)
+\sum_{\mathrm{next}}\theta_{1/2}(\lambda_j^1-\lambda_m^2)&=2\pi I_j^1,\\{}
\sum_{\mathrm{prev}}\theta_{1/2}(\lambda_j^a-\lambda_m^{a-1})
-\sum_{\mathrm{same}}\theta_1(\lambda_j^a-\lambda_m^a)
+\sum_{\mathrm{next}}\theta_{1/2}(\lambda_j^a-\lambda_m^{a+1})&=2\pi I_j^a.
\end{aligned}
```

The same-level self term is zero. The packed labels are
`I_j^a=j-(M_a-1)/2+shift_a`, with zero shift when $`M _{a -1}+M _{a +1}`$
is even and both shifts ±1/2 when it is odd. All independent displacement
choices are compared, modulo simultaneous reflection. Exact half-integer
labels use Uni20's `half_int`, not floating-point rounding. Energy and
momentum are

```math
E_{\mathrm{perm}}=L-\sum_j\frac1{(\lambda_j^1)^2+1/4},\qquad
\mathrm{momentum\_index}=\frac{LM_1-\sum_{a,j}2I_j^a}{2}\pmod L.
```

Full signed roots are solved together using an analytic Jacobian and
damped Newton corrections. A balanced thermodynamic density is only an
initial guess: neither that density nor a thermodynamic string
approximation replaces the finite-size equations.

The scan compares packed real-sea ground branches in all compatible
highest weights. It is not an enumeration of every eigenstate or a
general completeness proof for nested Bethe solutions. In particular,
`converged` means all branches required by this ground-state construction
were solved; it is not a certified spectral interval from diagonalization.

## Budgets, precision and library use

Both budgets apply to the **whole scan**, not separately to each singlet
sector: `--max-iterations` counts attempted Newton corrections and
`--max-branches` counts attempted independent sea branches (both default
to 10000). Each highest weight is solved once and reused in compatible
sectors. The number of four-row diagrams grows cubically with L; each
branch uses a dense Newton system with at most floor(3L/2) roots. This
first implementation favors moderate finite ladders, not thousand-rung
scans. Large calculations may need explicitly increased budgets.

An incomplete scan exits with status 2 and labels any retained energy a
**candidate upper bound**, not a sector minimum. Only converged branch
energies/roots are retained, never an unsolved seed. Zero Newton budget
can still find exact zero-root or zero-rapidity solutions; zero branch
budget gives no candidate, except for the direct product-state analytic
path. Upper-bound language is subject to ordinary numerical error, not
interval arithmetic. Invalid parameters exit with status 1.

The tolerance is the maximum logarithmic residual divided by L, default
`32*epsilon<Real>`; it is not an energy-error bound. All equations,
comparisons, root tables and energy shifts retain fp64, long-double or
fp128 precision. CPU time and scan counts appear in the overview.

```cpp
#include <bethe/ladder.hpp>

auto ground = bethe::ladder::ground_state<long double>(12, 1.0L);
auto sector = bethe::ladder::sector_ground_state<long double>(6, 4, 0.0L);
auto scan = bethe::ladder::sector_ground_states<long double>(12, -1.0L);
auto field_ground = bethe::ladder::ground_state<long double>(12, 5.0L, 3.0L);
auto field_sector = bethe::ladder::sector_ground_state<long double>(6, 4, 0.0L, -0.125L);
// Optional SolverOptions follows the field; zero-field overloads remain supported.
auto fixed_sz = bethe::ladder::magnetization_ground_state<long double>(12, 2, 1.0L);
auto sz_scan = bethe::ladder::magnetization_sector_ground_states<long double>(12, 2, 1.0L);
// Last optional argument fixes N_s too: (L, M, J_r, h, options, N_s).
auto both = bethe::ladder::magnetization_ground_state<long double>(6, 0, 0.0L, 0.0L, {}, 4);
// Check ground.converged / scan.complete before interpreting minima.
// State::energy is optional; State::highest_weight owns the representative roots.
```

Regression tests cover independent fixed-population permutation matrices
through seven rungs, literal spin-basis ladder matrices through four rungs,
the descendant counterexamples, original multiplicative Bethe equations,
momentum, analytic Jacobians, SU(2)/SU(3) reductions, both field signs,
magnetization/energy slopes, every fixed-Sz/singlet-count sector through six
rungs, literal spin-basis Sz blocks through four rungs, exact rung/one-triplet/polarized
limits and native-precision algebraic energies. Larger scans, budget
failures, parsing, formatting and generated literature citations are also
checked. See [the centralized bibliography](../CITATIONS.md#wang-1999).
