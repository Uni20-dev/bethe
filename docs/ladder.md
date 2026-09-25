# An integrable two-leg spin ladder

[Overview](../README.md) · [SU(3) chain](su3.md) · [Model catalogue](models.md)

`bethe-ladder-pbc` implements Wang's periodic ladder, with **a particular
four-spin interaction required for integrability**. It is not a solver for
an ordinary two-leg Heisenberg ladder at arbitrary couplings. There are L
rungs, two spin-1/2 operators S and T per rung, and longitudinal field h:

```text
H = sum_j [S_j.S_(j+1) + T_j.T_(j+1)
           + 4 (S_j.S_(j+1))(T_j.T_(j+1))] + J_r sum_j S_j.T_j
    - h sum_j (Sz_j+Tz_j).
```

The leg coefficient is 1, the four-spin coefficient is 4, and `J_r` can
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
`--sz M` additionally fixes total physical Sz to the integer M in `[-L,L]`;
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
even for odd L. There is no extra parity constraint on `L-N_s-M`, because
the triplet t0 has zero projection. Feasibility is `|M|<=L-N_s`.
With `--sectors`, only `N_s=0,...,L-|M|` are listed. Incompatible explicit
constraints are input errors, not empty successful calculations.
At `|M|=L` the unique polarized state is returned analytically, even when
it is *not* the unconstrained ground state and even with zero solve budgets.

These results are minima of **Sz sectors**, not states of prescribed total
spin S, and do not enumerate excited levels within a sector. At h=0 a
spin-S multiplet contributes to all `|M|<=S`; equal sector energies can
therefore be ordinary SU(2) degeneracies.

## Why a ladder becomes a four-color chain

Each rung has one singlet s and three triplets t+, t0, t-. Their rung
energies are `-3J_r/4`, `J_r/4-h`, `J_r/4`, `J_r/4+h`, respectively. The identity

```text
P_rung(j,k) = (2 S_j.S_k + 1/2)(2 T_j.T_k + 1/2)
```

turns the leg/four-spin part into `P_rung-1/4`. Thus

```text
H = sum_j P_rung(j,j+1) - L/4 + J_r*(L/4-N_s) - h*M,
E = E_perm - L/4 + J_r*(L/4-N_s) - h*M.
```

The singlet count N_s is conserved. Only a chemical-potential term changes
as J_r and h vary; a particular multiplet branch has roots independent
of both, although the minimizing multiplet can change. The field is a linear
combination of conserved color populations, so the same Bethe equations apply.
We use the permutation form in
[Wang's Eqs. (2)–(5)](https://arxiv.org/html/cond-mat/9901168), with
`J_r=2J` in those equations and the energy constants restored as above.

Two useful exact checks follow. The all-singlet product has
`E=3L(1-J_r)/4`. It is a global ground state for `J_r-|h|>=4`, so the default
calculation returns it analytically, even with zero iteration/branch
budgets. This sufficient condition follows from the operator inequality
`sum(P-1)>=-4*N_triplet`; it does not assert uniqueness on the boundary or locate
every finite-ring crossing. The lowest one-triplet energy above the product is
`J_r-|h|-4*sin²(pi*floor(L/2)/L)`. The fully polarized triplet product is
also analytic when `|h|>=4+max(J_r,0)`, with
`E=L*(3/4+J_r/4-|h|)` and `M=sign(h)*L`. The same bound applies to defects
of the favored triplet color: every other rung color then costs at least 4.
These are sufficient conditions, not all finite-ring phase boundaries.
In the N_s=0, h=0 sector the model reduces to an SU(3) triplet permutation
chain plus the constant `L*(J_r-1)/4`.

## Finite rings: highest weights are not physical populations

An SU(4) highest weight is a Young diagram with row lengths
`mu_1>=mu_2>=mu_3>=mu_4>=0`, summing to L. A multiplet contains a population
vector w exactly when mu dominates its sorted entries: every prefix sum
of mu is at least the corresponding prefix sum of sorted(w).

This matters on a ring. It is tempting to solve just the packed real sea
for the physical populations and call it the sector minimum. That is
wrong even for six rungs:

| Physical populations (s,t+,t0,t-) | Own highest-weight sea E_perm | Actual minimizing highest weight | Minimum E_perm |
| --- | --- | --- | --- |
| 4,1,1,0 | -1 | 4,2,0,0 | `1-sqrt(5)` |
| 2,2,1,1 | `-(5+sqrt(17))/2` | 2,2,2,0 | `-1-sqrt(13)` |

In both cases the lower state is an SU(4) **descendant**. Lowering operators
change its populations without changing the permutation energy or momentum;
the physical rung contribution still uses the requested N_s. The report
therefore keeps physical populations, highest-weight rows, and the
descendant flag separate. `--roots` prints the finite roots of that
highest-weight representative, **not fabricated finite roots for its
descendant**. With `--sectors --roots`, only the selected best state's
roots are printed.

At h=0, for N_t=L-N_s triplets, it suffices to use the population representative
`(N_s,q+[r>0],q+[r>1],q)`, where `N_t=3q+r`. Every SU(3) triplet multiplet
contains this balanced weight. Scanning all compatible SU(4) highest
weights therefore retains every triplet-multiplet choice at fixed N_s.
This weight-space statement should not be confused with assuming a
generic ordering of periodic multiplet energies:
[Hakobyan's periodic ordering theorem](https://arxiv.org/abs/cond-mat/0403587)
has an additional row-parity condition.

At nonzero h we instead minimize the Zeeman term within each multiplet.
Writing its rows as `lambda_0,...,lambda_3`, restriction to the three triplet
colors gives interlacing diagrams `lambda_i >= mu_i >= lambda_(i+1)` with
`sum(mu)=L-N_s`. A singlet count occurs precisely when
`lambda_3 <= N_s <= lambda_0`. For positive h the extremal triplet weight
maximizes `mu_0-mu_2`. With `T=L-N_s`, choose

```text
mu_0 = min(lambda_0, T-lambda_2-lambda_3)
mu_2 = max(lambda_3, T-mu_0-lambda_1)
mu_1 = T-mu_0-mu_2.
```

Assign `(mu_0,mu_1,mu_2)` to `(t+,t0,t-)`, swapping t+ and t- for negative
h. This finite-dimensional weight optimization is checked against exhaustive
dominance enumeration through 16 rungs, and the resulting energies against
independent color-word Hamiltonians through six rungs. For L=6, N_s=4 and
small positive h, the minimizing population is `(4,2,0,0)`, not `(4,1,1,0)`:
simply shifting the old balanced representative loses a unit of Zeeman energy.

For a fixed M, there is no need to solve every individual triplet population.
Dominance of a four-color weight by lambda is equivalent to bounding all
single entries between `lambda_3` and `lambda_0`, and all pair sums above
by `P=lambda_0+lambda_1`. At fixed N_s, with `T=L-N_s`, every triplet
population consequently lies in the same interval

```text
lower = max(lambda_3, T-P)
upper = min(lambda_0, P-N_s).
```

For nonnegative M write `(N_+,N_0,N_-)=(x+M,T-M-2x,x)` and intersect these
three bounds over integer x. Negative M swaps t+ and t-. An empty interval
means the multiplet is absent from this sector. Otherwise one representative
suffices: its permutation energy is independent of x. This constant-time
membership calculation is exhaustively checked against all population
vectors through 16 rungs. The shared scan still solves each relevant
highest-weight sea only once and reuses it across singlet-count sectors.

## Three nested real seas

Remove trailing zero rows and let k<=4 be the remaining number of colors.
There are k-1 root families, with
`M_0=L`, `M_a=sum_(b>a) mu_b`, and `M_k=0`. Define
`theta_w(x)=2 atan(x/w)`. The equations solved are

```text
L theta_(1/2)(lambda_j^1)
  - sum_same theta_1(lambda_j^1-lambda_m^1)
  + sum_next theta_(1/2)(lambda_j^1-lambda_m^2) = 2*pi*I_j^1,

sum_prev theta_(1/2)(lambda_j^a-lambda_m^(a-1))
  - sum_same theta_1(lambda_j^a-lambda_m^a)
  + sum_next theta_(1/2)(lambda_j^a-lambda_m^(a+1)) = 2*pi*I_j^a.
```

The same-level self term is zero. The packed labels are
`I_j^a=j-(M_a-1)/2+shift_a`, with zero shift when `M_(a-1)+M_(a+1)`
is even and both shifts ±1/2 when it is odd. All independent displacement
choices are compared, modulo simultaneous reflection. Exact half-integer
labels use Uni20's `half_int`, not floating-point rounding. Energy and
momentum are

```text
E_perm = L - sum_j 1/[(lambda_j^1)^2+1/4],
momentum_index = (L*M_1 - sum_(a,j) 2*I_j^a)/2 mod L.
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
