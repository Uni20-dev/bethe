# Two bound pairs scattering on a free-end chain

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [One pair plus real roots](biquadratic-pair-defect.md) · [Isolated bound pairs](biquadratic-bound-pairs.md)

Four TL singlet insertions can form two bound pairs. This selected-state
library solver follows that family in **ell=N-8**, for odd/even N>=8.
It retains both finite string deviations in fp64, native long double or
enabled fp128. CLI scans and table exports for this family are a later
integration step; `--pair-defects` still means **one** pair plus real roots.

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;

auto state = ferro::two_bound_pairs<long double>(128, {121, 122});
if (state.reference.converged) {
    auto gap = state.tl_energy; // E-E0, evaluated directly; E0=N-1
    auto energy = state.energy;
    auto centers = state.reference.centers;
    auto logarithmic_deviations = state.reference.log_deviations;
}
```

The arguments are **ordered Bethe labels** `1<=J1<J2<=N-6`, not energy
ranks or physical lattice momenta. The high-label corner `(N-7,N-6)`
approaches `E-E0=10/3`, the separated-pair threshold `5/3+5/3`.
It does not give the minimum of the entire four-defect module; other
string topologies are missing. These are TL insertions, not spin flips,
and no operator spectral weights or physical SU(2) labels are supplied.

## Coupling two existing pair systems

The reference remains quantum-group XXZ with opposite end fields,
Delta=3/2. Its roots are

```text
u_i,± = (eta+d_i)/2 ± i*a_i/2,      i=1,2,
0<a_1<a_2<pi,                     eta=acosh(Delta),
d_i = sigma_i*exp(-L_i),           |d_i|<eta,
sigma_i = (-1)^(N-J_i-i).          (i is one-based here)
```

The ordering-dependent sign is essential. Applying the isolated-pair
rule to both strings loses a pi phase in the second pair's equation.
`L_i` and `sigma_i` are authoritative even if exp(-L_i) underflows or the
deviation disappears when added to eta. Distinct ordered centers prevent
coincident pair roots; a seed that cannot resolve them at the selected
precision is rejected.

The solver reuses each pair's phase, log-modulus, analytic Jacobian and
direct energy contribution from the [two-string system](xxz-open-two-string.md).
Only the pair-pair coupling is new. The shared damped Newton driver solves
four real unknowns, with work and storage per iteration independent of N.
The ideal-string stage initializes the solve; convergence always uses the
finite-deviation equations and the original iteration/tolerance budget.

Here is our regularization of
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020). For pair i scattering
with j, set `t=(d_i-d_j)/2`, `s=(d_i+d_j)/2` and `beta=a_i±a_j`. Define

```text
Theta(beta;w) = 2 atan2(sin(beta/2), tanh(w)*cos(beta/2)),
C(beta;w) = 2 atan(tanh(w)/tan(beta/2)),
G(beta;w) = log|sinh(w+i*beta/2)|,

F(beta) = Theta(beta;eta+t) + Theta(beta;eta-t)
          + Theta(beta;2eta+s) + C(beta;s),
K(beta) = G(beta;eta+t) - G(beta;eta-t)
          + G(beta;2eta+s) - G(beta;s).
```

From each isolated pair's scaled phase residual subtract
`[F(a_i-a_j)+F(a_i+a_j)]/(2N)`. From its log-modulus residual subtract
the analogous sum of K. The complementary phase C is chosen continuously
for negative as well as positive beta; using the unadjusted `atan2`
complement would introduce a spurious 2pi jump at negative beta.

The original unsquared pair equation includes a `3pi` branch contribution
when `a_i>a_j`. This gives the deviation signs above. The upper center's
excluded a=pi endpoint has label N-5, leaving the largest admissible
candidate label N-6. The energies add through their direct auxiliary
shifts, then use the common ferro mapping `gap=-2*energy_shift`.

The lower API is
`bethe::xxz::quantum_group::two_pairs::solve(N,Delta,{J1,J2},options)`.
It accepts finite Delta>1, but the string ansatz is not a guarantee of
convergence or branch existence for every parameter. Failed iterates retain
consistent estimates, signs and residuals, not verified energies.

## Validation and remaining four-defect states

All `binomial(N-6,2)` two-pair candidates match distinct module ED levels
for N=8,...,10 at Delta=1.25, 1.5, 2 and 3. At the biquadratic value, combining
them with the existing four-real-root and one-pair-plus-two-real-root
families also matches ED **as a multiset**, so these are distinct subsets:

| N | Four real roots | One pair + two real roots | Two pairs | Module dimension | Still missing |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1 | 6 | 1 | 14 | 6 |
| 9 | 5 | 20 | 3 | 42 | 14 |
| 10 | 15 | 45 | 6 | 90 | 24 |

This is finite-size numerical evidence, **not** a general completeness
proof. A three-string plus a real root and a four-string droplet remain
separate targets. Tests additionally check the original complex equations
for both deviation parities, every analytic Jacobian column in ideal and
finite-deviation modes, failure diagnostics, and the 10/3 threshold through
N=100000, where both deviations underflow but their logarithms remain finite.
