# A bound pair scattering with one defect

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Real-root scattering](biquadratic-scattering.md) · [Bound triples](biquadratic-bound-triples.md)

Three TL singlet insertions need not all bind, or all remain separate. A
two-string plus one real root describes the intermediate **pair-plus-defect**
branch in ell=N-6. A selected-state library API now supports this branch on
odd/even N>=6, in fp64, native long double and enabled fp128. CLI scans and
exports for this family are the next integration step; `--excitations` still
means the purely real-root family.

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
namespace ferro = bethe::biquadratic::ferromagnetic;

std::size_t const n = 128;
auto state = ferro::pair_defect<long double>(n, n - 3, n - 5);
if (state.reference.converged) {
    auto gap = state.tl_energy; // E-E0, E0=N-1; evaluated directly
    auto energy = state.energy;
    auto alpha = state.reference.rapidities[0];
    auto a = state.reference.center;
    auto L = state.reference.log_deviation;
}
```

The two arguments after N are integer Bethe labels **I,J**, not energy
ranks or physical lattice momenta:

```text
I = 1,...,N-3       real root
J = 1,...,N-5       bound pair
```

The high-label corner `(N-3,N-5)` approaches the separated-cluster threshold
`5/3 + 1 = 8/3`. This is above the three-defect droplet edge 2 and below the
three-unbound-defect edge 3. It is not the global first excitation, and the
labels do not imply a universal ordering of finite-size levels. All of these
counts refer to TL insertions, **not physical spin flips**. Operator spectral
weights and physical-spin decomposition remain separate questions.

## Same equations, different labels

The auxiliary model is quantum-group XXZ with opposite end fields at
Delta=3/2, not the usual zero-field open XXZ chain. We reuse the
[two-string equations](xxz-open-two-string.md) with one real root:

```text
v = i*alpha/2,
u_± = (eta+d)/2 ± i*a/2,
eta = acosh(Delta),  d = sign*exp(-L),
0 < alpha,a < pi,   |d| < eta.
```

Replace the real-root right-hand side by `2*pi*I`, and the pair-phase
right-hand side by `2*pi*J`. The modulus equation is unchanged. Retaining
the unsquared pair equation fixes `sign(d)=(-1)^(N-J-1)`; the product-phase
equation alone loses this sign information. The real root may cross the
pair center: no additional ordering condition `alpha<a` is imposed.
In the continuous logarithmic convention, the alpha=pi endpoint has
I=N-2 and the a=pi endpoint has J=N-4; those excluded endpoint roots give
the label bounds above.

These are our regularized equations derived from
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020), not an ideal-string
energy formula. L and the separate deviation sign remain authoritative
when exp(-L) underflows. The existing analytic Jacobian, damped Newton
driver, iteration budget, native precision and direct energy-shift mapping
are shared with the AF singlet and isolated bound-pair solvers. There are
three unknowns, so storage and work per Newton step do not grow with N.

`reference.real_labels` stores I and `reference.string_label` stores J.
The common cluster result's `mode` is zero for this mixed branch: there is
no one-dimensional mode ordering. An unconverged result retains estimates
and residuals, **not verified energies**. No tolerance relaxation, hidden
extra iterations or change of precision is performed.

The lower API
`bethe::xxz::quantum_group::two_string::pair_defect(N,Delta,I,J,options)`
accepts finite Delta>1. Branch collapse or convergence failure can occur;
the biquadratic Delta=3/2 coverage below is not a guarantee for other Delta.

## Checks and limits of coverage

For N=6,...,10, independent exact diagonalization verifies every level of
the ell=N-6 module when the three families are combined:

| Family | Candidate count |
| --- | ---: |
| Three positive real roots | binomial(N-3,3) |
| One pair plus one real root | (N-3)(N-5) |
| One three-string | N-5 |

Their sum is the module dimension `binomial(N,3)-binomial(N,2)`. Matching
the dimension alone would not establish completeness: the regression
compares the **entire sorted spectrum including multiplicities**. This is
small-chain numerical evidence, not a completeness proof for arbitrary N.

Further tests check the original complex equations for both deviation signs
and odd/even chains, the analytic Jacobian, native-precision six-site
characteristic polynomial, zero/one-iteration failure diagnostics, and the
8/3 threshold up to N=100000. The generic finite-size string description
still permits numerical failure and does not provide physical SU(2)
labels, scattering amplitudes or form factors.
