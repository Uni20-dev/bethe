# Three-defect bound droplets on open chains

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Bound pairs](biquadratic-bound-pairs.md) · [Exports](output.md)

Three TL singlet insertions can bind into a droplet, just as two can form a
bound pair. `--bound-triples` follows the **one-three-string family in
ell=N-6** for `H=+sum(S.S)^2`, on odd/even N>=6. It solves three real
unknowns, independent of chain length, in fp64, native long double or enabled
fp128. This is a finite-size Bethe calculation, not an ideal-string energy
substituted for a finite chain.

```sh
build/bethe-biquadratic-obc 128 --ferromagnetic --bound-triples 8
build/bethe-biquadratic-obc 65 --ferromagnetic --bound-triples all --roots
build/bethe-biquadratic-obc 100000 --ferromagnetic --bound-triples 4 \
  --precision long-double --json triples.json --csv triples.csv
```

`COUNT` selects modes 1,...,min(COUNT,N-5); `all` selects all N-5 candidate
three-string labels, **not all states of the three-defect module**. The
`--max-candidates` row budget (default 10000) is checked before allocation.
This option requires `--ferromagnetic`, fixes ell=N-6, and cannot be mixed
with other state selectors, including `--bound-pairs`.

## Where this branch lies

The gap is measured from the entire exact ground space, `E0=N-1`. Its
lowest mode approaches **E-E0=2**, compared with 3 for three widely
separated lowest one-defect excitations, or 8/3 for a separated bound pair
and one defect. Two separated single defects also have threshold 2, but
belong to a different TL module. These are not physical spin counts.

The thermodynamic droplet limit follows from
[Nachtergaele, Spitzer and Starr, Theorem 2.1](../CITATIONS.md#nachtergaele-spitzer-starr-2007).
In our normalization, at eta=acosh(3/2),

```text
g_M = 2*sinh(eta)*tanh(M*eta/2),
g_1 = 1,    g_2 = 5/3,    g_3 = 2,    g_infinity = sqrt(5).
```

That theorem supplies a thermodynamic module-edge limit, not a proof of
our finite-N mode ordering. The lowest mode matches the module minimum in
independent ED for N=6,...,10. Other modules and scattering levels can
interleave these modes. The global positive gap still approaches 1.

At N=6 the three-defect module has dimension five. Its TL characteristic
polynomial is `(g-6)(g-7)(g^3-17g^2+80g-106)`. The targeted droplet is the
smallest cubic root, `g=2.28668884804195...`, with physical energy `5+g`.
This makes a useful native-precision check independent of the Bethe equations.
At N=128 the first gap is approximately `2.00008198366582`.

## Why the deviation needs a phase

The reference is quantum-group XXZ with opposite end fields, Delta=3/2,
not ordinary zero-field open XXZ. A canonical representative of its three
roots is

```text
u0 = i*a/2,
u+ = eta + Re(z) + i*(a/2 + Im(z)),    u- = conjugate(u+),
z = exp(-L + i*phi),                  eta = acosh(Delta),
J = N-4-mode,                        mode=1,...,N-5.
```

The central root and conjugate pair together make the three-string; there
is no additional real-root sea. Unlike a two-string's signed real deviation,
`z` generally has both real and imaginary parts. Setting its phase to zero
does not solve the finite-chain equations. `L=-log|z|` and `phi` remain
authoritative even when adding z to eta rounds away the correction, or
exp(-L) underflows.

Here is our regularization of the original complex equations
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020). Put `u=u+`,
`b=a+2 Im(z)` and define

```text
Theta(beta;w) = 2 atan2(sin(beta/2), tanh(w)*cos(beta/2)),
D = log sinh(u+eta/2) - log sinh(u-eta/2),
A = log sinh(2eta+z),
B = log sinh(2eta+z+i*a) - log sinh(z+i*a),
C = log sinh(eta+i*b) - log sinh(-eta+i*b),
R = log sinh(3eta+2 Re(z)) - log sinh(eta+2 Re(z)),
c = log(sinh(z)/z),                   c(0)=0.
```

With `wrap(t)=atan2(sin(t),cos(t))`, the three residuals are

```text
f0 = Theta(a;eta/2) + 2 Im(D)
     - [2 Im(B) + Theta(2b;eta) + pi*(J+1)]/N,
f1 = Re(D) - [Re(A) + L - Re(c) + Re(B) + R]/(2N),
f2 = wrap(2N Im(D) - Im(A) + phi + Im(c) - Im(B) - Im(C))/(2N).
```

The first equation is the product phase with the singular internal
scattering cancelled; the other two retain the original outer-root equation.
For tiny z, the regular correction c is evaluated by its Taylor series.
The search stays in `0<a<pi`, `0<b<pi`, `|z|<eta/4`.

An ideal-string stage initializes the solve. Its center satisfies
`N*Theta(a;3eta/2)-Theta(2a;eta)-Theta(2a;2eta)=pi*J`.
Final convergence always uses the finite-deviation residuals above; both
stages share one iteration budget. Column/row equilibration improves Newton
conditioning without relaxing the tolerance. Each Newton step costs O(1).

## Output, validation and reuse

Tables are `states`, `reference`, `string`, and optional `roots`. The string
table records a, L and phi, with `deviation_real`/`deviation_imag` null on
magnitude underflow. Root coordinates are rounded, not a substitute for
these logarithmic parameters. The center and mode are **not lattice momentum**.

`tl_energy` evaluates the gap contribution directly, avoiding subtraction
of extensive energies. Multiplicity is m_(N-6), the number of physical
states per TL eigenvector, not a physical SU(2) spin label; it is null on
uint64 overflow. Unconverged modes retain estimates but have no published
`gap` and give exit status 2. The exact ground reference remains valid.
Residuals are normalized equation errors, not certified energy-error bounds.
Unresolvable seeds are rejected; an arbitrarily small tolerance cannot add
digits beyond the selected arithmetic.

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
auto s = bethe::biquadratic::ferromagnetic::bound_triple<long double>(128, 1);
if (s.reference.converged) {
  auto gap = s.tl_energy;
}
```

The equations live in `xxz::quantum_group::three_string::bound_triple(N,Delta,mode)`.
This lower API returns auxiliary-sign XXZ energies; the biquadratic adapter
only applies the physical energy convention and representation multiplicity.
General finite Delta>1 is accepted, but individual branches can collapse or
fail away from 3/2; no general-anisotropy completeness claim is made.

Both string families and the AF two-string singlet use the common
[logarithmic-string Newton driver](../include/bethe/detail/open_string_solver.hpp),
including phase functions, ideal-stage handling, line search, scaling and
budget/status rules. Pair/triple CLI scans, state/reference tables, exports
and their output regression checks are also shared.

Tests cover original complex equations (including the central root), analytic
Jacobians, native-precision N=6 energies, all branch labels at small odd/even
sizes, and ED through N=10. XXZ ED checks also exercise Delta=1.25, 2 and 3.
Long-chain checks reach N=100000, including both family edges and deviation
underflow. Scattering branches, larger droplets, physical-spin resolution
and operator-dependent spectral weights remain separate extensions.
