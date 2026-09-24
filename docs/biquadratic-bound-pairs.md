# Two-defect bound pairs on long open chains

[Ferromagnetic overview](biquadratic-ferromagnetic.md) · [Signed-string equations](#equations-and-branch-labels) · [Exports](output.md)

For the spin-1 chain `H=+sum(S.S)^2`, two singlet insertions can bind.
The targeted solver follows the **one-two-string family in ell=N-4**, without
enumerating scattering states or searching for every Q polynomial.
It has only two unknowns, independent of chain length, and supports odd/even
N>=4 in fp64, native long double and enabled fp128.

```sh
# First eight modes of the bound-pair family, not the first eight levels
# of the entire physical spectrum:
build/bethe-biquadratic-obc 128 --ferromagnetic --bound-pairs 8
# Every supported two-string label at this size:
build/bethe-biquadratic-obc 65 --ferromagnetic --bound-pairs all --roots
# Long chain, native higher precision, with reusable exports:
build/bethe-biquadratic-obc 100000 --ferromagnetic --bound-pairs 4 \
  --precision long-double --json pairs.json --csv pairs.csv
```

`COUNT` selects modes 1,...,min(COUNT,N-3); `all` means all N-3 candidate
two-string labels, **not all states of the two-defect module**. A row budget
`--max-candidates` (default 10000) is checked before allocating results.
Selecting four modes on a large chain needs only four solves, not a scan of
all modes. `--bound-pairs` requires `--ferromagnetic` and cannot be combined
with other state selectors; it fixes the through-line sector automatically.

## Excitation structure

The gap is measured from the entire exact ground manifold, `E0=N-1`.
For N=4 there is one bound-pair mode:

```text
E = (15-sqrt(17))/2,
E-E0 = (9-sqrt(17))/2 = 2.438447187...
```

This is the complex-root singlet missed by the all-real scan. At larger N,
ell=N-4 need not be a physical singlet: each TL eigenvector carries the
representation multiplicity m_(N-4), with several possible physical spins.

The lowest two-defect branch approaches **5/3**, below the energy 2 needed
for two widely separated lowest one-defect excitations (each has limiting
gap 1). This is a bound-pair threshold, not the global positive gap, which
remains 1. The droplet theorem of
[Nachtergaele, Spitzer and Starr](../CITATIONS.md#nachtergaele-spitzer-starr-2007),
Theorem 2.1, gives the thermodynamic module-edge limit. Rescaling their
Hamiltonian by `2*Delta=3`, with eta=acosh(3/2), gives
`2*sinh(eta)*tanh(eta)=5/3` for two defects.

The lowest mode matches the module minimum in independent ED checks for
N=4,...,10. We do not use that finite validation as a proof of finite-N
ordering at arbitrary N. Other modules and scattering levels can interleave
the modes shown. The mode number and string center are **not physical
lattice momentum**; free ends have no translation quantum number.
No form factors or coupling strengths to a chosen ground state are supplied.

## Equations and branch labels

Use the auxiliary XXZ chain with Delta=3/2 and opposite end fields, not the
ordinary zero-field open XXZ chain. Our regularization of
[Bajnok et al., Eq. (5.12)](../CITATIONS.md#bajnok-2020) writes the two roots as

```text
u_± = (eta+d)/2 ± i*a/2,        0<a<pi, |d|<eta,
d = sigma*exp(-L),              L=-log|d|,
J = N-2-mode,                   mode=1,...,N-3,
sigma = (-1)^(mode+1).
```

There is **no real-root sea**. Using Theta, C and G from the
[two-string derivation](xxz-open-two-string.md), the two equations are

```text
2N [Theta(a;eta+d/2)+C(a;d)] - 2 Theta(2a;eta) = 2 pi J,

2N [G(a;eta+d/2)-G(a;d/2)] - log(sinh(2eta+d)) - L
 + log(sinh(d)/d) = 0.
```

Both signs of d are essential. The reflected scattering factor has a
negative denominator when d<0, changing its phase by pi. Consequently
positive d requires N-J odd, negative d requires N-J even. The lowest-edge
mode has positive deviation and the next has negative deviation. Keeping
only positive deviations would omit every second mode.

The limit `log(sinh(d)/d)=d^2/6+...` is regular for either sign. Retaining L
avoids losing the finite-size deviation when adding it to eta rounds it
away, or when exp(-L) underflows. An ideal-string stage supplies only the
initializer; convergence always checks the finite-deviation equations.
Both stages share the same iteration budget and selected precision.
The two-variable Newton correction is equilibrated for the very different
angle and logarithm scales; this does not rescale the convergence tolerance.
Work and storage per Newton step are O(1); output storage is O(COUNT).

## Diagnostics and library API

Output tables are `states`, `reference`, `string`, and optionally `roots`.
`string` records the sign and L separately; `deviation` is null on underflow,
not an exact zero. Rounded root coordinates must not be used to reconstruct
the deviation. In `states`, `tl_energy=E-E0` is evaluated directly from the
pair contribution, without subtracting extensive total energies.

Unconverged modes remain as estimates with `converged=false`, no published
`gap`, and exit status 2. The exact ground reference remains valid. Residuals
are phase/log-modulus errors divided by 2N, not energy-error bounds. Tight
tolerances do not certify more digits than the chosen arithmetic resolves.
Very large N can make distinct root coordinates unresolvable; such cases
must not be interpreted as a certified set of distinct modes.
An unresolvable initial center is rejected explicitly.

```cpp
#include <bethe/biquadratic_ferromagnetic.hpp>
auto s = bethe::biquadratic::ferromagnetic::bound_pair<long double>(128, 1);
if (s.reference.converged) {
  auto gap = s.tl_energy;
  auto energy = s.energy;
}
```

The lower API is `xxz::quantum_group::two_string::bound_pair(N,Delta,mode)`.
It returns auxiliary-sign energies, not physical ferro energies. Delta>1
is accepted, but away from 3/2 an individual finite-size branch can collapse
or fail to converge; no general-anisotropy completeness claim is made.

Validation covers independent ED on N=4,...,10, both deviation signs in the
original complex equations, analytic Jacobians and native-precision exact
values. Long-chain tests reach N=100000, including deviation underflow.
The next extensions are multi-defect droplets and scattering branches,
physical-spin resolution, and operator-dependent spectral weights.
