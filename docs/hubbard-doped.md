# Doped Hubbard spinon and charge lines

Use the same `bethe-hubbard-dispersion` frontend with `--density n` for
`0<n<1`, `U>0`, zero magnetic field, and hopping `t=1`. This is a
thermodynamic calculation, not a sequence of finite-ring shell fillings.
`fp64`, native `long-double`, and optional `fp128` are supported throughout.
Omitting `--density` retains the [half-filled calculation](hubbard-dispersion.md).

```sh
build/bethe-hubbard-dispersion --u 4 --density 0.75 --reference fermi --format csv
build/bethe-hubbard-dispersion --u 4 --density 0.5 --branch spinon --momentum 1
build/bethe-hubbard-dispersion --u 4 --density 0.75 --precision fp128 --points 101 --format tsv
```

## Which lines are included?

| Branch | DeltaN | Spin | Unwrapped momentum interval |
|---|---:|---:|---|
| `spinon` | 0 | 1/2 | `[0,pi*n]` |
| `holon` | -1 | 0 | `[-pi*n/2,3*pi*n/2]` |
| `charge-particle` | +1 | 0 | `[pi*n/2,2*pi-3*pi*n/2]` |

`--branch all` selects these three at `n<1`. The charge particle adds a real
charge root outside the occupied sea. It is **not** the gapped antiholon of
the half-filled Mott insulator; that branch is available only at `n=1`.
This distinction follows the real particle/hole construction in
[Luo–Pu–Guan, Sec. II](https://arxiv.org/html/2307.00890v3).
Gapped string excitations, densities above one, attraction, finite fields,
continuum thresholds and spectral weights are not implemented.

Grids include both endpoints and are uniform in dressed momentum, not rapidity.
`--momentum p` replaces the grid. With `--branch all`, it must lie in all
three intervals. Momenta are in radians for a one-site unit cell; fold modulo
`2*pi` or into an enlarged iMPS Brillouin zone as appropriate.

These fractional lines are not isolated finite-ring electron levels. Momentum
origins depend on the choice of Fermi-point spectators/domain-wall gauge;
the conventions below fix ours explicitly. Match this gauge, spin and charge
before comparing an iMPS line. No automatic spectral-threshold identification
is made.

## Hamiltonian energies versus Fermi-referenced energies

There are two independent choices:

- `--convention symmetric|unshifted` selects the interaction in the Hamiltonian,
  exactly as in the [half-filled guide](hubbard-dispersion.md#hamiltonian-and-energy-zero).
- `--reference hamiltonian|fermi` selects `DeltaE_H` (the default) or
  `DeltaE_H-mu_H*DeltaN`.

The code solves the background chemical potential and reports it in both
conventions:

```text
mu_symmetric = mu_unshifted - U/2
E_fermi = E_unshifted - mu_unshifted*DeltaN
        = E_symmetric - mu_symmetric*DeltaN.
```

All three doped lines have zero endpoint energy **in the Fermi reference**.
The Hamiltonian charge energies need not be positive. For example, at
`U=4, n=0.5`, `mu_unshifted` is approximately `-0.726930427436425` and
`mu_symmetric` is approximately `-2.726930427436425`. A charge particle right
at the Fermi point has that chemical-potential cost under the corresponding
Hamiltonian, while its Fermi-referenced energy is zero. Spinons are unchanged.

When comparing to an iMPS calculation that minimizes `H-mu*N`, choose `fermi`
and check the reported `mu` against the simulation. `symmetric_energy` in the
table remains the symmetric **Hamiltonian** difference even in this mode;
`fermi_energy` is always the subtracted energy.

There is a useful subtlety as `n` approaches one from below. The chemical
potential approaches the **lower edge** of the Mott plateau, whereas the
exactly half-filled tool chooses its **middle**. Hamiltonian hole and spinon
energies have the expected continuous limit; the two Fermi energy zeros do
not coincide. The real charge-particle momentum interval shrinks to a point,
not to the half-filled antiholon band.

## Equations and numerical method

We eliminate the infinite zero-field spin sea analytically and solve the
finite charge-sea Fredholm equations of
[Essler, Eqs. (103)–(106)](https://arxiv.org/html/1002.1671).
Let `u=U/4`, with integrals below over `[-Q,Q]`:

```text
R(x) = (1/pi) integral_0^infinity cos(w*x)/(1+exp(2*u*w)) dw
s(x) = sech(pi*x/(2*u))/(4*u)

rho(k) = 1/(2*pi) + cos(k) integral R(sin(k)-sin(k')) rho(k') dk'
eps_c(k) = -2*cos(k)-mu_unshifted
           + integral cos(k') R(sin(k)-sin(k')) eps_c(k') dk'
eps_s(lambda) = integral cos(k) s(lambda-sin(k)) eps_c(k) dk.
```

Fix `Q` by `integral rho=n`, and the chemical potential by `eps_c(+-Q)=0`.
Our `mu_unshifted` is `mu+2u` in Essler's Eq. (105). The Fermi energies are
`-eps_s`, `-eps_c` inside the sea, and `+eps_c` outside it.

To remove any momentum-origin ambiguity, our definitions are

```text
A(x) = integral_0^x R(y) dy
p_c(k) = k + 2*pi integral rho(k') A(sin(k)-sin(k')) dk'
p_s(lambda) = 2 integral rho(k) atan(exp(pi*(sin(k)-lambda)/(2*u))) dk
p_holon(k) = pi*n/2 - p_c(k)
p_particle(k) = p_c(k) - pi*n/2  (add 2*pi on the negative-k segment).
```

Here `p_c(+-Q)=+-pi*n` and `p_c(+-pi)=+-pi`. The particle grid passes
continuously through the zone boundary, even though its bare parameter jumps
from `pi` to `-pi`. Spinons have infinite rapidity at the two endpoints.

Evenness reduces the linear systems to `[0,Q]`. Gauss–Legendre nodes,
kernel evaluations, pivoted solves and momentum inversion all use the selected
scalar. `R` uses a convergent Euler series; `A` uses a paired log-gamma-ratio
expansion with rational coefficients. No fp64 kernel table or external Python
runtime is used. The background is solved once and reused across every point.

## Accuracy, budgets and failure

The default tolerance is `4096*epsilon` of the selected scalar. We double the
positive-half quadrature order from 16 up to 256 and require successive
backgrounds to agree in `Q`, chemical potential, ground energy, and sampled
dispersions. Each requested point is then inverted on **both** final meshes
at the same physical momentum. These are numerical error estimates, not
rigorous interval bounds.

- `energy_mesh_error` estimates the selected energy's mesh difference and
  roundoff, including the chemical-potential shift when appropriate. It
  excludes propagation of momentum-inversion error.
- `momentum_error` reports inversion residuals, mesh disagreement at the
  returned parameter, and the density residual contribution, in radians.
- Metadata includes the background status, `Q`, both chemical potentials,
  unshifted ground energy per site, mesh estimate, node count, work and CPU time.

`--tolerance` sets the dimensionless target in units `t=1`. Background mesh
agreement is absolute; point energies use `tol*max(1,abs(E))`. Charge momentum
inversion is absolute in radians; spinon inversion scales with the distance
to the nearer endpoint. Defaults allow 64 density solves **across all meshes**
and 160 momentum updates **across both meshes per point**.
Use `--initial-nodes`, `--max-nodes` (at most 512),
`--max-background-iterations` and `--max-iterations` to change these budgets.
Half-filled-only `--max-evaluations`/`--max-levels` are rejected for doping.

Weak coupling makes narrow kernels expensive to resolve; extreme filling,
huge interactions and extremely small endpoint distances can exhaust the
mesh or native precision. There is no guarantee that every positive `U` and
every representable density converges within these limits. Failed backgrounds
and points expose explicit statuses (`mesh_limit`, `density_limit`,
`linear_failure`, `momentum_limit`, `precision_limit`), omit energies and
return exit status 2. Invalid arguments return 1. Increasing precision does
not replace checking convergence.

## C++ interface and validation

Include `<bethe/hubbard_doped.hpp>`:

```cpp
namespace ht = bethe::hubbard::thermo;
ht::DopedSolver<double> solver(4.0, 0.75);
auto const& background = solver.background();
auto point = solver.at_momentum(ht::DopedBranch::holon, 1.0,
                               ht::Convention::symmetric, ht::EnergyReference::fermi);
```

Check `background.converged` and `point.converged` before reading the optional
energies. Construct a solver per interaction/density/precision; reuse it for
grids. `momentum_range(branch)` gives the supported unwrapped interval.

Tests cover all three precisions, independent 48-decimal full-interval
digamma/log-gamma reference calculations, thermodynamic energy identities,
`mu=de0/dn`, convention shifts, gapless endpoints, the large-U limit, approach
to half filling, invalid inputs, work budgets and native-precision CLI parsing.
The optional reference generator is
[`scripts/reference_hubbard_doped.py`](../scripts/reference_hubbard_doped.py).
