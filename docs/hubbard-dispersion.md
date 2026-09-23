# Hubbard spinon and charge dispersions

`bethe-hubbard-dispersion` evaluates elementary excitation lines of the infinite,
repulsive Hubbard chain at zero field, at or below half filling. Unlike the
[finite-ring](hubbard.md) and [free-end](hubbard-open.md) solvers, it takes no
length or boundary-condition argument. It requires `U>0`, hopping `t=1`, and
`0<n<=1`. The default is half filling (`--density 1`), described below.
For `--density n` with `n<1`, see the [doped dispersion guide](hubbard-doped.md).
Attraction, nonzero magnetic fields, spectral weights, and multiparticle
continuum thresholds are not implemented here.

```sh
build/bethe-hubbard-dispersion --u 4 --points 101 --format csv
build/bethe-hubbard-dispersion --u 4 --branch spinon --momentum 1
build/bethe-hubbard-dispersion --u 4 --convention unshifted --format tsv
build/bethe-hubbard-dispersion --u 4 --precision fp128 --points 33 --format csv
```

Calculations and output retain the selected `fp64` (default), `long-double`, or
optional `fp128` precision. Grids are uniform in **dressed physical momentum**,
not the auxiliary Bethe parameter. `--points` includes both endpoints; charge
endpoints `-pi` and `pi` are the same point in the Brillouin zone. `--momentum`
instead selects one point in radians. With `--branch all` it must belong to
every selected branch.

Use `--help` or `-h` for Uni20-rendered options and examples, and `--references`
for the literature and its applicability. `--version` and `--build-info`
report the compiled application and dependency
identity without calculating or opening files. This frontend accepts both
`--u 4` and `--u=4`, rejects repeated scalar options, and retains repeatable
file-export options. See [shared CLI behavior](command-line.md#uni20-help-and-option-parsing).

## Which excitation?

| Branch | Particle change | Total spin | Momentum interval |
|---|---:|---:|---|
| `spinon` | 0 | 1/2 | `[0,pi]` |
| `holon` | -1 | 0 | `[-pi,pi]` |
| `antiholon` | +1 | 0 | `[-pi,pi]` |

“Chargon” describes the charge excitations, not a fourth independent branch.
Spin projection does not split the spinon line at zero field. At the same bare
charge parameter, holon and antiholon have equal symmetric energies but
`p_antiholon = p_holon - pi (mod 2*pi)`.

Spinons are gapless at both endpoints. Holon and antiholon minima occur at
`p=-pi/2` and `p=pi/2` respectively. The common symmetric minimum is **half**
the charge gap `E(N+1)+E(N-1)-2E(N)`. At `U=4` this single-charge minimum is
approximately `0.64336351100645219736629488464487225`.

## Hamiltonian and energy zero

The default `--convention symmetric` uses the SO(4)-symmetric Hamiltonian

```text
H_symmetric = -sum(c†_j,s c_(j+1),s + h.c.)
              + U sum (n_j,up - 1/2)(n_j,down - 1/2).
H_symmetric = H_unshifted - U*N/2 + U*L/4.
```

`--convention unshifted` matches `U*n_up*n_down` in our finite Hubbard tools.
For excitation energies above the same background at fixed length,

```text
E_unshifted = E_symmetric + U*DeltaN/2.
```

The extensive background cancels. Spinons are unchanged, holons shift down by
`U/2`, and antiholons up by `U/2`. Negative unshifted removal energies are
legitimate: half filling is centered at chemical potential `U/2` in that
Hamiltonian, not at zero chemical potential.

The default `--reference hamiltonian` reports these Hamiltonian differences.
`--reference fermi` instead reports `E_H-mu_H*DeltaN`, independent of the chosen
interaction convention. At half filling we choose the middle of the Mott
plateau: `mu_symmetric=0`, `mu_unshifted=U/2`. Thus the Fermi-referenced energies
equal the symmetric energies here, including the nonzero charge gap. Both
chemical potentials appear in the metadata. For doping they are solved from
the requested density, not fixed at these half-filled values.

## iMPS comparisons

These are **elementary lines**, not necessarily sector minima at every momentum.
A holon with two spinons can lie below the elementary holon at the same total
momentum. No continuum minimization is performed here.

Fractional excitations need an appropriate domain-wall/topological iMPS ansatz;
they are not arbitrary isolated finite-ring states. Match hopping, interaction
convention, particle change and spin first. Our momenta are one-site momenta.
An enlarged unit cell folds the Brillouin zone, and a domain-wall gauge can
shift the momentum origin; no universal shift is imposed. The
[MPSKit Hubbard example](https://quantumkithub.github.io/MPSKit.jl/stable/examples/quantum1d/6.hubbard/#Excitations)
illustrates the comparison.

## Method and accuracy

The reference formulas are the Fourier-Bessel dispersions in
[Essler–Korepin, Eqs. (4)–(5)](https://arxiv.org/html/cond-mat/9808018).
Set `u=U/4` in that paper's notation. Our spin rapidity is
`lambda=2*u*beta/pi`; `k` denotes bare charge momentum.

Direct oscillatory integration is costly at weak coupling, and subtracting
order-one energies can erase the exponentially small charge gap. Instead we
use a nonoscillatory spinon convolution with `sech` and a positive-integrand
resummation of the exact modified-Bessel series in
[Melzer](https://arxiv.org/html/cond-mat/9410043). His attractive massive spin
branch maps to our repulsive charge branch: his interaction parameter becomes
our `u`, and energies multiply by two to change hopping from `1/2` to `1`.

For reproducibility, set `a=pi/(2*u)` and `A(z)=atanh(exp(-z))`. Our resummation is

```text
E_inner(k) = (4/pi) integral_0^infinity cosh(t)
             * [A(a*(cosh(t)-sin(k))) + A(a*(cosh(t)+sin(k)))] dt
Q(k)       = (2/pi) integral_0^infinity
             [A(a*(cosh(t)-sin(k))) - A(a*(cosh(t)+sin(k)))] dt
E_h(k)     = E_inner(k) + 4*max(cos(k),0).
```

Unwrapped holon momentum is `pi/2-2*k+Q` for `|k|<pi/2`, `-pi/2+Q` for
`k>=pi/2`, and `3*pi/2+Q` for `k<=-pi/2`. Stable trigonometric differences
retain the integrable logarithmic endpoint at `k=+-pi/2`.

Tanh-sinh quadrature compares successive meshes; infinite-domain tails are
bounded separately. Reported errors are **estimates**, not interval certificates:

- `energy_quad_error` estimates quadrature, tail and shift-roundoff error at
  the returned parameter. It excludes propagation of momentum-inversion error.
- `momentum_error` includes the dressed-momentum residual and quadrature error.
  The spinon target scales with distance to the nearest gapless endpoint;
  the charge target is absolute in radians.
- `parameter` is rapidity or bare momentum, not plotted momentum. Analytic
  spinon endpoints have infinite rapidity and zero energy.

Defaults: tolerance `256*epsilon` of the selected scalar, 1,000,000 quadrature
samples **per entire point**, 12 refinement levels and 160 inversion updates.
Use `--tolerance`, `--max-evaluations`, `--max-levels`, and `--max-iterations`
to change them. Weak coupling/fp128 can need substantially more work. Gaps and
endpoint distances below the scalar's representable range are not promised.

`quadrature_limit`, `momentum_limit`, or `precision_limit` mark failed points:
their energies are omitted, never replaced with a claimed converged zero.
Exit status is 0 for complete success, 2 for failed points, and 1 for invalid
arguments, runtime errors or output failures. Inspect status before increasing
a budget or changing precision.

## Output and C++ API

`auto`, `pretty`, and `plain` use Uni20 presentation. Add `--csv FILE`,
`--tsv FILE` or `--json FILE` to export alongside the screen report, or use
`--quiet` for files only. `--format` selects stdout, including machine formats.
The [output guide](output.md) explains streaming, memory retention, provenance,
overwrite protection, and precision-preserving JSON.

CSV/TSV have one header row, unwrapped native-precision tokens, decimal spins
(`0.5`), empty failed-energy cells, and initial `#` metadata. Status and solver
CPU time are **trailing** comments. `--no-preamble` produces strict CSV/TSV.
JSON has typed column metadata, decimal strings for real values, numerical
half-integers, and `null` for missing values.
`energy` follows the selected convention/reference;
`symmetric_energy` always means the symmetric **Hamiltonian** difference and
`fermi_energy` always means the chemical-potential-subtracted energy.
Citations appear in `--references` and [CITATIONS.md](../CITATIONS.md), not ordinary help or numeric output.

Include `<bethe/hubbard_thermo.hpp>` and use `bethe::hubbard::thermo`:

```cpp
auto point = dispersion(Branch::holon, 4.0, -1.0, Convention::symmetric);
auto spin = spinon_at_rapidity(4.0, 0.5);
auto charge = charge_at_bare_momentum(Branch::antiholon, 4.0, 1.0);
```

`Point<Real>` contains optional energies, exact particle change and half-integer
spin, parameters, error estimates, work counts and status. Tests use independent
Fourier-Bessel references at all three precisions, a weak-coupling relative-gap
reference, convention shifts, symmetry, budgets and the large-U XXX limit.
