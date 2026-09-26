# Thermodynamic XXZ spinons and continua

[Back to the overview](../README.md)

`bethe-xxz-dispersion` supplies reference curves for infinite-chain excitation
calculations, including symmetry-resolved MPS and iMPS ansätze. Unlike
[`bethe-xxz-pbc`](xxz.md), it does not solve a finite ring. It supports zero
field and zero temperature, positive exchange, and anisotropy $`\Delta\gt -1`$:

```math
H=J\sum_j\left(S_j^xS_{j+1}^x+S_j^yS_{j+1}^y+\Delta S_j^zS_{j+1}^z\right),
\qquad S^\alpha=\sigma^\alpha/2,\quad J\gt 0.
```

There is no subtraction of the ferromagnetic constant. The overview gives
the bulk ground-state energy **per site**; the plotted energies are excitation
energies **above** that ground state, not extensive energies.

## Run and export

```sh
build/bethe-xxz-dispersion --delta 2 --branch spinon --csv spinons.csv
build/bethe-xxz-dispersion --delta 2 --branch two-spinon --folded --tsv continuum.tsv
build/bethe-xxz-dispersion --delta 1.0078125 --precision long-double --json curves.json
build/bethe-xxz-dispersion --delta 0.5 --exchange 2 --momentum 1
build/bethe-xxz-dispersion --references
```

By default `--branch all` produces both families in the named `dispersion`
table. `--points` sets the number of equally spaced samples **per family**;
`--momentum` instead requests one momentum and excludes `--points`.

- `spinon`: `p` ranges from 0 to π; `energy` is populated.
- `two-spinon`: `p` denotes the total momentum Q, ranging from 0 to 2π;
  `lower` and `upper` are populated. `energy` is empty.
- `p_over_pi` is convenient for plotting; `cell_momentum` is the phase of
  translation by **two sites**, in [0,2π), not momentum in radians/site.

Metadata records the Hamiltonian, normalization, sectors, continuum convention,
precision, numerical controls and bulk result. Use `--branch` to get a single
family in a CSV file. Screen and file sinks, `--quiet`, `--stream`,
`--no-retain`, and overwrite protection follow the [common output contract](output.md).
The bibliography appears only with `--references`.

## Which excitation is an MPS calculation seeing?

In the massive antiferromagnet, $`\Delta\gt 1`$, there are two Néel vacua related
by translation by one site. A single spinon has $`S^z=\pm1/2`$ and interpolates
between these different asymptotic vacua. It belongs to a **topological**
excitation ansatz, not to the local excitation space above a fixed vacuum.
Two spinons return to the original asymptotic vacuum and can have
$`S^z=-1,0,+1`$. Their energies do not depend on these spin orientations at
zero field. These are Sz labels, not SU(2) total-spin labels at generic Δ.

A two-site MPS unit cell has translation eigenphase

```math
q_{\mathrm{cell}}=2p\pmod{2\pi}.
```

Thus p and p+π become indistinguishable under two-site translation. The
single-spinon interval [0,π] already spans this reduced Brillouin zone;
its endpoints are equivalent under this translation. One-site translation
exchanges the massive vacua, so assigning a one-site eigenmomentum to a
single domain wall requires a choice of translation-sector combinations.

For two spinons, the default **unfolded** table keeps the constituent
convention $`Q=p_1+p_2\pmod{2\pi}`$. With `--folded`, each Q includes the
union of this continuum and its π-shifted image. This is the appropriate
energy envelope for a two-site-cell comparison. We retain the full plotting
interval [0,2π], so the folded results repeat with period π.

For massive translation-sector combinations the one-site eigenmomentum is
$`P=p_1+p_2+s\pi\pmod{2\pi}`$, with $`s=0,1`$; the unfolded output specifies
the s=0 constituent branch, not the union of all translation sectors.
In the gapless regime, folding is simply Brillouin-zone folding and does
not imply two symmetry-broken Néel vacua.

Folding energies does **not** make a structure factor π-periodic: the two
translation branches can have different matrix elements. In particular,
some weights vanish in the isotropic limit. The tool computes kinematics,
not intensities, observable selection rules, or the complete multiparticle
spectrum. See [Caux–Mossel–Pérez Castillo](../CITATIONS.md#caux-mossel-perez-castillo-2008)
and the [topological MPS excitation framework](../CITATIONS.md#zauner-stauber-2018).

## Single-spinon dispersion

In the gapless regime, $`-1\lt\Delta\le1`$, the existing sine dispersion is reused:

```math
\epsilon(p)=v\sin p,\qquad
v=\frac{J\pi\sin\gamma}{2\gamma},\qquad
\Delta=\cos\gamma,\quad 0\le p\le\pi.
```

At Δ=1 the limiting coefficient is Jπ/2. In the massive regime set
$`\Delta=\cosh\eta`$ and let k be the elliptic **modulus**, not the parameter
k². With $`K'=K(\sqrt{1-k^2})`$:

```math
\frac{K'}{K}=\frac{\eta}{\pi},\qquad
I=\frac{JK\sinh\eta}{\pi},\qquad m=I\sqrt{1-k^2},\qquad
\epsilon(p)=\sqrt{m^2\cos^2p+I^2\sin^2p}.
```

Here m is the **single-spinon** gap; the minimum two-spinon energy is 2m.
The maximum single-spinon energy is I. As Δ→∞ at fixed J,
$`m/J=\Delta/2-1+1/(4\Delta)+O(\Delta^{-2})`$ and
$`I/J=\Delta/2+1+1/(4\Delta)+O(\Delta^{-2})`$.

Near Δ=1 the gap is exponentially small. Numerically we evaluate theta
constants with nome $`e^{-\eta}`$ or its modular transform
$`e^{-\pi^2/\eta}`$, reusing the native-precision elliptic infrastructure.
We **never** obtain the complementary modulus by subtracting rounded k²
from one. Endpoint sine values are handled explicitly so a rounded sin(π)
does not swamp the gap.

## Two-spinon edges

The continuum follows by minimizing and maximizing
$`\epsilon(p_1)+\epsilon(p_2)`$ at fixed $`p_1+p_2=Q`$ with each momentum in
[0,π]. Reflection supplies Q in [π,2π], so consider [0,π]. Define

```math
A(Q)=2\epsilon(Q/2),\qquad B(Q)=m+\epsilon(Q),\qquad
Q_* = 2\arctan\sqrt{m/I}.
```

For the gapless case use m=0 and I=v. The unfolded edges are

```math
E_{\mathrm{lower}}(Q)=
\begin{cases}
A(Q),&0\le Q\le Q_*,\\{}
(I+m)\sin Q,&Q_*\le Q\le\pi/2,\\{}
B(Q),&\pi/2\le Q\le\pi,
\end{cases}
\qquad E_{\mathrm{upper}}(Q)=\max\{A(Q),B(Q)\}.
```

The lower-edge stationary point changes from equal momenta, to an unequal
pair, to an endpoint pair. These formulas are derived directly from the
two-particle energy and tested against independent momentum scans. In
particular, the last lower-edge branch is **m+ε(Q)**: do not substitute
ω₊(Q) from the piecewise expression printed after Eq. (33) of
arXiv:0806.3069v1, which disagrees with its accompanying endpoint description
and with direct minimization. The upper edge can be evaluated as a maximum
without solving the paper's crossing quartic.

The folded lower/upper bounds are respectively the minimum/maximum of
these bounds at Q and Q+π. All of this is shared elliptic-band kinematics,
separate from the XXZ coupling-to-band conversion.

## Bulk energy and numerical limits

The zero-field bulk formulas, with our constant restored, are

```math
\frac{e_\infty}{J}=\frac{\Delta}{4}
-\sin\gamma\int_0^\infty
\frac{\sinh((\pi-\gamma)x)}{\sinh(\pi x)\cosh(\gamma x)}\,dx
\quad(-1\lt\Delta\lt 1),
```

```math
\frac{e_\infty}{J}=\frac{\Delta}{4}
-\sinh\eta\left(\frac12+2\sum_{n=1}^\infty\frac{1}{1+e^{2n\eta}}\right)
\quad(\Delta\gt 1).
```

We use explicit XX/XXX limits and the rapidly converging series for η≥1.
For η<1, Poisson resummation gives

```math
\rho(\lambda)=\frac{1}{2\eta}\sum_{r\in\mathbb Z}
\mathrm{sech}\!\left[\frac{\pi}{\eta}(\lambda-r\pi)\right].
```

Integration of the bare energy against this density, with λ=ηx, remains
well-conditioned and inexpensive near Δ=1. The gapless integral similarly
rescales its variable. Tanh–sinh mesh comparisons, analytic tail bounds and
a roundoff allowance provide an **error estimate**, not a certified interval.

`--tolerance` is an absolute target in units J max(1,|Δ|), defaulting to 256
machine epsilon. Its allowed range is [128 epsilon, 0.01]. `--max-evaluations`
bounds bulk integrand evaluations or Fourier terms; `--max-levels` bounds
quadrature refinements. Exact XX/XXX limits do not consume these budgets.

fp64, native long double and optional MPLAPACK fp128 retain their native
arithmetic. A massive gap below the scalar type's **normal** range produces
`precision_limit`, not a spurious zero gap. Long double and fp128 usually
extend that exponent range on x86; even these cannot represent every gap
arbitrarily close to Δ=1. Bulk energy can still succeed there. Conversely,
exhausting the bulk budget does not invalidate independently computed
analytic dispersion rows. Either failure gives exit status 2; missing table
values are JSON null or empty CSV/TSV cells.

## Library interface

```cpp
#include <bethe/xxz_dispersion.hpp>
#include <bethe/xxz_bulk.hpp>

bethe::xxz::SpinonDispersion<long double> band(2.0L);
auto gap = band.gap();
auto energy = band.energy(1.0L);
auto unfolded = band.continuum(1.0L);
auto folded = band.continuum(1.0L, bethe::SpinonMomentum::folded);
auto bulk = bethe::xxz::bulk_energy_density(2.0L);
if (bulk.converged) { /* *bulk.energy, *bulk.error */ }
```

Invalid inputs throw `std::invalid_argument`. The dispersion constructor
throws `std::underflow_error` for an unrepresentable massive gap and
`std::overflow_error` for overflowing energies. Bulk nonconvergence returns
absent energy/error fields. The older gapless-only `spinon_energy` function
in `<bethe/spinon.hpp>` retains its original domain and behavior.

`BulkEnergy::status` distinguishes `converged`, `evaluation_limit`,
`mesh_limit`, and `precision_limit`; the same labels appear in run metadata.

Finite field, dressed finite-magnetization dispersions, bound-string
branches, finite-size quantization and spectral weights are outside this
milestone. Ground-state formulas are recorded in
[Bortz–Göhmann](../CITATIONS.md#bortz-gohmann-2005); provenance and the full
bibliography are in [CITATIONS.md](../CITATIONS.md).
