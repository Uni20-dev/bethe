# XYZ thermodynamic spinons and bound branches

`bethe-xyz-dispersion` supplies infinite-chain excitation energies for the
[same XYZ Hamiltonian](xyz.md#convention-and-initial-target) as `bethe-xyz-pbc`.
It does not enumerate finite-ring eigenstates. All arithmetic and exports use
the selected fp64, native long-double, or optional fp128 precision.

## Hamiltonian and scope

```math
\begin{aligned}
H&=J\sum_j(J_xS_j^xS_{j+1}^x+J_yS_j^yS_{j+1}^y+J_zS_j^zS_{j+1}^z),\qquad S=\sigma/2,\\{}
J_x&=\frac{\theta_4(\eta\mid it)}{\theta_4(0\mid it)},\quad
J_y=\frac{\theta_3(\eta\mid it)}{\theta_3(0\mid it)},\quad
J_z=\frac{\theta_2(\eta\mid it)}{\theta_2(0\mid it)},\\{}
\theta_j(u\mid it)&=\vartheta_j(\pi u\mid it),\qquad
0\lt\eta\lt1,\quad t\gt0,\quad J\gt0.
\end{aligned}
```

The family has $`J_x\gt J_y\gt |J_z|`$ at finite interior parameters;
$`J_z`$ changes sign at $`\eta=1/2`$. The field and temperature are zero.
There are two symmetry-breaking x-Néel vacua. A single spinon connects
different asymptotic vacua; a local excitation of a fixed vacuum is not a
single spinon. The two domain-wall orientations share the displayed dispersion.
Generic XYZ has **no conserved total Sz**.

The elementary and bound curves are based on
[Johnson–Krinsky–McCoy (JKM), Sec. VII](https://doi.org/10.1103/PhysRevA.8.2526).
Energies, not spectral weights, are implemented. Arbitrary exchange inversion,
other coupling regions, finite fields, scattering amplitudes and form factors
remain outside this interface.

## Use and output

```sh
bethe-xyz-dispersion --eta 0.4 --t 1 --csv spinons.csv
bethe-xyz-dispersion --eta 0.75 --t 1 --branch bound --json bound.json
bethe-xyz-dispersion --eta 0.875 --t 0.75 --bound-state 2 --precision long-double
bethe-xyz-dispersion --eta 0.75 --t 1 --branch bound --momentum 0
bethe-xyz-dispersion --references
```

`--branch` selects `spinon`, `two-spinon`, `bound`, or `all` (default).
Every strictly existing bound index is included unless `--bound-state S`
selects one. Exactly at a merger the new index is **not** included.
A bound-only request in a region with no bound branch is an input error.
The enumeration limit is 4096 indices, and total output is capped at one
million rows; select an individual index or fewer points for larger requests.
These are output-resource limits, not a truncation disguised as completeness.

The single `dispersion` table supports the common [output options](output.md),
including simultaneous screen, CSV, TSV and JSON sinks, streaming and disabled
retention. It contains:

- `branch`, optional bound index `s`, and relative `delta_rx`, `delta_rz` parities;
- `p`, `p_over_pi`, `cell_momentum`, and the bound curve's `reduced_q`;
- `energy` for lines, or `lower`/`upper` for the two-spinon envelope;
- `status`, with missing rather than invented numerical values on failure.

`--points` defaults to 65 per branch/parity copy; `--momentum` selects one
physical momentum instead. Spinons use $`0\le p\le\pi`$ and other rows use
$`0\le p\le2\pi`$. When spinons are selected, a shared `--momentum` must lie
in their smaller interval. Momentum units are radians per site, not units of pi.
Endpoint-inclusive grids repeat periodic endpoints deliberately.

Metadata records the theta parameters, dimensionless couplings, overall J,
spinon gap/band maximum, conventions, precision, provenance and CPU time.
All excitation energies include J. Numerical range failures exit 2, retain
labels and export JSON null or empty CSV/TSV energy cells. Invalid input exits
1 before opening output files, including with `--force`.

## Convention audit and elementary band

Here is our explicit conversion to JKM's Hamiltonian
$`H=-\frac12\sum_a\mathcal J_a\sigma^a\sigma^a`$:

```math
(\mathcal J_x,\mathcal J_y,\mathcal J_z)
=-\frac J2(J_z,J_y,J_x),\qquad
l=\frac{\theta_2(0\mid it)^2}{\theta_3(0\mid it)^2},\quad
K=\frac\pi2\theta_3(0\mid it)^2,\quad K'=tK,\quad\zeta=K\eta.
```

The axis permutation is global, so it does not change translation momentum.
Using the [Jacobi/theta identities](https://dlmf.nist.gov/22.2), our own
coupling ratios give $`J_z/J_x=\mathrm{cn}(2K\eta,l)`$ and
$`J_y/J_x=\mathrm{dn}(2K\eta,l)`$. Thus JKM's
$`\mu=\pi\eta`$, $`\lambda=\pi\eta/t`$, and $`\tau=\pi/t`$.
Their tau is **not** our imaginary theta period.

Introduce the excitation nome parameter $`t_1=\eta/t`$. In our normalization,
the maximum I and mass m become

```math
\begin{aligned}
I&=\frac{J\pi}{2t}\frac{\theta_1(\eta\mid it)}{\theta_1'(0\mid it)}
       \theta_3(0\mid it_1)^2,\\{}
m&=\frac{J\pi}{2t}\frac{\theta_1(\eta\mid it)}{\theta_1'(0\mid it)}
       \theta_4(0\mid it_1)^2,\\{}
\epsilon(p)&=\sqrt{m^2\cos^2p+I^2\sin^2p},\qquad 0\le p\le\pi.
\end{aligned}
```

The prime differentiates our normalized u, not pi*u. This follows by
substituting the coupling map into JKM (7.8) and applying Jacobi's derivative
identity. The implementation evaluates the small mass directly, never by
subtracting a rounded elliptic modulus from one. Modular transformation and
scaled theta functions preserve exponentially small gaps. A gap below the
normal range of the selected scalar is rejected rather than silently set to zero.

Two-spinon rows use the [shared band kinematics](xxz-dispersion.md#two-spinon-edges)
and take the envelope of both pi-shifted translation copies. This is an
all-sector envelope, **not** a symmetry-resolved threshold or spectral support
for any particular operator. A two-site MPS comparison uses
$`k_{\rm cell}=2p\pmod{2\pi}`$.

## Bound curves and symmetry labels

For positive integer s the strict condition is

```math
s(1-\eta)\lt\eta.
```

For example, eta=0.75 has s=1,2, with s=3 at its excluded merger.
Let $`K_1'/K_1=\eta/t`$, $`r=m/I=k_1'`$, and define

```math
\begin{aligned}
a_s&=\mathrm{sn}\!\left(\frac{sK_1(1-\eta)}{t},k_1'\right),\qquad0\lt a_s\lt1,\\{}
E_s(Q)&=\frac{2I}{a_s}
\sqrt{\sin^2(Q/2)+r^2a_s^2\cos^2(Q/2)}
\sqrt{\sin^2(Q/2)+a_s^2\cos^2(Q/2)},\\{}
E_s(0)&=2ma_s.
\end{aligned}
```

Within a fixed asymptotic vacuum the excitation gap is 2m in the region without
bound branches, and $`E_1(0)`$ when the first bound branch exists. This excludes
the vanishing finite-ring splitting between the two ground-state combinations;
it is distinct from the reported **single-spinon** gap m.

We fix the positive energy prefactor directly from JKM (7.11a), and verify
the eliminated expression against the complex-rapidity parametrization in
the independent reference script below. The two square roots are important:
replacing m in the spinon curve with a bound-state mass gives the wrong
lattice dispersion. At a merger $`a_s\to1`$, the curve tends to
$`2\epsilon(Q/2)`$, not necessarily the lowest continuum edge at that Q.
A bound curve can run above scattering states at other momenta or in other
channels; “bound” does not promise a line below the all-sector envelope everywhere.

For even rings define $`R_a=\prod_j\sigma_j^a`$. The library and table use
**ratios to the chosen reference vacuum's parity**, not absolute ring labels:

```math
\delta r_x=(-1)^s,\qquad\delta r_z=\pm1,\qquad
Q=p-\pi\frac{1-\delta r_z}{2}\pmod{2\pi}.
```

This is the axis-mapped version of JKM's $`\Delta\nu'=s\bmod2`$ and
$`\Delta P=\Delta Q+\pi\Delta\nu''`$. Both z-parity copies are output;
their reduced dispersions coincide but their physical momenta differ by pi.
The two-site translation phases therefore coincide. Single-spinon rows leave
global ring parities empty: their two ends belong to different vacua.
We do not infer a finite-size degeneracy or state count from these curves.

## Library and validation

`bethe::xyz::SpinonDispersion<Real>(eta,t,J=1)` provides `gap()`,
`maximum_energy()`, `energy(p)`, `continuum(q,convention)`, `bound_exists(s)`,
`bound_x_parity(s)`, `bound_gap(s)`, and `bound_energy(s,Q)`.
Unlike the frontend envelope, the library continuum defaults to a fixed sum
of two spinon momenta; request `SpinonMomentum::folded` for the union.
The bound library method takes **reduced Q**, while the frontend's
`--momentum` takes physical p. Domain mistakes throw `invalid_argument`;
numerical failures throw a `runtime_error` subtype.

The independent oracle can be rerun with mpmath (not a build dependency):

```sh
python3 scripts/reference_xyz_dispersion.py --eta 0.75 --t 1 --digits 65
```

Checks include native-precision references, strict merger conditions,
exchange scaling, tiny masses, explicit underflow, and these independent limits:

- XY at eta=1/2: $`I=(J_x+J_y)J/2`$, $`m=(J_x-J_y)J/2`$.
- Gapless XXZ as t grows: $`I\to J\sin(\pi\eta)/(2\eta)`$.
- Massive AF XXZ as t decreases with $`\eta=t\lambda/\pi`$, after axis
  permutation and rescaling by Jy: $`\Delta=\cosh\lambda`$.
- Ferromagnetic XXZ with $`\eta=1-t\lambda/\pi`$: after the corresponding
  staggered rotation and rescaling, the s-magnon curve is
  $`\sinh\lambda[\cosh(s\lambda)-\cos Q]/\sinh(s\lambda)`$.

Literal spin-Hamiltonian diagonalization on N=6,8,10 checks the first bound
branch and both momentum/parity copies at eta=0.75,t=1. Their gaps approach
the thermodynamic 0.24378676516 from opposite sides; at N=10 they are
0.23078578101 and 0.25686154761. The finite-size splitting is real, not a
solver tolerance. This tests normalization and momentum assignment without
claiming that a thermodynamic formula equals a small-ring eigenvalue.
