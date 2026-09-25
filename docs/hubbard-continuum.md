# Half-filled Hubbard continuum edges

**Status: two-spinon and three charge-containing families implemented in
the library and the shared frontend.**

This extends the [elementary dispersion lines](hubbard-dispersion.md) to a
specified scattering family, not to the entire spectrum at fixed quantum
numbers. The current scope is the infinite half-filled repulsive chain,
zero magnetic field, hopping t=1, U>0.

```sh
bethe-hubbard-continuum --u 4 --points 33 --csv edges.csv
bethe-hubbard-continuum --u 4 --momentum -1 --precision long-double
bethe-hubbard-continuum --u 4 --points 9 --precision fp128 --json edges.json
bethe-hubbard-continuum --references
```

The default grid has 33 equally spaced total momenta on [-pi,pi], including
both endpoints (the same Brillouin-zone point). `--momentum` selects one point
and excludes `--points`. The `two_spinon` table reports both energies,
constituent momenta, separate error estimates, work counts and status.
`--tolerance`, `--max-evaluations`, `--max-levels` and `--max-iterations`
control the native solver. Precision defaults to fp64; fp128 requires MPLAPACK.

The [common output options](output.md) support screen output and additional
JSON/CSV/TSV files. Exports stream each completed point; `--no-retain`
discards delivered rows. Metadata includes conventions and numerical controls,
with CPU time and overall outcome in the final summary. Failed points have
missing energies (JSON null, empty CSV/TSV cells) and yield exit status 2;
other requested points are still attempted. Invalid arguments exit with
status 1 before output files are opened, even with `--force`.

## Two spinons

[Essler–Korepin, Eq. (18) and the following paragraph](https://arxiv.org/html/cond-mat/9808018)
give additive two-spinon energies and momenta, an endpoint spinon at the
lower boundary, and equal rapidities at the upper boundary. For a requested
total momentum P in [-pi,pi], put q=|P|. In terms of the existing spinon line:

```math
\begin{aligned}
E_{\mathrm{lower}}(P)&=\epsilon_s(q),\\
E_{\mathrm{upper}}(P)&=2\epsilon_s(q/2).
\end{aligned}
```

For P>=0, constituent momenta are (0,q) and (q/2,q/2).
For P<0 they are (pi,pi-q) and (pi-q/2,pi-q/2), respectively;
their sums equal P modulo 2*pi. Endpoints denote limiting scattering states.
No finite-ring multiplicity or spectral weight is implied.

Both edges have particle change zero. Consequently the symmetric and
unshifted Hubbard conventions, and Hamiltonian and Fermi references, give
the same excitation energies. The upper edge bounds **two** spinons only;
it is not an upper bound on all multiparticle excitations.

## API and numerical contract

```cpp
#include <bethe/hubbard_continuum.hpp>
auto edges = bethe::hubbard::thermo::two_spinon_continuum(4.0, 1.0);
if (edges.converged) {
  auto lower = *edges.lower_energy;
  auto upper = *edges.upper_energy;
}
```

The implementation reuses `hubbard_thermo.hpp` without another quadrature
implementation or a numerical optimizer. Native fp64, long-double and
fp128 follow the same path. `Options<Real>` has the elementary-line semantics,
except its evaluation budget is shared across both edges; inversion updates
remain limited per constituent point. Work counts include both solves.
At P=0 both edges are analytically zero and need no quadrature.

Neither energy is published if either solve fails. The underlying
`quadrature_limit`, `momentum_limit` or `precision_limit` status is retained.
Energy quadrature errors and total-momentum errors are separate: quadrature
error does **not** include propagation of momentum inversion error. Unknown
errors remain infinite. These are numerical estimates, not certified bounds.

Regression tests compare with independent Fourier–Bessel spinon values,
verify momentum witnesses and interior scattering energies, check the
large-U Heisenberg limit, and exhaust the shared budget after the lower-edge
solve to ensure a partial pair of edges is not published.

## Charge-containing channels

Select `--channel spinon-holon`, `spinon-antiholon` or `holon-antiholon`;
the default is `two-spinon`. Their particle changes are respectively -1,
+1 and 0, and their total spins are 1/2, 1/2 and 0. For example:

```sh
bethe-hubbard-continuum --u 4 --channel spinon-holon --momentum 0
bethe-hubbard-continuum --u 4 --channel holon-antiholon --points 9 --csv charge.csv
bethe-hubbard-continuum --u 4 --channel spinon-antiholon --momentum 0 \
  --convention unshifted --reference fermi
```

`--convention symmetric` (default) uses the SO(4)-symmetric interaction;
`unshifted` adds U*DeltaN/2 to the symmetric excitation energy.
`--reference hamiltonian` (default) reports that Hamiltonian's energy
change. `--reference fermi` subtracts the half-filled chemical potential
times DeltaN, recovering the symmetric energy for either convention.
Negative unshifted removal energies are therefore not a negative symmetric
gap. These switches do not change two-spinon energies.

The `charge_continuum` table includes selected and symmetric energies,
constituent momentum witnesses, separate vertical and horizontal error
estimates, quadrature and search counters, and distinct `search_status`
and `constituent_status`. All energies, witnesses and errors are missing
on failure. An `objective_failure` means a constituent failed; inspect
`constituent_status` for its reason. A search budget failure can occur
even when every sampled constituent converged.

Charge channels reuse the native elementary-line solver and a shared
bounded-extrema search. Controls are separated by layer:

- `--tolerance`, `--max-evaluations` (default 1000000), `--max-levels`
  and `--max-iterations` apply per constituent. Unlike two spinons,
  charge searches need many constituent evaluations.
- `--max-total-evaluations` (default 200000000) caps all quadrature work
  for one requested total momentum, across both edges.
- `--search-tolerance` defaults to 65536 epsilon and applies to E/max(1,U).
  `--position-tolerance` defaults to sqrt(epsilon), multiplied by the
  search domain's coordinate scale (pi) for the local bracket width.
- `--max-search-evaluations` (20000), `--max-search-iterations` (256 per
  local bracket), `--initial-intervals` (16) and `--max-intervals` (128)
  limit the extrema search. Mesh sizes must satisfy 4<=initial<=max<=4096.

Search-specific options are rejected for the two-spinon channel rather
than silently ignored. Budgets reset at each requested total momentum.
Native fp128 charge searches can be expensive: a default-tolerance U=4,
P=0 holon–antiholon point took about 15 minutes in GCC 13 Debug. Start
with a single point; use explicitly looser tolerances for exploratory grids,
for example `--tolerance 1e-10 --search-tolerance 1e-7 --position-tolerance 1e-5`.
No lower-precision fallback is used.

The [charge-continuum library and independent oracle](hubbard-charge-continuum-design.md)
now cover spinon–holon, spinon–antiholon and holon–antiholon families.
Their channels, conventions, search limitations and reference points are
documented there, including the independent reference calculations.

The charge-containing library searches the corresponding sums at fixed
total momentum, refining candidate stationary branches and endpoints and
propagating constituent failures. Its mesh-convergence evidence is not a
global certificate, nor is a chosen particle family's lower edge a minimum
over arbitrary additional particles. Doped continua and spectral weights
remain separate work.
