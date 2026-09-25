# Half-filled Hubbard continuum edges

**Status: two-spinon library and frontend implemented; charge-containing continua pending.**

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

```text
E_lower(P) = epsilon_s(q)
E_upper(P) = 2*epsilon_s(q/2).
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

## Subsequent slices

Spinon–holon/antiholon and holon–antiholon thresholds require constrained
minimization of the corresponding sums at fixed total momentum. That work
must inspect all stationary branches and endpoints, propagate constituent
failures, and distinguish the chosen particle family from a global sector
minimum. Doped continua and spectral weights remain separate work.
