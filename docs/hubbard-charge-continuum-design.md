# Charge-containing Hubbard continuum library

**Status: native charge-continuum library and CLI implemented.** The existing
[continuum tool](hubbard-continuum.md) selects these families with `--channel`.

## Channels and energy reference

Use the half-filled, zero-field repulsive chain with t=1 and our U>0.
[Essler–Korepin (1994)](../CITATIONS.md#essler-korepin-1994-scattering),
Eqs. (4)-(6) and the two-particle construction on preprint pp. 5-6,
give additive constituent energies and momenta. Their interaction is our U/4.

| Channel | Delta N | Total spin S | First momentum domain |
| --- | ---: | ---: | --- |
| spinon–holon | -1 | 1/2 | [0,pi] |
| spinon–antiholon | +1 | 1/2 | [0,pi] |
| holon–antiholon | 0 | 0 | [-pi,pi] |

For constituents a,b and total momentum P modulo 2*pi, minimize/maximize
`epsilon_a(p)+epsilon_b(wrap(P-p))` over the first domain. These are family
edges, not minima over arbitrary extra particles. They do not determine
spectral weights. Holon/antiholon energies agree at the same bare parameter,
but their momenta differ by pi. The two mixed-spin channels therefore map
under a pi shift of P in the symmetric energy convention.

Return symmetric energies first; apply `U*DeltaN/2` only for unshifted
Hamiltonian energies. At half filling, the middle-of-gap Fermi reference
agrees with symmetric excitation energies. Conventions never alter the
minimizing momenta within a fixed channel.

## Native library contract

```cpp
#include <bethe/hubbard_charge_continuum.hpp>
namespace hubbard = bethe::hubbard::thermo;
auto edges = hubbard::charge_continuum(
    hubbard::ChargeChannel::spinon_holon, 4.0, 1.0);
```

`ChargeChannel` selects `spinon_holon`, `spinon_antiholon` or
`holon_antiholon`. Inputs are finite U>0 and total momentum in [-pi,pi].
`Convention::symmetric` is the default; `Convention::unshifted` changes only
the reported Hamiltonian energies. Both symmetric and selected energies
are returned, with particle change and physical spin. Fermi-referenced
energies at half filling are the symmetric values.

`ChargeContinuumOptions<Real>` keeps three work controls distinct:

- `search`: the shared extrema options, acting on energy divided by
  `max(1,U)`. Its default value tolerance is `65536*epsilon`. Multiply
  by `max(1,U)` to obtain the corresponding physical energy tolerance.
- `constituent`: the existing dispersion solver's controls, including its
  per-point quadrature limit and relative momentum target.
- `max_quadrature_evaluations`: cumulative quadrature samples across the
  entire edge pair, default 200000000. Cached constituent points consume
  no additional quadrature budget.

Constituent spin reflection and the charge-doublet momentum shift are used
to reuse evaluations. Neither an interpolation table nor a double-precision
fallback is involved. Full native fp128 searches can be substantially more
expensive than a single elementary-line point. The default-tolerance U=4,
P=0 holon–antiholon regression took about 15 minutes in the GCC 13 debug
build used during development; this is not a cheap momentum-grid operation.
For exploratory plots, use fp64 or explicitly choose appropriate tolerances.

Both energies remain absent on any search or constituent failure. Inspect
`search_status` and `constituent_status` separately: an objective failure
can originate in quadrature, momentum inversion or representability.
Counters distinguish objective calls, local-search updates, mesh sizes,
quadrature samples and constituent inversion updates.

Successful results include both constituent momentum pairs. `lower_error`
and `upper_error` are vertical search/quadrature estimates in physical energy
units, including convention-shift roundoff. They **exclude propagation of
momentum-inversion error**; separate horizontal errors accompany each edge.
They are not rigorous global error bounds. A loose constituent tolerance
cannot be repaired merely by requesting a tighter outer search tolerance.

## Numerical design

Reuse the native elementary dispersion evaluator, not interpolated fp64
tables. A common one-dimensional extrema utility should own sampling,
bracketing and refinement, with the model supplying the objective and domain.
It must support both endpoints, periodic seams, multiple local extrema,
equal competing minima and failed objective samples.

Interior stationary points satisfy equal constituent velocities. An
endpoint-only prescription is insufficient, as the oracle below demonstrates.
A single bounded minimization over the whole domain assumes unimodality and
is also insufficient. Locate candidates on successively refined meshes,
refine every detected extremum and endpoint neighborhood, and require stable
edge energies across independent refinements. Report this as numerical
convergence evidence, not a rigorous global certificate. In particular,
agreeing meshes cannot prove the absence of an arbitrarily narrow unseen
well. Tests must deliberately include competing extrema and branch changes.

Keep independent limits for objective calls, mesh refinements, local
iterations and cumulative constituent quadrature work. Budget exhaustion or
any unresolved constituent must not publish a claimed converged edge.
Both requested edges should follow the two-spinon API's all-or-nothing
contract. Return witnesses, mesh/refinement diagnostics and constituent
quadrature/momentum errors; do not silently treat a horizontal momentum error
as a vertical energy bound. Error floors must prevent endless refinement
below the precision actually supplied by the constituent solves.

### Shared helper checkpoint

`bethe/detail/extrema.hpp` now supplies `bounded_extrema<Real>` independently
of the Hubbard model. Its callback returns an optional `ObjectiveSample<Real>`
containing value and estimated absolute uncertainty. It caches samples and
aborts on missing/nonfinite values or invalid uncertainty. The callback must
be deterministic at a fixed argument; model-specific work is counted by its
owner, while the helper counts uncached objective calls.

Defaults are absolute value tolerance `65536*epsilon`, coordinate tolerance
`sqrt(epsilon)*max(1,|lo|,|hi|)`, 16 initial intervals, at most 128 intervals,
20000 objective calls and 256 golden-section updates per local bracket.
Meshes double; two successive agreements require at least three meshes.
Local refinement requires both coordinate resolution and a small value
spread, with input uncertainty below `tolerance/16`. Endpoints and their
adjacent cells are always considered. Strict improvement on at least one
side avoids refining every point of an exactly flat mesh as a new extremum.

Both optional extrema remain absent on objective, evaluation, iteration,
mesh or precision failure. The selected witnesses have value-error estimates
combining local spread, input uncertainty and mesh-to-mesh change. These
remain **heuristic search estimates**, not certified global bounds or a
guarantee of accurate witness coordinates for nearly degenerate minima.
Tests cover native precision, non-grid minima, multiple periodic extrema,
the domain seam, competing wells, a constant function, uncertainty floors
and exhausted budgets. The native Hubbard wrapper now supplies this objective;
the shared continuum frontend exposes its controls, statuses and missing values.

## Independent developer oracle

`scripts/reference_hubbard_continuum.py` uses SciPy's direct Fourier–Bessel
quadrature, root inversion and bounded scalar optimization, independently of
the native nonoscillatory integration. It is explicitly double precision,
restricted to U>=1, and is not a production dependency or an fp128 oracle.
Its grid search is exploratory evidence, not a completeness proof.

```sh
python3 scripts/reference_hubbard_continuum.py --u 4 --intervals 32 64
python3 scripts/reference_hubbard_continuum.py --u 4 --intervals 64 128 \
  --channels spinon-antiholon --momenta-over-pi -1 -0.75 -0.5 0
```

At U=4, 32- and 64-interval searches agree within 1e-9 in edge energy:

| Channel | P/pi | Lower | Upper |
| --- | ---: | ---: | ---: |
| spinon–holon | 0 | 1.8856450004 | 3.4996123276 |
| spinon–holon | 0.25 | 1.4739946797 | 4.4768844283 |
| spinon–holon | 0.5 | 0.6433635110 | 5.2424524206 |
| spinon–holon | 1 | 3.4996123276 | 5.8856450004 |
| holon–antiholon | 0 | 1.2867270220 | 9.2867270220 |
| holon–antiholon | 0.25 | 2.3929572325 | 9.1469949420 |
| holon–antiholon | 0.5 | 4.0090497566 | 8.7248705500 |
| holon–antiholon | 1 | 5.2867270220 | 6.9992246552 |

The first lower edge has constituent momenta (pi/2,-pi/2), an interior
minimum. At P=pi/4 its minimizing spinon momentum is approximately
2.4552859, so equal sharing is not a general rule either. At P=pi/2 the
minimum instead reaches a zero-energy spinon endpoint and the holon gap.
The antiholon run agrees with the mixed-holon energies at the pi-shifted
momenta to better than 1e-9.

Native tests now compare the oracle at U=1, 4 and 16, check stationary
momentum sharing by independently differentiating constituent energy sums,
and verify momentum shifts, convention offsets, shared budgets and missing
outputs. The default-tolerance holon–antiholon edges at P=0 are compared to
the independent high-precision charge-gap reference in fp64, long-double and
fp128. Frontend regressions cover channel selection, convention and reference
mapping, budget forwarding, streaming JSON/CSV/TSV, native precisions and
validation before output-file creation.
