# Charge-containing Hubbard continua: implementation plan

**Status: source conventions and independent reference oracle; no native
charge-continuum solver or CLI selection yet.** The existing
[two-spinon tool](hubbard-continuum.md) remains restricted to two spinons.

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

Before publishing a native solver, verify these references across couplings
and meshes, add native-precision stationary-point checks, and test the
momentum shifts, convention offsets, shared budgets and missing-output rules.
