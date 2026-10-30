# Six-spinon sampling geometry

The internal `bethe::heisenberg::detail::SpinonJointProposal<Real>` in
`bethe/detail/xxx_spinon_proposal.hpp` generates four free rapidities for a
fixed-energy, fixed-momentum six-spinon integration. It supplies proposal
densities and balanced coarea weights. Physical form factors and the solution
of the two conservation equations are separate components.

## Fresh draws

The constructor takes a rapidity cutoff R, a soft fraction strictly between
zero and one, and positive mixture weights that sum to one. Runtime-precision
types require an explicit finite precision. Five open uniform coordinates
select a component and its four continuous coordinates. No saved sample file,
Python process or global random-number state is needed.

The mixture contains a background, independent Cauchy rapidities, one soft pair,
or two soft pairs. A pair uses its centre B and separation d with unit Jacobian:

```math
\beta_i=B+d/2,\qquad \beta_j=B-d/2,\qquad
|B|<R,\quad |d|<2(R-|B|).
```

Both variables have normalized truncated Cauchy distributions. One-pair draws
choose uniformly from all six free-label pairs; two-pair draws choose among
the three matchings. The background provides full support. Each returned
`log_density` evaluates the **complete mixture**, including components that
could also have produced the draw. Its measure is product(d beta), not
product(dp). The parameters use the convention p=2 atan(exp(beta)).

Random-number generation, independent replicates and sample sharing across
frequencies belong to the caller. A spectral estimator must preserve those
replicate identities when reporting uncertainties and covariances.

## Balanced coarea charts

For each possible solved pair (i,j), let q_ij be the complete proposal density
of the four remaining rapidities. At J=1 its pushforward intensity is

```math
h_{ij}=q_{ij}\prod_{a\ne i,j}\cosh\beta_a\,
\frac{\pi}{2}\frac{|\sinh(\beta_i-\beta_j)|}
{\cosh\beta_i\cosh\beta_j}.
```

This convention emits **both labelled solutions** of the conservation equations.
For a permutation-invariant physical density and integration domain, symmetry
allows replacing the original chart's inverse intensity by the inverse mean
of all fifteen intensities:

```math
J_{\rm balanced}=\frac{15}{\sum_{i<j}h_{ij}}.
```

`balanced_log_jacobian` returns its logarithm per labelled solution. The caller
sums both solutions and all allowed momentum sections. The physical density
retains its single state-counting factorial; this balancing factor does not
change that normalization. An ordered-domain or hardest-pair restriction is
a different estimator and must not be stacked on this full-domain balance.

Logarithmic assembly keeps the separate cosh and proposal factors from
overflowing. It does not establish that a subsequent momentum map or energy
constraint is resolved at the same precision. Unresolved draws fail explicitly;
callers must recover them or report an incomplete estimate rather than discard
them. A vanishing individual chart is permitted; an entirely degenerate chart
orbit requires a separate limiting calculation.

## Validation

`test_xxx_spinon_proposal.cpp` checks independent high-precision density,
coarea and draw references, all pair choices, permutation symmetry and the
literal 720-label orbit measure. It also checks invalid inputs, inverse CDFs,
deep-soft configurations and explicit precision growth beyond binary128.
The optional `tests/reference_xxx_spinon_proposal.py` regenerates the small
fixtures with mpmath using direct positive formulas and a bisection inverse.
