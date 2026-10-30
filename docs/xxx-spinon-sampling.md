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

`SpinonSampleStream<Dimension>` keys an ordinary Monte Carlo stream by master
seed, replicate id and a declared uniform bit count. It uses standardized
`seed_seq` and `mt19937_64` engines. A draw retains its binary integer words;
`SpinonUniforms::at_precision` constructs the exact midpoint `(integer+1/2)/2^b`
at the requested precision. At least b+1 significand bits are required. A retry
reuses those words and never redraws or silently lowers the grid resolution.
This stream is not the historical campaign's scrambled quasi-Monte Carlo rule.

`SpinonReplicateStatistics<Real>` collects complete, equal-size replicate
means. It reduces them in replicate-id order and supplies the full covariance
of the mean, preserving correlations between frequencies or moment orders.
One replicate supplies no sampling-error estimate. Missing or failed replicates
leave the result incomplete/failed with no successful mean over the survivors.
Genuine support zeros count as observations; numerical failures do not.
Underflow of a nonzero variance is an explicit numerical failure, not a zero
sampling error. These statistical diagnostics do not include cutoff or numerical
integration bias, and the uniform-grid resolution is a separate input choice.

## Energy and momentum constraints

`SpinonEnergySurface<Real>` in `bethe/detail/xxx_spinon_surface.hpp` accepts
the spinon count (2, 4 or 6), cutoff, momentum in [0,2 pi), momentum-section
index and positive frequency. Its `solve` method takes n-2 free rapidities and
returns both labelled solutions for the remaining pair. The section fixes
the total spinon momentum to k+2 pi s; sum every admissible section.

The map's `log_coarea_jacobian` is the inverse velocity difference per root,
before proposal or free-rapidity Jacobians. Use either that unbalanced chart
measure with the complete proposal PDF or the balanced measure below, not
their product. The map has no physical form-factor normalization.

Resolved support exclusions return `inside=false`. Ambiguous boundaries,
unresolved cutoffs and numerical collisions of the solved roots require higher
precision or a separate limiting calculation. In particular, resolving the
momentum cutoff for R=1040 needs much more arithmetic precision than evaluating
the proposal's logarithmic PDF. The tests exercise this geometry at 1800 and
2048 bits; successful low-precision proposal evaluation is not a geometry check.

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

`test_xxx_spinon_surface.cpp` independently checks the inverse map against
forward kinematics and velocity differences at prescribed rapidities, including
all three six-spinon momentum sections. Its optional reference generator uses
720 decimal digits, sufficient to retain the deep-soft inputs tested here.
