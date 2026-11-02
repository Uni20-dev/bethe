# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |
| `bethe::experimental::xxz::VertexEvaluator<Real>` | Infinite massive XXZ two/four-spinon complex vertex amplitudes, with independent product/residue controls |
| `bethe::experimental::xxz::ThermodynamicVertexSpectrum<Real>` | Fixed-order full-zone moments and Gaussian frequency cuts, including coherent vacuum sectors |

See [finite-q XXZ spectra](../docs/experimental-xxz-vertex.md) for the reference
scripts, normalized measures and explicit convergence limitations.

Finite complex-root solvers and multiparticle form-factor consumers are not
provided by the singlet foundation. Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.
