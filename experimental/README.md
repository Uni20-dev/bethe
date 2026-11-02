# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::scan_four_spinons_async` | Bounded finite-four requests with deterministic records, explicit failures, partial spectra and table exports |
| `bethe::experimental::xxx::FourSpinonState<Real>` | Selected finite narrow two-string triplets, explicit deviations and stabilized raising weights |
| `bethe::experimental::xxx::complex_roots` | Regular complex-root equations, Gaudin matrices and unreduced raising form factors |
| `reference/qtm/` | Optional isotropic QTM quadrature/hypergeometric and vertex correspondence controls; [guide](../docs/experimental-xxx-qtm.md) |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |
| `bethe::experimental::xxz::VertexEvaluator<Real>` | Infinite massive XXZ two/four-spinon complex vertex amplitudes, with independent product/residue controls |
| `bethe::experimental::xxz::ThermodynamicVertexSpectrum<Real>` | Fixed-order full-zone moments and Gaussian frequency cuts, including coherent vacuum sectors |
| [`bethe::experimental::xxx::SpinonRecoveryStore<Low, High>`](../docs/experimental-spinon-recovery-store.md) | Bounded POSIX restart storage for two-, four- and six-spinon density controls |

See [finite-q XXZ spectra](../docs/experimental-xxz-vertex.md) for the reference
scripts, normalized measures and explicit convergence limitations.

Regular finite-root equations and unreduced form factors are described in
[the complex-root guide](../docs/experimental-xxx-complex-roots.md). They are
building blocks for the [selected finite-four solver](../docs/experimental-xxx-four-spinon.md).
Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.

The [open XXX boundary guide](../docs/experimental-xxx-boundary-form-factors.md)
covers normalized local raising weights and a fresh standalone example.

The [boundary-dynamics guide](../docs/experimental-xxx-boundary-dynamics.md)
adds local/end-to-end time sums, an owned Async scan and the ground-only
boundary density matrix.

The [boundary-magnetization guide](../docs/experimental-xxx-boundary-magnetization.md)
adds signed same-sector endpoint elements and magnetization after a boundary spin flip.
