# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::scan_four_spinons_async` | Bounded finite-four requests with deterministic records, explicit failures, partial spectra and table exports |
| `bethe::experimental::xxx::FourSpinonState<Real>` | Selected finite narrow two-string triplets, explicit deviations and stabilized raising weights |
| `bethe::experimental::xxx::complex_roots` | Regular complex-root equations, Gaudin matrices and unreduced raising form factors |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |

Regular finite-root equations and unreduced form factors are described in
[the complex-root guide](../docs/experimental-xxx-complex-roots.md). They are
building blocks for the [selected finite-four solver](../docs/experimental-xxx-four-spinon.md).
Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.
