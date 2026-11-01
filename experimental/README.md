# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::complex_roots` | Regular complex-root equations, Gaudin matrices and unreduced raising form factors |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |
| `bethe::experimental::xxx::six_spinon` | Selected two-pair and three-string triplets with resolved deviations and explicit failed results |

Regular finite-root equations and unreduced form factors are described in
[the complex-root guide](../docs/experimental-xxx-complex-roots.md). They are
building blocks used by the [selected six-spinon solver](../docs/experimental-xxx-six-spinon.md).
Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.
