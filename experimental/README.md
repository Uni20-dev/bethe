# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |

Finite complex-root solvers and multiparticle form-factor consumers are not
provided by this first foundation. Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.

The [open XXX boundary guide](../docs/experimental-xxx-boundary-form-factors.md)
covers normalized local raising weights and a fresh standalone example.
