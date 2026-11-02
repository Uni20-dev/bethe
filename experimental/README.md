# Experimental Bethe components

These opt-in components have explicit numerical and API limitations. They are
tested, but are not installed with the supported command-line programs.
See [the build and API guide](../docs/experimental.md).

| Component | Current scope |
| --- | --- |
| `bethe::experimental::xxx::PreparedSinglet<Real>` | Validated even-ring XXX singlet preparation, reusable ground factors and energy, J=1 |
| [`bethe::experimental::xxx::SpinonRecoveryStore<Low, High>`](../docs/experimental-spinon-recovery-store.md) | Bounded POSIX restart storage for two-, four- and six-spinon density controls |

Finite complex-root solvers and multiparticle form-factor consumers are not
provided by this first foundation. Experimental status does not permit silent
failed-state deletion, altered physical normalization or an undocumented
ideal-string fallback.
