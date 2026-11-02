# Experimental numerical components

Experimental components live under `experimental/` and have explicitly limited
API or scientific coverage. They use the supported library and Uni20; supported
headers do not depend on them. They are disabled by default and are not
installed with Bethe's command-line programs.

## Build and use

From your chosen build directory, enable the additional include target:

```sh
cmake /path/to/bethe -DBETHE_BUILD_EXPERIMENTAL=ON -DBETHE_BUILD_TESTS=ON
cmake --build . --target bethe_experimental_singlet_tests
ctest -L experimental --output-on-failure
```

For a CMake consumer, set the option before adding Bethe and link
`bethe::experimental`. Linking `bethe::bethe` alone does not expose the
experimental headers. Library-only builds can leave applications and tests
disabled; the experimental library target itself does not require GoogleTest
or CLI parsing. Its headers currently use namespace `bethe::experimental`.

The CI presets opt in and run bounded experimental tests alongside the existing
native, wide-precision and sanitizer coverage. Ordinary source builds retain
the default OFF setting. General test dependencies and Uni20's optional
sanitizer policy are shared between the two test directories.
CI also removes the explicit opt-in temporarily and checks that default-OFF
removes only experimental tests, retaining every supported test, then restores
the original configuration.

## Prepared finite-XXX singlet

```cpp
#include <bethe/experimental/xxx_singlet.hpp>

auto ground = bethe::heisenberg::ground_state<double>(8);
bethe::experimental::xxx::PreparedSinglet<double> prepared(8, ground);
if (!prepared.ready()) {
    // Retain the failed request and inspect prepared.status()/ground_residual().
    throw std::runtime_error("singlet preparation failed");
}
auto energy = *prepared.energy();
auto standard_roots = prepared.aba_rapidities();
auto zero = *prepared.working_zero();
```

This snapshot is for the even periodic XXX chain with J=1 and spin operators
sigma/2. Its input is the existing real-root `RealState`, whose coordinate is
z=2 lambda. It exposes standard ABA lambda. The energy is recomputed from the
validated roots rather than trusting the caller's saved energy.

Preparation reuses `PreparedRaisingFormFactors`: the ground Gaudin determinant
and ground-only products are calculated once. The read-only
`ground_log_prefactor()` is

```math
\log N-\log|\det G_0|+\sum_a\log(\lambda_a^2+1/4)
-\sum_{a\gt b}\log\left(1+(\lambda_a-\lambda_b)^2\right).
```

It is the ground-only contribution to a Fourier raising weight, not a complete
norm or form factor. Excited-state determinants and products are still needed.
The four-site analytic singlet provides an independent check: E=-2,
det(G_0)=108 and this prefactor is -log(324). Its pi-momentum triplet has
Fourier raising weight 4/3. See [finite-chain conventions](xxx-structure-factor.md).

The same snapshot can evaluate existing real-root triplets through
`real_raising_form_factor(excited)`. It does not yet solve complex-root states
or supply their form factors. The general numerical formulas and physical
normalization must be validated for each added family.

Use `working_zero()` to create scratch with the prepared precision. MPFR inputs
must already have one explicit finite precision; mixed, exact or unset root
precision is rejected. The context owns its snapshot, so later input mutation
or destruction cannot invalidate it. Const evaluations use local scratch.
Failed and moved-from contexts expose no usable energy, roots, prefactor,
pivot diagnostic or working exemplar. Failed preparations retain a measured
residual when available; moved-from contexts report an absent residual and
unavailable status.
The pivot value is a rejection indicator, not a condition number or error bound.

This interface intentionally has no disk-cache reader. Historical serialized
factors with a different determinant method or precision contract must not be
silently substituted for validated preparation.

## Independent QTM references

The [isotropic QTM reference package](experimental-xxx-qtm.md) compares JMP
quadrature with an analytic PRL hypergeometric limit, then checks complete
complex partition sums against the public spinon vertex controls. These are
optional Python source-checkout tools, separate from the C++ build option.
They need mpmath, but no private data or C++ build. Detailed limiting algebra
is in the [derivation](experimental-xxx-qtm-derivation.md).

## Scope and promotion

Each experimental feature must document supported inputs, numerical controls,
precision coverage and failure behavior. A limited but validated family can
later become a supported API; this requires independent physical checks,
precision/conditioning evidence, documentation and normal packaging. Supported
code cannot acquire a hidden dependency on experimental headers during that
promotion. Known wrong results are not acceptable experimental output.
