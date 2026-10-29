# Developing Bethe

[Embed the C++ library](library.md) · [Testing and CI](testing.md) ·
[Benchmarks](benchmarks.md) · [Model roadmap](model-roadmap.md)

<a id="source-layout-and-development"></a>

## Examples and source layout

Set `BETHE_BUILD_EXAMPLES=ON` to build library composition examples in the
`examples/` subdirectory. The first example prepares one Lieb–Liniger background
asynchronously and submits both elementary branches before observing it; see
[shared backgrounds](async-execution.md#shared-backgrounds). A second example
evaluates a bounded [ground-state parameter sweep](async-examples.md#ground-state-parameter-sweeps),
using `mpreal` when Uni20 MPFR scalars are configured. Examples are checked
in CI and are not installed with the command-line tools.
The [ASEP example](asep.md#reusing-a-seed-across-rate-queries) compares independent
rate solves with queries sharing one TASEP seed, using direct or Async execution.
The [TASEP/XYZ sweep example](async-examples.md#independent-tasep-and-finite-xyz-solves)
shows independent complex-root solves with owned requests and bounded admission.
The [nested sweep example](async-examples.md#owning-nested-model-requests) demonstrates
independent TB sizes, SU(N) couplings and t–J fillings with complete diagnostics.
The [ladder example](async-examples.md#prepared-ladder-field-and-rung-sweeps) compares
repeated scans with shared root preparation and dependent field/rung queries.

### Where code belongs

| Location | Responsibility |
| --- | --- |
| `include/bethe/` | Scalar-templated numerical API, independent of CLI parsing and formatting |
| `include/bethe/detail/` | Shared numerical implementation helpers |
| `apps/` | Command-line frontends, model-specific arguments, formatting and report metadata |
| `tests/` | Numerical regressions, independent reference checks and frontend tests |
| `data/citations.json` | Literature metadata used to generate the C++ registry and bibliography |

Formatting stays in `apps/`, separate from the numerical API and any future
Python bindings. All executables link the private `bethe_cli` helper target;
numerical headers and data-output unit tests remain parser-independent.

The repository's `.clang-format` is copied from Uni20; use `clang-format -i`
on changed C++ files to apply the shared style.

### Linear solves and failure policy

We link `uni20_core` for scalar facilities, `uni20_common` for `half_int` and
presentation facilities, and `uni20_linalg` for square and least-squares solves.
Frontends additionally use `uni20_cli` for option parsing.

Square Newton and dressing-equation solves use Uni20's recoverable
`solve_inplace_with_info`, with LAPACK/MPLAPACK for compatible column-major
workspaces and native CPU fallback where needed. The row-major vector
adapter in `detail/newton.hpp` rearranges its owned coefficient copy in place
and exposes non-owning column-major views; the caller's coefficients are
preserved without an additional matrix allocation.

That adapter retains Bethe's explicit `64 * epsilon<Real>()` relative-pivot
threshold. Models already using Uni20 matrix workspaces retain the default
zero threshold and their own residual/condition checks. Neither policy is a
condition-number guarantee. Singular, rejected-pivot and nonfinite failures
return to the model's stalled/ill-conditioned or continuation-rejection path;
failed workspaces never become accepted iterates. Shape errors remain distinct
from numerical failures. No process-global error policy is changed.

Rectangular Newton corrections and continuation tangents use Uni20's
`least_squares_with_info` through read-only row-major views in `detail/newton.hpp`.
Normal dispatch selects LAPACK/MPLAPACK QR or Uni20's CPU Givens implementation;
native long double remains supported without MPLAPACK. Bethe retains its
`64 * epsilon` relative-diagonal rejection threshold, which is neither a
numerical-rank determination nor a condition estimate. No normal equations
are formed. A failed correction leaves the RHS unchanged; success replaces
it with the solution. An inconsistent least-squares fit can succeed: physical
residual checks and continuation acceptance still belong to the model.

The rectangular adapter also accepts an explicit finite `uni20::Precision`
for runtime-precision scalars. At this low-level boundary Uni20 converts inputs
at the requested working precision. The Richardson and central-spin model APIs
are stricter: inputs and tolerance must already have the options' finite
precision; mixed, exact and unset runtime values require explicit conversion.
Their sums, predictor/corrector state, backtracking limits and output diagnostics
retain that precision, including zero-budget and analytic cases. Numerical rank
via SVD, condition estimates and
refinement remain coordinated with
[Uni20's numerical infrastructure plan](https://github.com/Uni20-dev/uni20/blob/717c798c047c9e844535d9b21a09c152b49ea9af/docs/linalg/precision_numerics_plan.md).

### Shared iteration and continuation

`detail/newton_backtracking.hpp` shares trial construction, normalization
and damping. Models retain their pivot policy,
admissible domain, residual acceptance and iteration-counter meaning.
The string solver still requires strict decrease even below tolerance.

`detail/eigenvalue_continuation.hpp` shares the Richardson/central-spin adaptive
predictor-corrector controller and trust-box correction, including attempted-step
budgets and rejected-stage handling. Tangents, trust-radius selection, physical
energy bounds and final diagnostics remain model-specific. The tests exercise
these contracts in each enabled precision; continuation policy remains in
Bethe while square factorization and numerical solve diagnostics belong to Uni20.

### Frontend helpers

Use the shared lifecycle and output helpers; keep scientific arguments,
convergence criteria and missing-value policy in the model frontend.
Paths below are relative to `apps/`.

| Header | Responsibility |
| --- | --- |
| `cli-common.hpp`, `precision.hpp` | Common CLI utilities, checked precision requests and dispatch |
| `program-options.hpp` | Parse/information/error lifecycle, exact option adapters, Bethe identity and citations |
| `report-common.hpp`, `excitation-report.hpp` | Shared report rendering |
| `data-output-options.hpp` | Export flags for typed-table frontends |
| `run-metadata.hpp` | Uni20 run metadata projected into human overviews and legacy export keys |
| `result-output.hpp` | A frozen numerical summary shared across batch result tables |

Each frontend declares its own `TableDescriptor` list when registering output
options. The shared `table-descriptors.hpp` defines only the descriptor type and
generic validation: it has no executable-name registry or model-specific tables.
Declarations supply help, parser choices, screen flags, availability notes and
explicit primary-table eligibility. Result construction supplies the actual
available tables, which must match those declarations before exports are opened.

Keep display flags separate from production requirements. Use
`output.needs("--roots", roots)` when deciding whether root data is needed;
the unchanged `roots` flag controls display, while a named export can request
the same producer independently. Shared producers may construct sibling tables,
but only requested destinations receive them. Never manufacture an empty table
to satisfy an unavailable name. `Output.table_contracts` checks frontend help,
selection errors and auxiliary-file/screen separation.

The biquadratic frontend has private headers in `apps/biquadratic/`:

- `options.hpp`: option declarations and validation.
- `report.hpp`: common metadata, spin-content and real-root table writers.
- `real.hpp`, `analytic.hpp`, `qsystem.hpp`, `singlet.hpp`, `clusters.hpp`:
  individual calculation families.

These are frontend implementation details, not public numerical APIs.

### Citation maintenance

Literature metadata is centralized in `data/citations.json`. See
[maintaining citations](citations.md) to regenerate the C++ registry and the
bibliography in `CITATIONS.md`. Generated files are checked in, so Python is
needed only for regeneration and optional maintainer tests, not normal builds.

### Documentation equations

Use GitHub's backtick-protected inline math and fenced `math` display blocks.
In all math, use `\lt` and `\gt` instead of literal angle brackets: these can
survive Markdown processing but fail in GitHub's browser renderer. Use `\mathrm`
for named functions such as atan2; GitHub rejects `\operatorname` even though
MathJax itself supports it. End multiline display rows with `\\{}`: GitHub's Markdown
renderer adds an unwanted backslash to a bare `\\` at the end of a line.
The empty group preserves the TeX row break without relying on trailing spaces.

Check the conventions locally, or verify that every expression survives
GitHub's Markdown processing unchanged (requires authenticated `gh` and network
access; renders only, without publishing):

```sh
python3 scripts/check_markdown_math.py
python3 scripts/check_markdown_math.py --github
```

This checks the Markdown-to-TeX boundary and known browser restrictions, not
the full MathJax grammar or mathematical content. The rendering API does not
run GitHub's browser math component: confirm unfamiliar notation in a browser
too. Testing raw TeX alone does not catch Markdown escaping or GitHub's restrictions.

## Design and investigation records

- [Model roadmap](model-roadmap.md): candidate families, references and proposed extensions.
- [Completed MPS excitation sequence](mps-excitation-goal.md): the original milestones and their outcomes.
- [Hubbard charge-continuum design](hubbard-charge-continuum-design.md): the implemented construction and validation rationale.
- [Dependency sanitizer investigation](async-sanitizer-follow-up.md): dated findings, reproductions and upstream resolutions.
- [Precision/execution validation records](validation-records.md) and
  [performance measurements](benchmarks.md): checkpoint evidence with original configurations.

Current usage belongs in model/API guides. Historical records keep the date and
conditions of the original investigation, with a link to current instructions.
