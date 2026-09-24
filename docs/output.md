# Result tables and file exports

All 19 frontends share these options. Their numerical libraries remain
independent of command-line parsing and output sinks.

One calculation produces typed Uni20 data tables. Screen output and file
exports are independent views of those same values:

```sh
# See the usual final report and save two machine-readable files.
build/bethe-hubbard-dispersion --u 4 --points 101 \
  --csv half-filled.csv --json half-filled.json

# Doped curves, without a screen report or retained row history.
build/bethe-hubbard-dispersion --u 4 --density 0.75 --reference fermi \
  --points 1001 --quiet --no-retain --tsv doped.tsv

# Display points as they are computed and also save the full-precision table.
build/bethe-hubbard-dispersion --u 4 --stream --csv live.csv
```

`--csv FILE`, `--tsv FILE` and `--json FILE` are repeatable: several destinations
can be attached to one calculation. `--format auto|pretty|plain|csv|tsv|json`
controls **stdout only**. Use `--format csv > result.csv` for shell redirection,
not `--csv -`. `--quiet` suppresses stdout, but not exports or stderr warnings.

Existing files are refused unless `--force` explicitly allows replacement.
Duplicate destinations, including symlink/hardlink aliases and aliases of
redirected stdout, are rejected. Files must be regular files. Do not point an
export at an input, source file or other valuable data with `--force`.

## Named tables

All except Hubbard dispersion use a named-table document: Haldane–Shastry has `levels`;
Sutherland has `states` and, with `--pseudomomenta`, `pseudomomenta`.
The zero-based `state_id` links auxiliary rows to state records. Sutherland's
old space-separated pseudomomenta cell is replaced by typed `state_id,index,label,k`
rows; this avoids parsing numeric lists out of strings.

The remaining models use these tables:

| Frontend | Primary table | Optional auxiliary tables |
| --- | --- | --- |
| `bethe-su3-pbc` | `states` | `first_roots`, `second_roots` with `--roots` |
| `bethe-tb-pbc` | `states` | `strings`, `roots` with `--roots`; complex roots have separate real/imaginary columns |
| `bethe-richardson` | `states` | `variables` with `--variables`; blocked levels have null eigenvalue variables |
| `bethe-central-spin` | `states` | `variables` with `--variables`; null coupling identifies the central spin |
| `bethe-gaudin-yang-pbc` | `states` | `charge_roots`, `spin_roots`, or `free_up`, `free_down` with `--roots` |
| `bethe-tj-pbc` | `states` | `first_roots`, `second_roots`, or `free_modes` with `--roots` |
| `bethe-sun-fermions-pbc` | `states` | always `components`; `roots` or `free_modes` with `--roots` |
| `bethe-ladder-pbc` | `states` | always `representations`; selected-state `roots` with `--roots` |
| `bethe-hubbard-pbc`, `bethe-hubbard-obc` | `states` | `charge_roots`, `spin_roots`, or `free_modes` with `--roots` |
| `bethe-xxx-pbc`, `bethe-xxx-obc` | `states` | `roots` with `--roots`; periodic `--spinons` also has `spinons` |
| `bethe-xxz-pbc`, `bethe-xxz-obc` | `states` | `roots` with `--roots`; open ground-state modes also have `boundary_roots` |
| `bethe-lieb-liniger-pbc` | `states` | `roots` with `--roots` |
| `bethe-biquadratic-obc` | `states` | real modes: `quantum_numbers`; Q-system: `reference`, `q_coefficients`; two-string singlet / ferro bound pairs: `reference`, `string`; ferro analytic modes: `reference`; numerical `roots` with `--roots` |

Excitation scans additionally write `reference` (the ground state used for gaps)
and, if a candidate fails, `failed` (the first unranked estimate). State IDs are
unique across these tables: ranked states come first, then the reference, then
the failed candidate. Root rows may refer to any of them. XXX/XXZ scans also
write a `quantum_numbers` table, preserving exact labels even without `--roots`.
`state_id` is zero-based, replacing the old spin-chain one-based display level.

Failed candidates never enter ranked `states`; a failed reference leaves their
`gap` values null. Sector scans still retain unconverged estimates with explicit
status. Biquadratic `through_lines` is a TL module label, **not** physical spin;
its nullable integer `multiplicity` is null on overflow, not zero.
Its Q-system searches report count-matched numerical completeness separately
from document completion. Incomplete `--q-spectrum` rows are sorted discoveries,
not guaranteed lowest levels; a failed `--q-seed` row is an unverified estimate
with no gap. Q-system root columns use `x=cosh(2u)`, not the real solver's alpha.
The targeted `--singlet-excitation` mode retains `string.log_deviation=L=-log(d)`;
its rounded complex roots cannot always resolve d. The optional `deviation`
column is null on underflow. Failed targeted estimates remain visible with
`converged=false` and null gap, never as ranked eigenvalues.

XXZ `roots.lambda` is populated only for negative anisotropy; do not reconstruct
it from rounded scaled rapidities. Massive open XXZ boundary coordinates live
in `boundary_roots`, never in the bulk `roots` table. Its `inverse_square`,
`log_distance`, and `kind` retain the real/imaginary/infinite-root distinction.
`root_delta` records the reached continuation stage; energy and residual still
refer to the requested Delta. An empty boundary table means no distinguished
boundary root, not a failed solve.

SU(n) component indices preserve input order; `nesting_rank` is null for an
empty component. Root level zero contains charge momenta; higher levels contain
spin rapidities. Energies and momenta are null if no coupling stage converged,
and otherwise belong to `reached_c`. Gaudin–Yang retains unconverged root
estimates and reports `root_c`; its residuals are evaluated at the requested c.

Ladder `states` marks the selected row. Only that row has roots when scanning
sectors; `representations` distinguishes physical populations from SU(4)
highest-weight rows. An incomplete scan energy is a candidate upper bound,
not a certified minimum. Hubbard exports identify auxiliary-sector roots after
symmetry mappings; free modes have no Bethe labels. The open-chain state table
deliberately has no momentum column: standing-wave k is not lattice momentum.

Richardson energies belong to `reached_g`, not necessarily `requested_g`.
Central-spin energies similarly belong to `reached_field`; both energy and
reached field are null if no finite-field stage was reached. A null central-spin
variable denotes an analytic state, not a missing occupation. These are
eigenvalue variables, not occupations or Bethe rapidities.

The final human report keeps the model's compact overview once, followed by
its tables. Each exported table independently carries full provenance and
physical metadata so a selected auxiliary CSV file remains interpretable.

Human output shows all requested tables. JSON stores them in a single object:
`{"tables":{"states":{...},"pseudomomenta":{...}},"status":"complete"}`.
Each table has the Uni20 schema described below. The document status records
output completion, **not** numerical convergence or spectral completeness;
inspect the individual table summaries and the command's exit status.
Hubbard dispersion retains its existing single-table JSON shape.

The human-readable layout has changed from the old hand-written whitespace
rows and interleaved root blocks. For scripts, use CSV/TSV or JSON and select
columns by their schema identifiers, rather than parsing the screen report.

CSV/TSV always contain exactly one rectangular table. By default `--csv FILE`,
`--tsv FILE` and delimited stdout select the primary table. `--table NAME` changes
that selection. Repeatable `--csv-table NAME=FILE` / `--tsv-table NAME=FILE` export
additional tables without mixing schemas:

```sh
build/bethe-sutherland-pbc 3 --length 4 --lambda 2 --levels all --window 2 \
  --pseudomomenta --csv states.csv --tsv-table pseudomomenta=roots.tsv --json all.json
build/bethe-sutherland-pbc 3 --length 4 --lambda 2 \
  --pseudomomenta --table pseudomomenta --format csv
```

Selecting an unavailable table is an error before any export is opened; e.g.
`--table pseudomomenta` requires `--pseudomomenta`. The `--table` selector affects
only CSV/TSV, not the screen or JSON table collection.

## Metadata and reproducibility

CSV/TSV start with MPToolkit-style `# key: value` comments: program and version,
build-time Bethe and Uni20 revisions, UTC timestamp, shell-quoted command line,
Hamiltonian/energy conventions, physical parameters, scalar precision and
effective solver controls. Doped calculations also include the solved
background. Revisions describe the actual build checkouts, including a local
Uni20 override; `-dirty` means tracked modifications at build time. Untracked
files are not fingerprinted. A source without Git information reports
`unavailable` (or `unknown` through Uni20's provenance provider) rather than
claiming the pinned revision was used.

The rectangular header and data rows follow. Status and timing appear in
**trailing comments**; Hubbard dispersion also reports its accepted row count.
Readers should ignore `#` lines throughout the file, not only at its beginning.

All frontends use Uni20's typed run context for scalar metadata, provenance and
timing. Real values retain their selected precision, counts remain integers,
and optional values remain distinguishable from empty text until export. The
human overview and exported metadata use the same native values; metadata is
no longer copied out of a formatted report. Existing export keys and JSON
structures are retained. Half-integer metadata is now consistently decimal
(`0.5`), even where the human overview keeps a fraction (`1/2`). Model-specific
missing-value explanations, such as overflow or an unreached continuation
stage, are preserved. Metadata also records the compiler, build type and
platform. Its UTC timestamp marks the start of the selected-precision
calculation, before solving.

Every summary includes these timings:

- `CPU time`: numerical solve/enumeration CPU time, excluding report rendering
  and export work, with six fractional digits and an `s` suffix.
- `Run CPU seconds`: process CPU from entry into the selected-precision
  calculation until its numerical summary is frozen, including validation and
  setup within that function.
- `Elapsed seconds`: monotonic wall time over that same interval.

The new seconds fields contain round-trip decimal text without unit suffixes;
unavailable timing is `unavailable`. CLI parsing is outside this interval.
Batch calculations freeze one summary **before emitting their result tables**,
then copy it to every table: these are run timings, not per-table timings, and
do not include table rendering or export. Hubbard dispersion instead solves
and streams points in one loop, so its run CPU and elapsed time include the
interleaved output. Neither includes the final summary's rendering or flush.

`Outcome` is set explicitly by the model: `success` for converged results or
exact spectral rules, `partial` for incomplete scans, failed references,
unreached targets or enumeration-budget refusals. A completed JSON document
can therefore have `"status":"complete"` while its table summary says
`"Outcome":"partial"`; successful serialization does not establish scientific
convergence. Hubbard dispersion can also record `failed` if it aborts before
freezing its summary. Errors before output begins may produce only a diagnostic.
A later write/flush failure still causes exit 1 without rewriting a frozen
numerical outcome; named JSON documents report an independent transport abort.

This migration retains the existing CLI, output coordinator and final-report
layout. Staged human preambles, shared output sessions, configuration sources,
and a multi-table JSON format change are separate, not enabled here.

`--no-preamble` removes both initial and trailing CSV/TSV comments, leaving
strict rectangular data for readers that do not support comments. It does not
remove metadata from JSON or the screen report. Comments escape backslashes
and control characters (`\\`, `\n`, `\r`, `\t`, `\xHH`) so an unusual filename
cannot introduce extra data lines. Decode these escapes before interpreting
the shell-quoted `Command` field. JSON handles those characters with ordinary
JSON escaping.

## Numeric encoding

Machine output uses round-trip precision for the selected scalar, without
narrowing through `double` or wrapping numeric tokens. Half-integer quantities
are written as mathematical decimal values: spin one-half is `0.5`, not `1/2`
or the doubled integer `1`.

CSV/TSV use empty cells for unavailable optional values. JSON is a single
object containing `title`, `metadata`, `first_row`, `columns`, `rows`, and
`summary`. Rows are arrays in column order. Each column describes its identifier,
type, nullability, units and, for real numbers, radix and precision in bits.

- Real values are **decimal strings** in JSON, including fp64. This avoids
  silently rounding long-double/fp128 when a JSON reader uses binary64.
- Half-integers are JSON numbers such as `0.5`, with type `half_int` in the
  schema. Small signed integers such as particle changes are numbers too.
- Wide integer counters are decimal strings, as recorded by their column's
  `encoding`. They are not restricted to JavaScript's exact integer range.
- Missing values are `null`. Nonfinite real values are strings such as `"inf"`
  and `"-inf"`; infinite spinon endpoint rapidities are legitimate, not failed
  energies. Check the row's `status`, not just its `parameter`.

For example, load a JSON export without discarding the precision of an energy:

```python
import json
from decimal import Decimal

with open("half-filled.json") as source:
    table = json.load(source, parse_float=Decimal)
columns = [column["id"] for column in table["columns"]]
for values in table["rows"]:
    row = dict(zip(columns, values))
    if row["status"] == "converged":
        energy = Decimal(row["energy"])
```

## Retention, streaming and failures

Rows are retained in memory by default. Human stdout is a final report built
from this retained table. `--stream` instead displays human-readable rows as
they are delivered; narrow screens can use vertical records to keep numbers
intact. Files and machine stdout always receive rows incrementally, with normal
stream buffering and an explicit final flush. JSON becomes a complete document
only at finalization. The background solve precedes all output.

`--no-retain` discards accepted rows after delivery, useful for large grids.
It requires `--stream`, a machine stdout format, or `--quiet`; a final human
snapshot cannot be reconstructed without stored rows. This is the same typed
table API with a retention option, not a separate streaming data model.
Spectral scans still retain and sort their solver results before delivering rows;
`--no-retain` does not remove that solver storage or make sorted levels available early.
At the C++ level Uni20 also supports attaching sinks later and replaying
retained rows. The dispersion frontend attaches sinks before point solves;
finite-state and spectral frontends finish solving before delivering their tables.

Exit 0 means successful computation and output. Exit 2 means incomplete
numerical results: failed dispersion points have missing energies, finite
solves can retain labelled estimates, and excitation scans exclude failed
candidates from ranked levels. Inspect the model-specific status and reference
tables. Exit 1 means invalid arguments, a runtime error, or a required output
failure. Invalid scientific arguments are checked before opening exports.
For the exact-rule models, exit 2 instead means enumeration exceeded its budget;
no partial selection is presented as the lowest spectrum.

Output is **not transactional**. A write/open failure can leave partial or empty
files; abrupt termination can leave JSON incomplete. A handled failure stops
the calculation, reports the failed destination, and attempts to finalize
healthy sinks with `Status: aborted`. Accepted rows are never retried, so a
healthy sink does not receive duplicates. If failure occurs during final flush
or screen rendering, the numerical `converged` summary already written to
healthy files remains valid, but the command still exits 1. Always check the
exit status as well as the data summary.

## Implementation boundary

`apps/data-output.hpp` owns output configuration, streams, provenance and the
commented CSV/TSV adapters. `apps/data-output-options.hpp` declares the CLI flags
using Uni20/CLI11, without opening files during parsing. Uni20 owns typed columns,
retention/replay, numeric encoding and sinks. `apps/result-output.hpp` adapts a
model overview and several native tables to one output document. Model-specific
schemas and physical metadata stay in the frontends; solvers in `include/bethe/`
do not depend on output policy.
This division leaves the numerical API usable by future Python bindings.
