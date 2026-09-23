# Result tables and file exports

At present these options are implemented by `bethe-hubbard-dispersion`,
`bethe-haldane-shastry-pbc` and `bethe-sutherland-pbc`. Other
frontends retain their existing output options; they can migrate to the same
shared adapter without changing their numerical libraries.

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

The newer frontends use a named-table document: Haldane–Shastry has `levels`;
Sutherland has `states` and, with `--pseudomomenta`, `pseudomomenta`.
The zero-based `state_id` links auxiliary rows to the primary table. Sutherland's
old space-separated pseudomomenta cell is replaced by typed `state_id,index,label,k`
rows; this avoids parsing numeric lists out of strings.

Human output shows all requested tables. JSON stores them in a single object:
`{"tables":{"states":{...},"pseudomomenta":{...}},"status":"complete"}`.
Each table has the Uni20 schema described below. The document status records
output completion, **not** numerical convergence or spectral completeness;
inspect the individual table summaries and the command's exit status.
Hubbard dispersion retains its existing single-table JSON shape.

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
`unavailable` rather than claiming the pinned revision was used.

The rectangular header and data rows follow. Completion status, accepted row
count and CPU time appear in **trailing comments**, since they are not known
when streaming starts. Readers should ignore `#` lines throughout the file,
not only at its beginning. CPU time measures the numerical background and
point solves, excluding rendering and export work; it is not wall-clock time.

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
retained rows; the current CLI attaches all requested sinks before point solves.

Exit 0 means successful computation and output. Exit 2 means some points did
not converge: their rows still appear, with missing energies and diagnostic
statuses. Exit 1 means invalid arguments, a runtime error, or a required output
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
retention/replay, numeric encoding and sinks. Hubbard-specific column schemas and
physical metadata stay in its frontend; solvers in `include/bethe/` do not depend on output policy.
This division leaves the numerical API usable by future Python bindings.
