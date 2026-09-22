# Maintaining literature citations

[Back to the overview](../README.md)

The source of truth is [data/citations.json](../data/citations.json). It contains
stable reference IDs, bibliographic metadata, links, and each tool's selection
of references with a short explanation of their relevance. A reference may
also be library-only: the thermodynamic XXZ spinon notes, for example, belong
in the bibliography without appearing as a mode of `bethe-xxz-pbc`.
References used only by the [model catalogue](models.md) follow the same rule:
add bibliographic records, but no tool selection for an unimplemented model.
Presence in the registry or bibliography does not imply solver support.

The registry deliberately separates three concerns:

- **Reference metadata:** authors, title, publication, year, and source links,
  recorded once under a stable ID.
- **Selection and context:** which methods or tool modes use each reference.
  Citing classification background does not imply the corresponding strings,
  descendants, or other algorithms are implemented.
- **Presentation:** plain help text, a Markdown bibliography, or a future
  BibTeX export. Numerical routines do not print citations or load files.

## Editing and checking

1. Edit `data/citations.json`; reuse existing IDs rather than duplicating a work.
2. Update the hand-written provenance in [CITATIONS.md](../CITATIONS.md) if the
   equations, normalization, or applicability changed. Link to the stable ID's
   heading, e.g. `[lieb-wu-2003](#lieb-wu-2003)` within that document.
3. Regenerate and check:

   ```sh
   python3 scripts/generate_citations.py
   python3 scripts/generate_citations.py --check
   ```

The generator updates `include/bethe/citations.hpp` and only the marked
bibliography section of `CITATIONS.md`; it preserves prose outside the markers.
Both generated outputs are checked in, so normal builds and executables need
neither Python, JSON parsing, internet access, nor a source-tree data file.
Regeneration uses only the Python standard library. `--check` is read-only
and fails if either generated output is stale. With Python available, CTest
runs the freshness check and generator regression tests; otherwise CMake
reports that those maintainer checks are skipped. CLI citation tests still run.

Each front end calls the shared help renderer with a typed tool ID. Citations
appear with `--help` and in no-argument usage, on the same stream as the usage
text. Normal calculations keep their existing output and exit-status contract.
Help lists relevant references for the tool's supported modes, not an assertion
that every listed method was used in a particular run.

The generated data is also available to C++ callers independently of the CLI:

```cpp
#include <bethe/citations.hpp>

auto uses = bethe::citations::for_tool(bethe::citations::Tool::hubbard_pbc);
for (auto const& use : uses) {
    // use.reference has metadata and links; use.context explains applicability.
}
auto const* ref = bethe::citations::find("lieb-wu-2003"); // nullptr if unknown
```

## A possible shared facility with Uni20

This is initially a Bethe-local registry, not a new Uni20 dependency or an
upstream API commitment. The reusable part is the small metadata/ID/selection
interface and format-specific rendering, not the Bethe tool enumeration.

For Uni20, a natural extension is to expose reference IDs from algorithms and
backends, then let an application collect and deduplicate the references for
the choices actually used. Keep such a collection explicit and scoped to a
calculation; avoid global registration and unsolicited output from numerical
kernels. A compile-time list of available backends is not evidence that all
of them participated in a calculation.

A shared facility could support Markdown, plain text, and BibTeX output from
the same metadata, with project-specific provenance kept alongside the code.
If two projects share a work, use an agreed stable ID (and DOI where available)
to merge it. This keeps citation policy with the consuming application while
allowing Uni20 to supply authoritative citations for its own methods.
