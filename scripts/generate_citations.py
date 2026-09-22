#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Generate checked-in citation data and bibliography; stdlib only, no network.

Normal C++ builds consume the generated header and do not require Python.
The hand-written prose outside CITATIONS.md's markers is never regenerated.
"""

import argparse
import json
from pathlib import Path
import re
import sys
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]
BEGIN = "<!-- BEGIN GENERATED BIBLIOGRAPHY -->"
END = "<!-- END GENERATED BIBLIOGRAPHY -->"


def fields(record, expected):
    if not isinstance(record, dict) or set(record) != set(expected.split()):
        raise ValueError(f"expected fields: {expected}; got {record!r}")


def line(value):
    if not isinstance(value, str) or not value.strip() or any(ord(c) < 32 for c in value):
        raise ValueError(f"expected nonempty single-line text: {value!r}")


def records(value):
    if not isinstance(value, list) or not value:
        raise ValueError("expected a nonempty list")
    return value


def validate(data):
    fields(data, "schema_version references tools")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise ValueError("unsupported citation schema_version")
    ids = set()
    for ref in records(data["references"]):
        fields(ref, "id authors title publication year links")
        for key in ("id", "authors", "title", "publication"):
            line(ref[key])
        if not re.fullmatch(r"[a-z][a-z0-9-]*", ref["id"]) or ref["id"] in ids:
            raise ValueError(f"invalid or duplicate reference ID: {ref['id']}")
        ids.add(ref["id"])
        if ref["year"] is not None and (type(ref["year"]) is not int or not 1 <= ref["year"] <= 9999):
            raise ValueError(f"invalid year: {ref['year']!r}")
        urls = set()
        for link in records(ref["links"]):
            fields(link, "label url")
            line(link["label"])
            line(link["url"])
            url = urlsplit(link["url"])
            if url.scheme != "https" or not url.netloc or any(c.isspace() for c in link["url"]):
                raise ValueError(f"expected an absolute HTTPS URL: {link['url']}")
            if link["url"] in urls:
                raise ValueError(f"duplicate link in {ref['id']}")
            urls.add(link["url"])
    tool_ids, executables = set(), set()
    for tool in records(data["tools"]):
        fields(tool, "id executable references")
        line(tool["id"])
        line(tool["executable"])
        if not re.fullmatch(r"[a-z][a-z0-9]*_[a-z][a-z0-9_]*", tool["id"]) or tool["id"] in tool_ids:
            raise ValueError(f"invalid or duplicate tool ID: {tool['id']}")
        if not re.fullmatch(r"bethe-[a-z0-9]+(?:-[a-z0-9]+)*-(pbc|obc)", tool["executable"]) or tool["executable"] in executables:
            raise ValueError(f"invalid or duplicate executable: {tool['executable']}")
        tool_ids.add(tool["id"])
        executables.add(tool["executable"])
        selected = set()
        for use in records(tool["references"]):
            fields(use, "id context")
            line(use["id"])
            line(use["context"])
            if use["id"] not in ids or use["id"] in selected:
                raise ValueError(f"unknown or duplicate reference {use['id']} in {tool['id']}")
            selected.add(use["id"])


def cpp(text):
    # JSON and C++ use compatible escapes for these validated single-line strings.
    # Retain UTF-8 author names rather than losing accents during generation.
    return json.dumps(text, ensure_ascii=False)


def header(data):
    result = ["""// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Generated from data/citations.json by scripts/generate_citations.py. Do not edit.
// clang-format off
#pragma once

#include <array>
#include <span>
#include <stdexcept>
#include <string_view>

namespace bethe::citations
{
struct Link { std::string_view label, url; };
struct Reference
{
  std::string_view id, authors, title, publication;
  int year; // 0 means unspecified, e.g. undated web notes.
  std::span<Link const> links;
};
struct Use { Reference const* reference; std::string_view context; };
"""]
    refs = data["references"]
    for i, ref in enumerate(refs):
        result.append(f"inline constexpr std::array<Link, {len(ref['links'])}> links_{i}{{{{\n")
        for link in ref["links"]:
            result.append(f"  {{{cpp(link['label'])}, {cpp(link['url'])}}},\n")
        result.append("}};\n")
    result.append(f"\ninline constexpr std::array<Reference, {len(refs)}> references{{{{\n")
    for i, ref in enumerate(refs):
        values = ", ".join(cpp(ref[key]) for key in ("id", "authors", "title", "publication"))
        result.append(f"  {{{values}, {ref['year'] or 0}, links_{i}}},\n")
    result.append("}};\n\n")
    indices = {ref["id"]: i for i, ref in enumerate(refs)}
    for tool in data["tools"]:
        result.append(f"inline constexpr std::array<Use, {len(tool['references'])}> uses_{tool['id']}{{{{\n")
        for use in tool["references"]:
            result.append(f"  {{&references[{indices[use['id']]}], {cpp(use['context'])}}},\n")
        result.append("}};\n")
    result.append("\nenum class Tool { " + ", ".join(t["id"] for t in data["tools"]) + " };\n")
    result.append("\n[[nodiscard]] constexpr std::span<Use const> for_tool(Tool tool)\n{\n  switch (tool)\n  {\n")
    for tool in data["tools"]:
        result.append(f"    case Tool::{tool['id']}: return uses_{tool['id']};\n")
    result.append("""  }
  throw std::invalid_argument("unknown citation tool");
}

[[nodiscard]] constexpr Reference const* find(std::string_view id)
{
  for (auto const& reference : references)
    if (reference.id == id) return &reference;
  return nullptr;
}
} // namespace bethe::citations
// clang-format on
""")
    return "".join(result)


def markdown(text):
    return re.sub(r"([\\`*_\[\]<>])", r"\\\1", text)


def bibliography(data):
    result = [BEGIN + "\n", "\n## Bibliography\n\n",
              "Generated from [data/citations.json](data/citations.json); edit the registry, not this section.\n"]
    for ref in data["references"]:
        result.append(f"\n### {ref['id']}\n\n")
        year = f" ({ref['year']})" if ref["year"] is not None else ""
        result.append(f"{markdown(ref['authors'])}. *{markdown(ref['title'])}*.\n")
        result.append(f"{markdown(ref['publication'])}{year}.\n\n")
        result.append(", ".join(f"[{markdown(link['label'])}](<{link['url']}>)" for link in ref["links"]) + ".\n")
        uses = [(tool, use) for tool in data["tools"] for use in tool["references"] if use["id"] == ref["id"]]
        if uses:
            result.append("\nRelevant tool modes:\n\n")
            for tool, use in uses:
                result.append(f"- `{tool['executable']}`: {markdown(use['context'])}\n")
    result.append("\n" + END)
    return "".join(result)


def update_markdown(original, data):
    if original.count(BEGIN) != 1 or original.count(END) != 1 or original.index(BEGIN) >= original.index(END):
        raise ValueError("CITATIONS.md must contain one ordered pair of bibliography markers")
    start, end = original.index(BEGIN), original.index(END) + len(END)
    return original[:start] + bibliography(data) + original[end:]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail on stale generated files without writing")
    args = parser.parse_args()
    try:
        data = json.loads((ROOT / "data/citations.json").read_text(encoding="utf-8"))
        validate(data)
        doc = ROOT / "CITATIONS.md"
        outputs = {
            ROOT / "include/bethe/citations.hpp": header(data),
            doc: update_markdown(doc.read_text(encoding="utf-8"), data),
        }
        stale = []
        for path, content in outputs.items():
            if not path.exists() or path.read_text(encoding="utf-8") != content:
                stale.append(str(path.relative_to(ROOT)))
                if not args.check:
                    path.write_text(content, encoding="utf-8")
        if stale and args.check:
            print("Stale citation files: " + ", ".join(stale) + "; run python3 scripts/generate_citations.py", file=sys.stderr)
            return 1
        return 0
    except (OSError, ValueError) as error:
        print(f"citations: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
