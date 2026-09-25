#!/usr/bin/env python3
"""Check GitHub-safe math markup; optionally verify GitHub's Markdown output.

Run from the repository root. --github requires authenticated gh and network
access; it submits Markdown to the rendering API without publishing anything.
This checks the Markdown-to-TeX boundary, not mathematical correctness or the
complete MathJax grammar.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess


class MathReader(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.expressions = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        if tag == "math-renderer":
            self.current = ""

    def handle_data(self, data):
        if self.current is not None:
            self.current += data

    def handle_endtag(self, tag):
        if tag == "math-renderer" and self.current is not None:
            self.expressions.append(self.current)
            self.current = None


def source_expressions(path, text):
    expressions = []
    fence = None
    block = []
    for number, line in enumerate(text.splitlines(), 1):
        if line.startswith("```"):
            if fence == "math":
                expressions.append("$$" + "\n".join(block) + "$$")
                block = []
            fence = line[3:] if fence is None else None
        elif fence == "math":
            if line.endswith("\\\\"):
                raise ValueError(f"{path}:{number}: end math rows with \\\\{{}}, not bare \\\\")
            block.append(line)
        elif fence is None:
            for match in re.finditer(r"\$([^$\n]+)\$", line):
                content = match[1]
                if not (content.startswith("`") and content.endswith("`")):
                    raise ValueError(f"{path}:{number}: use protected $`...`$ inline math")
                content = content[1:-1]
                if "<" in content or ">" in content:
                    raise ValueError(f"{path}:{number}: use \\lt or \\gt in inline math")
                expressions.append("$" + content + "$")
    if fence is not None:
        raise ValueError(f"{path}: unclosed code fence")
    return expressions


def check_github(item):
    path, text, expected = item
    result = subprocess.run(
        ["gh", "api", "markdown", "--input", "-"],
        input=json.dumps({"text": text, "mode": "gfm"}),
        text=True, capture_output=True, check=True,
    )
    parser = MathReader()
    parser.feed(result.stdout)
    if parser.expressions != expected:
        for index, (source, rendered) in enumerate(zip(expected, parser.expressions), 1):
            if source != rendered:
                raise ValueError(f"{path}: expression {index} changed\nsource: {source!r}\nGitHub: {rendered!r}")
        raise ValueError(f"{path}: expected {len(expected)} expressions, got {len(parser.expressions)}")
    return len(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--github", action="store_true")
    args = parser.parse_args()
    paths = subprocess.check_output(["git", "ls-files", "-z", "*.md"], text=True).split("\0")
    items = []
    for path in filter(None, paths):
        text = Path(path).read_text()
        expressions = source_expressions(path, text)
        if expressions:
            items.append((path, text, expressions))
    if args.github:
        with ThreadPoolExecutor(max_workers=4) as pool:
            count = sum(pool.map(check_github, items))
    else:
        count = sum(len(item[2]) for item in items)
    print(f"Checked {count} expressions in {len(items)} files"
          + (" through GitHub's Markdown renderer." if args.github else "."))


if __name__ == "__main__":
    main()
