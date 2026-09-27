#!/usr/bin/env python3
"""Build Pages from Git-tracked sources, adapting GitHub math only in staging.

Run from any directory. New documentation must be git-added before building;
untracked local drafts are deliberately never published. Dependencies are in
docs/requirements-site.txt. Numerical data/figures are NOT regenerated here.
"""

import argparse
import json
import posixpath
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import quote, urlsplit, unquote

ROOT = Path(__file__).resolve().parents[1]


def math_entities(text):
    # Protect TeX from Markdown's escape/emphasis handling, including inside
    # HTML blocks that follow prose containing literal less-than signs.
    return ''.join(c if c == '\n' else f'&#{ord(c)};' for c in text)


def adapt_math(text):
    """Preserve code fences/spans; translate protected inline and fenced math."""
    result = []
    fence = None
    is_math = False
    for line in text.splitlines(keepends=True):
        match = re.match(r"^ {0,3}(`{3,}|~{3,})(.*)\n?$", line)
        if match and fence is None:
            fence = match[1]
            is_math = match[2].strip() == "math"
            result.append('\n<div class="arithmatex">' + math_entities('\\[\n') if is_math else line)
        elif (match and fence and match[1][0] == fence[0]
              and len(match[1]) >= len(fence) and not match[2].strip()):
            result.append(math_entities('\\]') + '</div>\n\n' if is_math else line)
            fence = None
            is_math = False
        elif fence:
            result.append(math_entities(line) if is_math else line)
        else:
            # Ordinary inline code wins as a complete token, so examples that
            # demonstrate math markup in a code span stay literal.
            result.append(re.sub(
                r"\$`([^`\n]+)`\$|(`+)(?!`)(.*?)(?<!`)\2(?!`)",
                # Entity-encode the inline TeX: Python-Markdown otherwise
                # processes escapes and emphasis even inside an HTML span.
                lambda m: ('<span class="arithmatex">'
                           + math_entities('\\(' + m[1] + '\\)')
                           + '</span>')
                if m[1] is not None else m[0], line))
    return "".join(result)


def site_sources(root):
    paths = subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=root, text=True).split("\0")
    return [Path(p) for p in paths if p and (
        p in {"README.md", "CITATIONS.md"} or p.startswith("docs/"))]


def adapt_links(text, relative, included):
    """Source/scripts/papers not published as site pages remain GitHub links."""
    def replace(match):
        url = urlsplit(match[1])
        if url.scheme or url.netloc or not url.path or url.path.startswith('/'):
            return match[0]
        target = posixpath.normpath(str(relative.parent / unquote(url.path)))
        if Path(target) in included:
            return match[0]
        suffix = ('#' + url.fragment) if url.fragment else ''
        return '](https://github.com/Uni20-dev/bethe/blob/main/' + quote(target) + suffix + ')'
    # Repository links have simple, unquoted destinations; code fences are
    # deliberately left alone, as are inline code examples.
    lines = []
    fence = None
    for line in text.splitlines(keepends=True):
        match = re.match(r'^ {0,3}(`{3,}|~{3,})(.*)\n?$', line)
        if match and fence is None:
            fence = match[1]
        elif (match and fence and match[1][0] == fence[0]
              and len(match[1]) >= len(fence) and not match[2].strip()):
            fence = None
        elif fence is None:
            line = re.sub(r'\]\(([^\s)]+)\)|(`+)(?!`)(.*?)(?<!`)\2(?!`)',
                          lambda m: replace(m) if m[1] is not None else m[0], line)
        lines.append(line)
    return ''.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path,
                        default=ROOT / "build_codex/docs-site")
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="staging-", dir=work) as temp:
        staging = Path(temp)
        sources = site_sources(ROOT)
        included = set(sources)
        for relative in sources:
            source = ROOT / relative
            target = staging / "source" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.suffix == ".md":
                target.write_text(adapt_math(adapt_links(source.read_text(), relative, included)))
            else:
                shutil.copyfile(source, target)
        # Preserve README.md as a link target for existing GitHub-first docs.
        shutil.copyfile(staging / "source/README.md", staging / "source/index.md")
        config = (ROOT / "zensical.toml").read_text().replace(
            'docs_dir = "source"', 'docs_dir = ' + json.dumps(staging.name + '/source'))
        (work / "zensical.toml").write_text(config)
        subprocess.run([sys.executable, "-m", "zensical", "build"],
                       cwd=work, check=True)
    print(f"Site: {work / 'site'}")


if __name__ == "__main__":
    main()
