#!/usr/bin/env python3
"""Maintain links between GitHub tutorials and their published pages."""

import argparse
from pathlib import Path
import subprocess
import tomllib
from urllib.parse import quote


ROOT = Path(__file__).resolve().parents[1]
BEGIN = '<!-- BEGIN GENERATED PAGE LINK -->'
END = '<!-- END GENERATED PAGE LINK -->'


def is_tutorial(relative):
    return relative.suffix == '.md' and relative.is_relative_to('docs/tutorials')


def tutorial_urls(relative, project):
    """Zensical uses directory URLs; index.md maps to its parent directory."""
    page = relative.parent if relative.name == 'index.md' else relative.with_suffix('')
    published = project['site_url'].rstrip('/') + '/' + quote(page.as_posix()) + '/'
    source = project['repo_url'].rstrip('/') + '/blob/main/' + quote(relative.as_posix())
    return published, source


def render_tutorial_link(text, relative, project, *, for_site=False):
    """Replace only the generated block below the title, preserving the prose."""
    title, newline, body = text.partition('\n')
    if not title.startswith('# ') or not newline:
        raise ValueError(f'{relative}: expected a Markdown title on the first line')
    if text.count(BEGIN) > 1 or text.count(END) > 1:
        raise ValueError(f'{relative}: duplicate generated page link markers')
    body = body.lstrip('\n')
    if body.startswith(BEGIN + '\n'):
        _, end, body = body.partition('\n' + END + '\n')
        if not end:
            raise ValueError(f'{relative}: unterminated generated page link')
        body = body.lstrip('\n')
    if BEGIN in body or END in body:
        raise ValueError(f'{relative}: misplaced or duplicate generated page link markers')
    published, source = tutorial_urls(relative, project)
    label, url = (('View this page on GitHub', source) if for_site else
                  ('Read this page on the documentation site', published))
    return f'{title}\n\n{BEGIN}\n[{label}]({url})\n{END}\n\n{body}'


def sync_tutorial_links(root, *, check=False):
    """Check/update tracked tutorials only; --check never writes source files."""
    project = tomllib.loads((root / 'zensical.toml').read_text())['project']
    tracked = subprocess.check_output(
        ['git', 'ls-files', '-z', '--', 'docs/tutorials'], cwd=root, text=True)
    updates = []
    for name in tracked.split('\0'):
        relative = Path(name)
        if not name or not is_tutorial(relative):
            continue
        original = (root / relative).read_text()
        updated = render_tutorial_link(original, relative, project)
        if updated != original:
            updates.append((relative, updated))
    if not check:
        for relative, text in updates:
            (root / relative).write_text(text)
    return [relative for relative, _ in updates]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true',
                        help='fail on missing or stale links without modifying files')
    args = parser.parse_args()
    try:
        changed = sync_tutorial_links(ROOT, check=args.check)
    except ValueError as error:
        parser.exit(1, f'{error}\n')
    if args.check and changed:
        parser.exit(1, 'Missing or stale tutorial page links:\n'
                    + ''.join(f'  {path}\n' for path in changed)
                    + 'Run python3 scripts/tutorial_links.py to refresh them.\n')
    print('Tutorial page links are up to date.' if args.check else
          f'Updated page links in {len(changed)} tutorials.')


if __name__ == '__main__':
    main()
