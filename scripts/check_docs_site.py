#!/usr/bin/env python3
"""Check built local links/assets and exact source-to-HTML math preservation."""

import argparse
from collections import Counter
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

from build_docs import ROOT, site_sources
from check_markdown_math import source_expressions


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.math = []
        self.current_math = None
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        for key in ('href', 'src'):
            if attrs.get(key):
                self.links.append(attrs[key])
        if 'arithmatex' in attrs.get('class', '').split():
            self.current_math = ''

    def handle_data(self, data):
        if self.current_math is not None:
            self.current_math += data

    def handle_endtag(self, tag):
        if self.current_math is not None and tag in ('span', 'div'):
            self.math.append(self.current_math.strip())
            self.current_math = None


def bibliography_link_errors(page, references):
    """Check final hrefs against the registry, not just local link existence."""
    expected = Counter(link['url'] for ref in references for link in ref['links'])
    missing = expected - Counter(page.links)
    return [f'CITATIONS.md: missing {count} rendered link(s) to {url}'
            for url, count in missing.items()]


def check_site(site):
    site = site.resolve()
    pages = {p: Page(p.read_text()) for p in site.rglob('*.html')}
    errors = []
    for path, page in pages.items():
        for link in page.links:
            # Zensical's generated 404 template has a skip link but no article.
            if path.name == '404.html' and link == '#__skip':
                continue
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            target = unquote(url.path)
            if target.startswith('/bethe/'):
                dest = site / target.removeprefix('/bethe/')
            elif target.startswith('/'):
                dest = site / target.lstrip('/')
            else:
                dest = path.parent / target if target else path
            if dest.is_dir():
                dest /= 'index.html'
            dest = dest.resolve()
            if not dest.is_file():
                errors.append(f'{path.relative_to(site)}: missing {link}')
            elif url.fragment and dest in pages and unquote(url.fragment) not in pages[dest].ids:
                errors.append(f'{path.relative_to(site)}: missing anchor {link}')
    references = json.loads((ROOT / 'data/citations.json').read_text())['references']
    bibliography = pages.get((site / 'CITATIONS/index.html').resolve())
    if bibliography is not None:
        errors.extend(bibliography_link_errors(bibliography, references))
    count = 0
    for relative in site_sources(ROOT):
        if relative.suffix != '.md':
            continue
        expected = source_expressions(relative, (ROOT / relative).read_text())
        expected = [('\\[\n' + s[2:-2] + '\n\\]') if s.startswith('$$')
                    else ('\\(' + s[1:-1] + '\\)') for s in expected]
        path = (site / relative.with_suffix('') / 'index.html').resolve()
        if relative.name == 'index.md':
            path = site / relative.parent / 'index.html'
        if relative == Path('README.md'):
            path = site / 'index.html'
        if path not in pages or pages[path].math != expected:
            errors.append(f'{relative}: source and rendered math differ')
        count += len(expected)
    if errors:
        raise ValueError('\n'.join(errors))
    links = sum(len(ref['links']) for ref in references)
    print(f'Checked links/assets in {len(pages)} pages, {links} bibliography URLs, and preserved {count} equations.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site', type=Path, nargs='?', default=ROOT / 'build_codex/docs-site/site')
    check_site(parser.parse_args().site)
