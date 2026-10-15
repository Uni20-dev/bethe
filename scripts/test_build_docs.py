import subprocess
import html
import tempfile
from pathlib import Path
import unittest

from build_docs import adapt_links, adapt_math, site_sources
from check_docs_site import Page, bibliography_link_errors
from tutorial_links import BEGIN, END, render_tutorial_link, sync_tutorial_links, tutorial_urls


class SiteTests(unittest.TestCase):
    def test_inline_and_block_math(self):
        self.assertEqual(html.unescape(adapt_math('Energy $`E=\\sum k^2`$.\n')),
                         'Energy <span class="arithmatex">\\(E=\\sum k^2\\)</span>.\n')
        source = '```math\n\\begin{aligned}\na&=b,\\\\{}\nc&=d.\n\\end{aligned}\n```\n'
        result = html.unescape(adapt_math(source))
        self.assertIn('a&=b,\\\\{}', result)
        self.assertTrue(result.startswith('\n<div class="arithmatex">\\[\n'))
        self.assertTrue(result.endswith('\\]</div>\n\n'))

    def test_code_examples_unchanged(self):
        for source in ('```sh\necho "$`x`$"\n```\n',
                       '````md\n```math\nx\n```\n````\n',
                       '~~~text\n$`x`$\n~~~\n',
                       'Use ``$`x`$`` for inline math.\n'):
            with self.subTest(source=source):
                self.assertEqual(adapt_math(source), source)

    def test_only_tracked_documentation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            subprocess.run(['git', 'init', '-q', temp], check=True)
            (root / 'docs').mkdir()
            (root / 'docs/public.md').write_text('Public')
            (root / 'docs/private.md').write_text('Never publish')
            (root / 'source.cpp').write_text('// not documentation')
            subprocess.run(['git', 'add', 'docs/public.md', 'source.cpp'], cwd=root, check=True)
            self.assertEqual(site_sources(root), [Path('docs/public.md')])

    def test_links_to_repository_assets(self):
        included = {Path('README.md'), Path('docs/output.md')}
        text = ('[Output](output.md#tables) [Home](../README.md) '
                '[Source](../scripts/example.py) [Web](https://example.com/)\n'
                '```md\n[Example](elsewhere.md)\n```\n'
                '`[Inline](example.md)`\n')
        result = adapt_links(text, Path('docs/page.md'), included)
        self.assertIn('[Output](output.md#tables)', result)
        self.assertIn('[Home](../README.md)', result)
        self.assertIn('[Source](https://github.com/Uni20-dev/bethe/blob/main/scripts/example.py)', result)
        self.assertIn('```md\n[Example](elsewhere.md)\n```', result)
        self.assertIn('`[Inline](example.md)`', result)

    def test_external_angle_bracket_links_are_unchanged(self):
        for destination in ('https://arxiv.org/abs/cond-mat/0012439',
                            'https://doi.org/10.1016/0550-3213(95)00105-2',
                            'https://example.org/a(b(c))?q=x&v=2#section',
                            '//example.org/paper', 'mailto:author@example.org'):
            source = f'[Reference](<{destination}>).\n'
            with self.subTest(destination=destination):
                self.assertEqual(adapt_links(source, Path('CITATIONS.md'), set()), source)

    def test_angle_bracket_local_links_and_code(self):
        source = ('[Guide](<docs/output.md#tables>) [Here](<#bibliography>) '
                  '[Paper](<papers/a (2020).pdf?raw=1#page=2>)\n'
                  '`[Code](<example.py>)`\n```md\n[Code](<example.py>)\n```\n')
        result = adapt_links(source, Path('CITATIONS.md'), {Path('docs/output.md')})
        self.assertIn('[Guide](<docs/output.md#tables>) [Here](<#bibliography>)', result)
        self.assertIn('[Paper](<https://github.com/Uni20-dev/bethe/blob/main/papers/a%20%282020%29.pdf?raw=1#page=2>)', result)
        self.assertIn('`[Code](<example.py>)`\n```md\n[Code](<example.py>)\n```', result)

    def test_bibliography_checker_rejects_rewritten_external_urls(self):
        url = 'https://doi.org/10.1016/0550-3213(95)00105-2'
        refs = [{'links': [{'url': url}]}]
        self.assertEqual(bibliography_link_errors(Page(f'<a href="{url}">DOI</a>'), refs), [])
        for wrong in ('https://github.com/Uni20-dev/bethe/blob/main/%3Chttps%3A/doi.org/10.1016/0550-3213%2895',
                      'https://doi.org/10.1016/0550-3213(95', ''):
            with self.subTest(url=wrong):
                self.assertTrue(bibliography_link_errors(Page(f'<a href="{wrong}">DOI</a>'), refs))
        # Repeated registry URLs must not silently lose one of their links.
        self.assertTrue(bibliography_link_errors(Page(f'<a href="{url}">DOI</a>'), refs * 2))


class TutorialLinkTests(unittest.TestCase):
    project = {'site_url': 'https://example.org/bethe/',
               'repo_url': 'https://github.com/example/bethe'}

    def test_page_and_index_urls(self):
        cases = [
            ('docs/tutorials/xxz-spinons.md', 'docs/tutorials/xxz-spinons/'),
            ('docs/tutorials/index.md', 'docs/tutorials/'),
            ('docs/tutorials/nested/index.md', 'docs/tutorials/nested/'),
            ('docs/tutorials/a (b).md', 'docs/tutorials/a%20%28b%29/'),
        ]
        for relative, published in cases:
            with self.subTest(relative=relative):
                site, source = tutorial_urls(Path(relative), self.project)
                self.assertEqual(site, 'https://example.org/bethe/' + published)
                self.assertEqual(source, 'https://github.com/example/bethe/blob/main/'
                                 + relative.replace(' ', '%20').replace('(', '%28').replace(')', '%29'))

    def test_source_banner_preserves_prose_and_is_idempotent(self):
        body = 'Energy $`E=\\sum k^2`$.\n\n```math\na=b\n```\n\n[Figure](figures/example.svg)\n'
        original = '# Example\n\n' + body
        relative = Path('docs/tutorials/example.md')
        result = render_tutorial_link(original, relative, self.project)
        self.assertEqual(result, '# Example\n\n' + BEGIN + '\n'
                         '[Read this page on the documentation site]'
                         '(https://example.org/bethe/docs/tutorials/example/)\n'
                         + END + '\n\n' + body)
        self.assertEqual(render_tutorial_link(result, relative, self.project), result)

    def test_staging_reverses_link_to_exact_source_file(self):
        for name in ('example.md', 'index.md'):
            relative = Path('docs/tutorials') / name
            source = render_tutorial_link('# Tutorial\n\nBody.\n', relative, self.project)
            staged = render_tutorial_link(source, relative, self.project, for_site=True)
            self.assertIn('[View this page on GitHub]'
                          f'(https://github.com/example/bethe/blob/main/docs/tutorials/{name})', staged)
            self.assertNotIn('https://example.org/', staged)
            self.assertTrue(staged.endswith('\n\nBody.\n'))
            self.assertEqual(render_tutorial_link(staged, relative, self.project), source)

    def test_renaming_refreshes_existing_link(self):
        old = render_tutorial_link('# Tutorial\n\nBody.\n', Path('docs/tutorials/old.md'), self.project)
        updated = render_tutorial_link(old, Path('docs/tutorials/new.md'), self.project)
        self.assertIn('https://example.org/bethe/docs/tutorials/new/', updated)
        self.assertNotIn('/old/', updated)
        self.assertEqual(updated.count(BEGIN), 1)

    def test_malformed_blocks_are_rejected_without_consuming_prose(self):
        cases = [
            '# Title\n\n' + BEGIN + '\nUnterminated.\n',
            '# Title\n\n' + END + '\nOrphan end marker.\n',
            '# Title\n\nBody.\n' + BEGIN + '\nLink\n' + END + '\n',
            '# Title\n\n' + BEGIN + '\n' + BEGIN + '\nLink\n' + END + '\n',
            'No heading\n\nBody.\n',
        ]
        for source in cases:
            with self.subTest(source=source), self.assertRaises(ValueError):
                render_tutorial_link(source, Path('docs/tutorials/example.md'), self.project)

    def test_check_does_not_write_and_generator_ignores_untracked_drafts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            subprocess.run(['git', 'init', '-q', temp], check=True)
            (root / 'zensical.toml').write_text('[project]\nsite_url = "https://example.org/bethe/"\n'
                                              'repo_url = "https://github.com/example/bethe"\n')
            (root / 'docs/tutorials').mkdir(parents=True)
            tracked = Path('docs/tutorials/public.md')
            draft = root / 'docs/tutorials/private.md'
            original = '# Tutorial\n\nBody.\n'
            (root / tracked).write_text(original)
            draft.write_text(original)
            subprocess.run(['git', 'add', str(tracked)], cwd=root, check=True)
            self.assertEqual(sync_tutorial_links(root, check=True), [tracked])
            self.assertEqual((root / tracked).read_text(), original)
            self.assertEqual(sync_tutorial_links(root), [tracked])
            self.assertEqual(sync_tutorial_links(root, check=True), [])
            self.assertEqual(draft.read_text(), original)
            # A copied banner still points at the old filename until refreshed.
            new = Path('docs/tutorials/renamed.md')
            subprocess.run(['git', 'mv', str(tracked), str(new)], cwd=root, check=True)
            self.assertEqual(sync_tutorial_links(root, check=True), [new])


if __name__ == '__main__':
    unittest.main()
