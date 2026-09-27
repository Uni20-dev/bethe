import subprocess
import html
import tempfile
from pathlib import Path
import unittest

from build_docs import adapt_links, adapt_math, site_sources


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


if __name__ == '__main__':
    unittest.main()
