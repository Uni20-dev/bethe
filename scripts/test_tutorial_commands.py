"""Keep runnable examples independent of the developer's build layout."""

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from tutorial_common import solver_executable


ROOT = Path(__file__).resolve().parents[1]


class TutorialCommandTests(unittest.TestCase):
    def test_solver_on_path_outside_working_directory(self):
        # Python is a real executable outside a fresh calculation directory.
        executable = Path(sys.executable).absolute()
        with tempfile.TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                with patch.dict(os.environ, {'PATH': str(executable.parent)}):
                    self.assertEqual(solver_executable(executable.name), executable.resolve())
            finally:
                os.chdir(previous)

    def test_explicit_absolute_and_dot_slash_paths(self):
        executable = Path(sys.executable).absolute()
        previous = Path.cwd()
        try:
            os.chdir(executable.parent)
            with patch.dict(os.environ, {'PATH': ''}):
                self.assertEqual(solver_executable(str(executable)), executable.resolve())
                self.assertEqual(solver_executable('./' + executable.name), executable.resolve())
        finally:
            os.chdir(previous)

    def test_missing_solver_is_an_argument_error(self):
        with patch.dict(os.environ, {'PATH': ''}):
            with self.assertRaises(argparse.ArgumentTypeError):
                solver_executable('no-such-bethe-executable')

    def test_existing_nonexecutable_is_rejected(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            solver_executable(str(ROOT / 'README.md'))

    def test_model_examples_do_not_assume_a_build_directory(self):
        tracked = subprocess.check_output(
            ['git', 'ls-files', 'README.md', 'docs/*.md'], cwd=ROOT, text=True).splitlines()
        for relative in tracked:
            text = (ROOT / relative).read_text()
            with self.subTest(path=relative):
                self.assertNotRegex(text, r'\bbuild(?:_codex)?/bethe-')
                self.assertNotIn('build_codex', text)

    def test_single_solver_scripts_use_path_aware_argument_type(self):
        for script in (ROOT / 'scripts').glob('plot_*tutorial.py'):
            text = script.read_text()
            if re.search(r'add_argument\([\"\x27]--solver[\"\x27]', text):
                with self.subTest(script=script.name):
                    self.assertRegex(text, r'add_argument\([\"\x27]--solver[\"\x27], type=solver_executable')


if __name__ == '__main__':
    unittest.main()
