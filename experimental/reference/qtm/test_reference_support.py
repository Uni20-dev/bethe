# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Isolated import/provenance checks; no numerical quadrature is needed."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REFERENCE_DIR = Path(__file__).resolve().parent
REPO_ROOT = REFERENCE_DIR.parents[2]


class ReferenceSourceTests(unittest.TestCase):
    def run_isolated(self, directory, program):
        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join((str(directory), str(REFERENCE_DIR)))
        return subprocess.run([sys.executable, '-c', program], cwd=directory,
                              env=env, text=True, capture_output=True, timeout=20)

    def test_foreign_search_path_cannot_shadow_numerical_dependencies(self):
        with tempfile.TemporaryDirectory() as directory:
            for name in ('reference_xxx_vertex_density', 'reference_hypergeometric_tail'):
                (Path(directory) / (name + '.py')).write_text(
                    "raise AssertionError('foreign reference executed')\n")
            result = self.run_isolated(directory, '''
import json
import reference_xxx_qtm_vertex
from _reference_support import source_hashes
print(json.dumps(source_hashes()))
''')
            self.assertEqual(result.returncode, 0, result.stderr)
            sources = json.loads(result.stdout)
            for name in ('reference_xxx_vertex_density', 'reference_hypergeometric_tail'):
                path = 'scripts/' + name + '.py'
                self.assertIn(path, sources)
            self.assertIn('experimental/reference/qtm/reference_xxx_qtm.py', sources)
            import hashlib
            for name, digest in sources.items():
                self.assertEqual(digest, hashlib.sha256((REPO_ROOT / name).read_bytes()).hexdigest())

    def test_already_imported_foreign_dependencies_are_rejected(self):
        # Include a direct control, a transitive control and the QTM evaluator.
        for name in ('reference_xxx_vertex_density', 'reference_hypergeometric_tail',
                     'reference_xxx_qtm'):
            with self.subTest(module=name), tempfile.TemporaryDirectory() as directory:
                (Path(directory) / (name + '.py')).write_text('foreign = True\n')
                result = self.run_isolated(directory, f'''
import {name}
import reference_xxx_qtm_vertex
''')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('ImportError: QTM reference module ' + name, result.stderr)
                self.assertIn('same Bethe checkout', result.stderr)

    def test_provenance_rejects_foreign_module_loaded_after_initial_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_isolated(directory, '''
import sys
from types import ModuleType
import reference_xxx_qtm_vertex
from _reference_support import source_hashes
foreign = ModuleType('reference_xxx_vertex_density')
foreign.__file__ = 'reference_xxx_vertex_density.py'
sys.modules[foreign.__name__] = foreign
source_hashes()
''')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('ImportError: QTM reference module reference_xxx_vertex_density', result.stderr)


if __name__ == '__main__':
    unittest.main()
