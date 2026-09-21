# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Generator tests use isolated temporary files, never modifying the checkout."""

import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location(
    "generate_citations", Path(__file__).resolve().parents[1] / "scripts/generate_citations.py"
)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


class CitationTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads((generator.ROOT / "data/citations.json").read_text(encoding="utf-8"))

    def test_registry(self):
        generator.validate(self.data)
        output = generator.header(self.data)
        self.assertIn("Müller", output)
        self.assertIn("uses_hubbard_pbc", output)
        bibliography = generator.bibliography(self.data)
        # A library-only citation is still part of the registry and bibliography.
        self.assertIn("### caux-xxz-spinons", bibliography)
        self.assertNotIn("caux-xxz-spinons", str(self.data["tools"]))

    def test_validation(self):
        mutations = [
            lambda d: d.update(schema_version=2),
            lambda d: d["references"].append(copy.deepcopy(d["references"][0])),
            lambda d: d["references"][0].update(year=True),
            lambda d: d["references"][0].update(title="line\nbreak"),
            lambda d: d["references"][0].update(links=[]),
            lambda d: d["references"][0]["links"][0].update(url="javascript:alert(1)"),
            lambda d: d["tools"][0]["references"][0].update(id="missing"),
            lambda d: d["tools"][0]["references"].append(copy.deepcopy(d["tools"][0]["references"][0])),
            lambda d: d["tools"][0].update(id="invalid-id"),
            lambda d: d["tools"][1].update(executable=d["tools"][0]["executable"]),
            lambda d: d["references"][0].update(unknown="typo"),
        ]
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                data = copy.deepcopy(self.data)
                mutate(data)
                with self.assertRaises(ValueError):
                    generator.validate(data)

    def test_escaping(self):
        title = 'A "quoted" \\ path with *stars* and [brackets]'
        self.data["references"][0]["title"] = title
        generator.validate(self.data)
        self.assertIn('A \\"quoted\\" \\\\ path', generator.header(self.data))
        self.assertIn(r"\*stars\* and \[brackets\]", generator.bibliography(self.data))

    def test_provenance_preserved_and_idempotent(self):
        original = "Hand-written prefix.\n" + generator.BEGIN + "\nstale\n" + generator.END + "\nHand-written suffix.\n"
        updated = generator.update_markdown(original, self.data)
        self.assertTrue(updated.startswith("Hand-written prefix.\n"))
        self.assertTrue(updated.endswith("\nHand-written suffix.\n"))
        self.assertNotIn("stale", updated)
        self.assertEqual(updated, generator.update_markdown(updated, self.data))
        for malformed in ("no markers", generator.END + generator.BEGIN, original + generator.BEGIN):
            with self.assertRaises(ValueError):
                generator.update_markdown(malformed, self.data)

    def test_check_mode_is_read_only_and_detects_both_outputs(self):
        with tempfile.TemporaryDirectory(prefix="bethe-citation-test-") as directory:
            root = Path(directory)
            (root / "data").mkdir()
            (root / "include/bethe").mkdir(parents=True)
            (root / "data/citations.json").write_text(json.dumps(self.data), encoding="utf-8")
            doc = root / "CITATIONS.md"
            header = root / "include/bethe/citations.hpp"
            doc.write_text("Provenance\n" + generator.BEGIN + "\n" + generator.END + "\n", encoding="utf-8")
            with patch.object(generator, "ROOT", root), contextlib.redirect_stderr(io.StringIO()):
                with patch.object(sys, "argv", ["generator", "--check"]):
                    before = doc.read_bytes()
                    self.assertEqual(generator.main(), 1)
                    self.assertEqual(doc.read_bytes(), before)
                    self.assertFalse(header.exists())
                with patch.object(sys, "argv", ["generator"]):
                    self.assertEqual(generator.main(), 0)
                with patch.object(sys, "argv", ["generator", "--check"]):
                    self.assertEqual(generator.main(), 0)
                    for path in (doc, header):
                        original = path.read_bytes()
                        path.write_bytes(original.replace(b"M\xc3\xbcller", b"stale-author"))
                        stale = path.read_bytes()
                        self.assertEqual(generator.main(), 1)
                        self.assertEqual(path.read_bytes(), stale)
                        path.write_bytes(original)
                    self.assertEqual(generator.main(), 0)


if __name__ == "__main__":
    unittest.main()
