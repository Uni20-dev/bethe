"""Regression tests for GitHub-specific math restrictions (no network needed)."""

import unittest

from check_markdown_math import source_expressions


class MarkdownMathTest(unittest.TestCase):
    def test_display_angle_brackets_rejected(self):
        for relation in ("<", ">"):
            with self.subTest(relation=relation), self.assertRaisesRegex(ValueError, "all math"):
                source_expressions("test.md", f"```math\n\\sum_{{i{relation}j}} x_i\n```\n")

    def test_operatorname_rejected_in_both_forms(self):
        for source in (r"$`\operatorname{atan2}(x,y)`$",
                       "```math\n" + r"\operatorname{atan2}(x,y)" + "\n```"):
            with self.subTest(source=source), self.assertRaisesRegex(ValueError, "GitHub rejects"):
                source_expressions("test.md", source)

    def test_safe_display_and_inline(self):
        equation = "\n".join((r"\begin{aligned}",
                              r"H&=\sum_{i\lt j}x_i,\\{}",
                              r"E&=\mathrm{atan2}(x,y).", r"\end{aligned}"))
        source = "```math\n" + equation + "\n```\n" + r"$`x\gt 0`$"
        self.assertEqual(source_expressions("test.md", source),
                         ["$$" + equation + "$$", r"$x\gt 0$"])

    def test_code_examples_are_not_math(self):
        self.assertEqual(source_expressions("test.md", "```cpp\nif (i<j) {}\n```"), [])

    def test_bare_row_break_rejected(self):
        with self.assertRaisesRegex(ValueError, "end math rows"):
            source_expressions("test.md", "```math\na&=b" + "\\\\" + "\nc&=d\n```")


if __name__ == "__main__":
    unittest.main()
