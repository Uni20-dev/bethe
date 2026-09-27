import math
import unittest

from plot_xxz_tutorial import DATA, POINTS, load_cases, read_export


class TutorialDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_isotropic_normalization(self):
        for row in self.cases['xxz-delta1']['spinon']:
            self.assertAlmostEqual(row['energy'], math.pi / 2 * math.sin(row['p']), places=13)

    def test_gap_and_two_particle_threshold(self):
        rows = self.cases['xxz-delta2']
        gap = rows['spinon'][0]['energy']
        self.assertAlmostEqual(gap, 0.19490113571001469, places=13)
        self.assertAlmostEqual(min(r['lower'] for r in rows['two-spinon']), 2 * gap, places=13)

    def test_folding_is_union_of_shifted_branches(self):
        raw = self.cases['xxz-delta2']['two-spinon']
        folded = self.cases['xxz-delta2-folded']['two-spinon']
        for i, row in enumerate(folded):
            shifted = raw[(i + (POINTS - 1) // 2) % (POINTS - 1)]
            self.assertAlmostEqual(row['lower'], min(raw[i]['lower'], shifted['lower']), places=12)
            self.assertAlmostEqual(row['upper'], max(raw[i]['upper'], shifted['upper']), places=12)

    def test_failed_missing_and_wrong_convention_data_rejected(self):
        text = (DATA / 'xxz-delta2.csv').read_text()
        lines = text.splitlines()
        row_index = next(i for i, line in enumerate(lines) if line.startswith('spinon,'))
        missing_energy = lines.copy()
        fields = missing_energy[row_index].split(',')
        fields[4] = ''
        missing_energy[row_index] = ','.join(fields)
        nan_energy = lines.copy()
        fields[4] = 'nan'
        nan_energy[row_index] = ','.join(fields)
        for invalid in (text.replace('# Outcome: success', '# Outcome: failure'),
                        text.replace(',converged', ',precision_limit', 1),
                        '\n'.join(lines[:row_index] + lines[row_index + 1:]),
                        '\n'.join(missing_energy), '\n'.join(nan_energy)):
            with self.subTest(invalid=invalid[:70]), self.assertRaises(ValueError):
                read_export(invalid, 2, False)
        with self.assertRaises(ValueError):
            read_export(text, 2, True)


if __name__ == '__main__':
    unittest.main()
