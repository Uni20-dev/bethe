"""Reject plausible-looking corruptions of the continuum-density tutorial."""
import unittest
from plot_xxx_thermodynamic_structure_factor_tutorial import DATA, NQ, NW, read_case


class XXXThermodynamicTutorial(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = {name: (DATA/f'xxx-thermo-{name}.csv').read_text() for name in ('spectrum', 'continuum')}

    def test_grid(self):
        values = read_case(self.original)
        self.assertEqual(len(values), NQ)
        self.assertEqual(len(values[0]), NW)
        self.assertIsNone(values[(NQ-1)//2][0])

    def test_wrong_metadata(self):
        for old, new in [('# Channel: zz', '# Channel: raising'), ('# Failed points: 0', '# Failed points: 1'),
                         ('# Outcome: success', '# Outcome: partial')]:
            texts = dict(self.original)
            texts['spectrum'] = texts['spectrum'].replace(old, new)
            with self.assertRaises(ValueError): read_case(texts)

    def test_missing_and_corrupted_rows(self):
        for mode in ('missing', 'outside_weight', 'status', 'threshold', 'negative'):
            texts = dict(self.original)
            lines = texts['spectrum'].splitlines()
            start = next(i for i, line in enumerate(lines) if line.startswith('q,omega,'))+1
            if mode == 'missing': del lines[start+100]
            else:
                index = start + ((NQ-1)//2)*NW if mode == 'threshold' else start
                cells = lines[index].split(',')
                if mode == 'outside_weight': cells[2] = '1'
                elif mode == 'status': cells[-1] = 'numerical_failure'
                elif mode == 'threshold': cells[2] = '0'
                else: cells[2] = '-1'
                lines[index] = ','.join(cells)
            texts['spectrum'] = '\n'.join(lines)
            with self.assertRaises(ValueError): read_case(texts)


if __name__ == '__main__': unittest.main()
