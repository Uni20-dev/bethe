import math
import unittest

from plot_spin1_tutorial import DATA, POINTS, load_cases, read_export


class SpinOneTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_uls_maxima_and_thresholds(self):
        rows = self.cases['su3-unfolded']
        amplitude = 4*math.pi/(3*math.sqrt(3))
        self.assertAlmostEqual(max(r['energy'] for r in rows['3']), 1.5*amplitude, places=12)
        self.assertAlmostEqual(max(r['energy'] for r in rows['bar3']), .5*amplitude, places=12)
        mid = (POINTS-1)//2
        self.assertAlmostEqual(rows['two-soliton'][mid]['lower'], amplitude, places=12)
        self.assertAlmostEqual(rows['four-soliton'][mid]['lower'], amplitude/2, places=12)
        for branch in ('3', 'bar3'):
            self.assertEqual(rows[branch][0]['energy'], 0)
            self.assertEqual(rows[branch][-1]['energy'], 0)

    def test_tb_normalization_and_spectator_threshold(self):
        rows = self.cases['tb-unfolded']
        for r in rows['spinon']:
            self.assertAlmostEqual(r['energy'], 2*math.pi*math.sin(r['p']), places=12)
        for two, four in zip(rows['two-spinon'], rows['four-spinon']):
            self.assertAlmostEqual(two['lower'], four['lower'], places=12)
        self.assertEqual(rows['two-spinon'][0]['upper'], 0)
        self.assertAlmostEqual(rows['four-spinon'][0]['upper'], 8*math.pi, places=12)

    def test_folded_envelopes(self):
        for model, images, suffix in [('su3', 3, 'soliton'), ('tb', 2, 'spinon')]:
            raw, folded = self.cases[f'{model}-unfolded'], self.cases[f'{model}-folded']
            for particle in ('3', 'bar3') if model == 'su3' else ('spinon',):
                self.assertEqual(raw[particle], folded[particle])
            for count in ('two', 'four'):
                branch = f'{count}-{suffix}'
                for i, row in enumerate(folded[branch]):
                    aliases = [raw[branch][(i+j*(POINTS-1)//images) % (POINTS-1)] for j in range(images)]
                    self.assertAlmostEqual(row['lower'], min(r['lower'] for r in aliases), places=11)
                    self.assertAlmostEqual(row['upper'], max(r['upper'] for r in aliases), places=11)

    def test_invalid_data_rejected(self):
        text = (DATA / 'su3-unfolded.csv').read_text()
        lines = text.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('3,3,'))
        bad = [text.replace('# Outcome: success', '# Outcome: failure'),
               text.replace(',converged', ',precision_limit', 1),
               '\n'.join(lines[:index] + lines[index+1:])]
        for column, value in [(1, 'bar3'), (2, '1'), (4, '1'), (5, ''), (5, 'nan'), (6, '0')]:
            copy = lines.copy()
            row = copy[index].split(','); row[column] = value
            copy[index] = ','.join(row)
            bad.append('\n'.join(copy))
        for text in bad:
            with self.assertRaises(ValueError):
                read_export(text, 'su3', False)


if __name__ == '__main__':
    unittest.main()
