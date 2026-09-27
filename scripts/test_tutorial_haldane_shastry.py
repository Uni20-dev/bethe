from collections import Counter
import math
import unittest

from plot_haldane_shastry_tutorial import DATA, load_cases, read_exports


class HaldaneShastryTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_ground_states_and_normalization(self):
        for n, (rows, content) in self.cases.items():
            ground = [r for r in rows if r['gap'] == 0]
            expected = -math.pi**2*(n*n + (5 if n % 2 == 0 else -1))/(24*n)
            self.assertEqual(len(ground), 1 if n % 2 == 0 else 2)
            self.assertEqual(sum(r['degeneracy'] for r in ground), 1 if n % 2 == 0 else 4)
            for r in ground:
                self.assertAlmostEqual(r['energy'], expected, places=12)
            polarized = next(r for r in rows if not r['motif'])
            self.assertAlmostEqual(polarized['energy'], math.pi**2*(n*n-1)/(24*n), places=12)
            self.assertEqual(polarized['degeneracy'], n+1)
            for r in rows:
                self.assertAlmostEqual(r['energy'], polarized['energy']-(math.pi/n)**2*sum(m*(n-m) for m in r['motif']), places=12)

    def test_full_spin_word_count(self):
        for n, (_, content) in self.cases.items():
            totals = sum(content.values(), Counter())
            # Independent SU(2) count from differences of spin-projection sectors.
            for down in range(n//2+1):
                s = n/2-down
                expected = math.comb(n, down)-(math.comb(n, down-1) if down else 0)
                self.assertEqual(totals[s], expected)
            self.assertEqual(sum(int(2*s+1)*m for s, m in totals.items()), 2**n)
        rows, content = self.cases[8]
        i = next(i for i, r in enumerate(rows) if r['motif'] == (2, 4, 6))
        self.assertEqual(content[i], {0: 1, 1: 1})
        self.assertEqual(rows[i]['degeneracy'], 4)

    def test_partial_and_corrupt_exports_rejected(self):
        texts = {name: (DATA / f'hs-n8-{name}.csv').read_text() for name in ('levels', 'spin_content')}
        bad = []
        for table, old, new in [('levels', '# Outcome: success', '# Outcome: partial'),
                                ('levels', '# Motifs enumerated: 34', '# Motifs enumerated: 33'),
                                ('spin_content', ',true,complete', ',false,work_limit')]:
            changed = texts.copy(); changed[table] = changed[table].replace(old, new, 1); bad.append(changed)
        for table, column, value in [('levels', 1, '1 2'), ('levels', 4, '1'), ('levels', 8, '2'),
                                     ('spin_content', 0, '999'), ('spin_content', 1, '.5'),
                                     ('spin_content', 2, '2')]:
            changed = texts.copy(); lines = changed[table].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            row = lines[i].split(','); row[column] = value; lines[i] = ','.join(row)
            changed[table] = '\n'.join(lines); bad.append(changed)
        for table in texts:
            changed = texts.copy(); lines = changed[table].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            changed[table] = '\n'.join(lines[:i]+lines[i+1:]); bad.append(changed)
        for changed in bad:
            with self.assertRaises(ValueError):
                read_exports(changed, 8)


if __name__ == '__main__':
    unittest.main()
