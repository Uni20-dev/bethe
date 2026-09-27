import csv
import io
import math
import unittest

from plot_potts_qg_tutorial import CASES, DATA, EVEN, ODD, POTTS_SIZES, casimir, load_cases, read_exports, schema


class PottsQuantumGroupTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.cases = load_cases()

    def test_potts_exact_small_ring_and_selected_count(self):
        rows = self.cases['potts-l2'][1]['levels']
        self.assertAlmostEqual(rows[0]['energy'], -2-2*math.sqrt(3), places=12)
        self.assertAlmostEqual(rows[1]['energy'], -(1+math.sqrt(57))/2, places=12)
        for n in POTTS_SIZES:
            rows = self.cases[f'potts-l{n}'][1]['levels']
            self.assertEqual(len(rows), 2*n+1)
            self.assertLess(len(rows), 3**n)
            self.assertEqual(sum(int(r['k']) == n//2 for r in rows), 2)

    def test_potts_cft_and_thermodynamic_branch(self):
        cs = [casimir(self.cases[f'potts-l{n}']) for n in POTTS_SIZES]
        self.assertTrue(all(a > b > .8 for a, b in zip(cs, cs[1:])))
        self.assertLess(cs[-1]-.8, .0001)
        rows = [r for r in self.cases['potts-l64'][1]['levels'] if r['charge'] == '1']
        self.assertLess(abs(rows[0]['x_scaled']-2/15), .0003)
        self.assertLess(abs(rows[1]['x_scaled']-17/15), .001)
        self.assertLess(abs(rows[32]['gap']-3*math.sqrt(3)), .04)

    def test_endpoint_parity_minima_and_signed_sectors(self):
        for n in (*EVEN, *ODD):
            meta, tables = self.cases[f'qg-free-n{n}']
            exact = (1-(1/math.tan(math.pi/(2*n)) if n % 2 == 0 else 1/math.sin(math.pi/(2*n))))/2
            self.assertAlmostEqual(min(r['energy'] for r in tables['blocks']), exact, places=12)
            self.assertEqual(int(meta['Size-two blocks']), math.comb(n-2, n//2-1) if n % 2 == 0 else 0)
        self.assertLess(abs(casimir(self.cases['qg-free-n12'])+2), .003)
        self.assertLess(abs(casimir(self.cases['qg-free-n13'])-1), .002)
        spectra = []
        for name in ('qg-free-plus', 'qg-free-minus'):
            spectra.append(sorted(round(r['energy'], 12) for r in self.cases[name][1]['blocks']
                                  for _ in range(r['block_size'])))
        self.assertEqual(*spectra)
        for name in ('qg-free-up', 'qg-free-down'):
            rows = self.cases[name][1]['blocks']
            self.assertEqual([(r['energy'], r['block_size']) for r in rows], [(0, 1)])

    def test_two_site_defect_not_two_eigenvectors(self):
        rows = self.cases['qg-free-n2'][1]['blocks']
        self.assertEqual([(r['energy'], r['block_size']) for r in rows], [(0, 2)])
        # Independent two-site spin Hamiltonian: rank one, nonzero, square zero.
        h = [[-.5j, .5], [.5, .5j]]
        self.assertNotEqual(h[0][0], 0)
        self.assertEqual(h[0][0]*h[1][1]-h[0][1]*h[1][0], 0)
        self.assertTrue(all(sum(h[i][k]*h[k][j] for k in range(2)) == 0
                            for i in range(2) for j in range(2)))

    def test_regular_family_is_not_complete_and_one_root_is_analytic(self):
        rows = self.cases['qg-scan'][1]['levels']
        self.assertEqual(len(rows), 10)
        self.assertLess(len(rows), math.comb(8, 2))
        self.assertEqual(rows[0]['gap_from_sea'], 0)
        rows = self.cases['qg-one'][1]['levels']
        self.assertEqual(len(rows), 5)
        for row in rows:
            i = int(row['numbers'])
            self.assertAlmostEqual(row['energy'], 7*.6/4-.6-math.cos(math.pi*i/8), places=12)

    def test_corrupt_observables_joins_and_blocks_are_rejected(self):
        mutations = {
            'potts-l4': [('levels', 'gap', 'nan'), ('levels', 'x_scaled', '1'),
                         ('levels', 'k', '4'), ('levels', 'residual', '1')],
            'qg-free-n6': [('blocks', 'block_size', '2'), ('blocks', 'zero_occupation', '1'),
                           ('blocks', 'modes', '1,2,3'), ('blocks', 'energy', '')],
            'qg-scan': [('roots', 'source', 'failed'), ('roots', 'lambda', '2'),
                        ('levels', 'numbers', '1,6'), ('reference', 'energy_shift', '0')],
        }
        for name, edits in mutations.items():
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schema(CASES[name])}
            for table, key, value in edits:
                lines = texts[table].splitlines()
                reader = csv.DictReader(line for line in lines if not line.startswith('#'))
                rows = list(reader); rows[0][key] = value
                out = io.StringIO(); writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                writer.writeheader(); writer.writerows(rows)
                corrupt = dict(texts)
                corrupt[table] = '\n'.join(line for line in lines if line.startswith('#'))+'\n'+out.getvalue()
                with self.subTest(name=name, column=key):
                    with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            table = next(iter(texts)); corrupt = dict(texts)
            corrupt[table] = texts[table].replace('# Outcome: success', '# Outcome: partial')
            with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])
            # A missing row is not a valid smaller family or smaller sector.
            lines = texts[table].splitlines()
            i = next(i for i, line in enumerate(lines) if not line.startswith('#'))+1
            corrupt = dict(texts); corrupt[table] = '\n'.join(lines[:i]+lines[i+1:])+'\n'
            with self.assertRaises(ValueError): read_exports(corrupt, CASES[name])


if __name__ == '__main__': unittest.main()
