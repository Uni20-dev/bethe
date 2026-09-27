from collections import Counter
import itertools
import math
import unittest

from plot_biquadratic_tutorial import CASES, DATA, load_cases, multiplicity, read_exports


class BiquadraticTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_complete_four_site_levels_and_physical_spins(self):
        expected = {0: [-(15+math.sqrt(17))/2, -(15-math.sqrt(17))/2],
                    2: [-6-math.sqrt(2), -6, -6+math.sqrt(2)], 4: [-3]}
        magnetizations = Counter(sum(word) for word in itertools.product((-1, 0, 1), repeat=4))
        for sign in ('af', 'ferro'):
            physical_states = 0
            spin_counts = Counter()
            for ell in (0, 2, 4):
                _, rows, content = self.cases[f'bq-{sign}-q-l{ell}']
                energies = sorted((-1 if sign == 'ferro' else 1)*e for e in expected[ell])
                for i, (row, energy) in enumerate(zip(rows, energies)):
                    self.assertAlmostEqual(row['energy'], energy, places=12)
                    physical_states += row['multiplicity']
                    spin_counts.update(content[i])  # Exclude duplicate ground-reference rows.
            self.assertEqual(physical_states, 3**4)
            for spin in range(5):
                self.assertEqual(spin_counts[spin], magnetizations[spin]-magnetizations[spin+1])
        self.assertEqual(self.cases['bq-af-q-l2'][2][0], {1: 1, 2: 1})

    def test_sign_reversal_and_real_root_omission(self):
        for ell in (0, 2, 4):
            af = self.cases[f'bq-af-q-l{ell}'][1]
            ferro = self.cases[f'bq-ferro-q-l{ell}'][1]
            for a, f in zip(af, reversed(ferro)):
                self.assertAlmostEqual(a['energy'], -f['energy'], places=12)
                self.assertAlmostEqual(a['reference_energy'], f['reference_energy'], places=12)
                self.assertAlmostEqual(a['gap']+f['gap'], (9+math.sqrt(17))/2, places=12)
        for sign in ('af', 'ferro'):
            real = self.cases[f'bq-{sign}-real-l0'][1]
            full = self.cases[f'bq-{sign}-q-l0'][1]
            self.assertEqual(len(real), 1)
            self.assertEqual(len(full), 2)
            self.assertAlmostEqual(real[0]['energy'], full[0 if sign == 'af' else 1]['energy'], places=12)
        self.assertAlmostEqual(self.cases['bq-af-singlet-n4'][1][0]['gap'], math.sqrt(17), places=12)

    def test_band_and_droplet_edges(self):
        _, band, _ = self.cases['bq-ferro-band-n32']
        for row in band:
            self.assertAlmostEqual(row['wave_number'], math.pi*row['mode']/32, places=14)
            self.assertAlmostEqual(row['gap'], 3+2*math.cos(row['wave_number']), places=13)
        self.assertAlmostEqual(band[0]['gap'], 3-2*math.cos(math.pi/32), places=13)
        for kind, limit in [('pair', 5/3), ('triple', 2)]:
            gaps = [self.cases[f'bq-ferro-{kind}-n{n}'][1][0]['gap'] for n in (16, 32, 64, 128)]
            self.assertTrue(all(a > b for a, b in zip(gaps, gaps[1:])))
            self.assertTrue(all(g > limit for g in gaps))
            self.assertLess(gaps[-1]-limit, .001)
        self.assertIsNone(self.cases['bq-ferro-pair-n64'][1][0]['multiplicity'])
        self.assertGreater(multiplicity(60), 2**64-1)
        self.assertLess(self.cases['bq-ferro-pair-n128'][1][0]['gap'], 2)
        self.assertLess(self.cases['bq-ferro-triple-n128'][1][0]['gap'], 8/3)

    def test_invalid_q_and_spin_tables_rejected(self):
        name = 'bq-af-q-l0'
        texts = {table: (DATA / f'{name}-{table}.csv').read_text() for table in ('states', 'spin_content')}
        bad = []
        for table, old, new in [('states', '# Outcome: success', '# Outcome: partial'),
                                ('states', '# Discovered levels: 2', '# Discovered levels: 1'),
                                ('states', ',true,converged;', ',false,converged;'),
                                ('spin_content', ',exact', ',unavailable')]:
            changed = texts.copy(); changed[table] = changed[table].replace(old, new, 1); bad.append(changed)
        for table, column, value in [('states', 2, '8'), ('states', 4, '1'), ('states', 6, 'nan'),
                                     ('states', 8, '1'), ('spin_content', 2, '1'), ('spin_content', 3, '2')]:
            changed = texts.copy(); lines = changed[table].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            row = lines[i].split(','); row[column] = value; lines[i] = ','.join(row)
            changed[table] = '\n'.join(lines); bad.append(changed)
        for table in texts:
            changed = texts.copy(); lines = changed[table].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            changed[table] = '\n'.join(lines[:i]+lines[i+1:]); bad.append(changed)
        for candidate in bad:
            with self.assertRaises(ValueError): read_exports(candidate, CASES[name])

    def test_failed_reference_and_droplet_rejected(self):
        name = 'bq-af-singlet-n16'
        texts = {table: (DATA / f'{name}-{table}.csv').read_text() for table in ('states', 'reference')}
        texts['reference'] = texts['reference'].replace(',true,converged', ',false,failed')
        with self.assertRaises(ValueError): read_exports(texts, CASES[name])
        name = 'bq-ferro-pair-n64'
        text = (DATA / f'{name}-states.csv').read_text()
        for old, new in [('# Outcome: success', '# Outcome: partial'),
                         ('# Selected modes: 1', '# Selected modes: 2'),
                         (',true,converged', ',false,stalled')]:
            with self.assertRaises(ValueError):
                read_exports({'states': text.replace(old, new, 1)}, CASES[name])


if __name__ == '__main__':
    unittest.main()
