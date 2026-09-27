import unittest

from plot_hubbard_tutorial import CASES, DATA, POINTS, load_cases, read_export


class HubbardTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_half_filled_gap_and_spin_band(self):
        rows = self.cases['hubbard-half-symmetric'][1]
        self.assertEqual(rows['spinon'][0]['energy'], 0)
        self.assertEqual(rows['spinon'][-1]['energy'], 0)
        self.assertAlmostEqual(max(r['energy'] for r in rows['spinon']), 1.2422814894195624, places=12)
        for branch, momentum in (('holon', -.5), ('antiholon', .5)):
            minimum = min(rows[branch], key=lambda r: r['energy'])
            self.assertEqual(minimum['p_over_pi'], momentum)
            self.assertAlmostEqual(minimum['energy'], .643363511006452197, places=12)

    def test_holon_antiholon_momentum_shift(self):
        rows = self.cases['hubbard-half-symmetric'][1]
        for i, holon in enumerate(rows['holon']):
            antiholon = rows['antiholon'][(i + (POINTS - 1) // 2) % (POINTS - 1)]
            self.assertAlmostEqual(holon['energy'], antiholon['energy'], places=11)

    def test_half_filled_convention_shifts(self):
        sym = self.cases['hubbard-half-symmetric'][1]
        raw = self.cases['hubbard-half-unshifted'][1]
        fermi = self.cases['hubbard-half-fermi'][1]
        for branch in sym:
            for s, r, f in zip(sym[branch], raw[branch], fermi[branch]):
                self.assertAlmostEqual(r['energy'], s['energy'] + 2 * s['delta_n'], places=12)
                self.assertAlmostEqual(f['energy'], s['energy'], places=12)

    def test_energy_columns_and_chemical_potentials(self):
        for name, (metadata, branches) in self.cases.items():
            with self.subTest(case=name):
                mu = float(metadata['Mu symmetric'])
                self.assertAlmostEqual(float(metadata['Mu unshifted']) - mu, 2, places=13)
                for rows in branches.values():
                    for row in rows:
                        self.assertAlmostEqual(row['fermi_energy'],
                                               row['symmetric_energy'] - mu * row['delta_n'], places=12)
                        selected = row['fermi_energy'] if CASES[name][2] == 'fermi' else (
                            row['symmetric_energy'] + (2 * row['delta_n'] if CASES[name][1] == 'unshifted' else 0))
                        self.assertAlmostEqual(row['energy'], selected, places=12)

    def test_doped_reference_and_endpoints(self):
        meta, sym = self.cases['hubbard-doped-symmetric']
        raw = self.cases['hubbard-doped-unshifted'][1]
        self.assertAlmostEqual(float(meta['Mu unshifted']), .37043977496999747, places=11)
        self.assertNotIn('antiholon', sym)
        for branch in sym:
            self.assertEqual(sym[branch][0]['fermi_energy'], 0)
            self.assertEqual(sym[branch][-1]['fermi_energy'], 0)
            for s, r in zip(sym[branch], raw[branch]):
                self.assertAlmostEqual(s['energy'], r['energy'], places=12)

    def test_bad_exports_are_rejected(self):
        name = 'hubbard-doped-symmetric'
        text = (DATA / f'{name}.csv').read_text()
        lines = text.splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith('spinon,'))
        invalid = [text.replace('# Outcome: success', '# Outcome: failure'),
                   text.replace('# Background status: converged', '# Background status: mesh_limit'),
                   text.replace(',converged,', ',precision_limit,', 1),
                   '\n'.join(lines[:index] + lines[index + 1:])]
        for column, value in ((3, ''), (3, 'nan'), (3, 'inf'), (6, '1'), (1, '0.1'),
                              (7, 'nan'), (0, 'antiholon')):
            modified = lines.copy()
            row = modified[index].split(',')
            row[column] = value
            modified[index] = ','.join(row)
            invalid.append('\n'.join(modified))
        for bad in invalid:
            with self.subTest(sample=bad[index:index + 40]), self.assertRaises(ValueError):
                read_export(bad, *CASES[name])
        with self.assertRaises(ValueError):
            read_export(text, .75, 'unshifted', 'fermi')


if __name__ == '__main__':
    unittest.main()
