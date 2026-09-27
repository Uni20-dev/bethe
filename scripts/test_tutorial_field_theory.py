import math
import unittest

from plot_field_theory_tutorial import CASES, DATA, SCHEMAS, SG_LENGTHS, LY_LENGTHS, LY_EXCITED_LENGTHS, load_cases, read_exports


class FieldTheoryTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_particle_masses_threshold_boundary_and_charges(self):
        branches = self.cases['sg-particles-p0.4'][1]['branches']
        self.assertAlmostEqual(branches['B1'][0]['energy'], 1.1755705045849463, places=14)
        self.assertAlmostEqual(branches['B2'][0]['energy'], 1.9021130325903071, places=14)
        self.assertLess(branches['B2'][0]['energy'], branches['soliton+antisoliton'][0]['energy'])
        self.assertEqual(branches['soliton+antisoliton'][0]['charge'], 0)
        self.assertEqual(branches['soliton+soliton'][0]['charge'], 2)
        for p, particles in [(.4, ['soliton', 'antisoliton', 'B1', 'B2']),
                             (.5, ['soliton', 'antisoliton', 'B1']), (1, ['soliton', 'antisoliton'])]:
            b = self.cases[f'sg-particles-p{p}'][1]['branches']
            self.assertEqual([k for k in b if '+' not in k], particles)
        # Unequal masses share velocity at the continuum minimum, not momentum.
        a, b, k = 1, branches['B1'][0]['rest_energy'], 2.5
        common_velocity = math.hypot(a, k*a/(a+b))+math.hypot(b, k*b/(a+b))
        threshold = branches['soliton+B1'][32]['energy']
        self.assertAlmostEqual(threshold, common_velocity, places=14)
        self.assertGreater(math.hypot(a, k/2)+math.hypot(b, k/2), threshold)

    def test_mass_length_rescaling(self):
        for first, second, table, keys in [
                ('sg-particles-p0.4', 'sg-particles-scaled', 'dispersion', ['rest_energy', 'momentum', 'energy']),
                ('sg-nlie-l1', 'sg-scaled-nlie', 'levels', ['casimir_energy']),
                ('ly-vacuum-l1', 'ly-vacuum-scaled', 'vacuum', ['casimir_energy']),
                ('ly-excited-l5', 'ly-excited-scaled', 'levels', ['casimir_energy'])]:
            a, b = self.cases[first][1][table], self.cases[second][1][table]
            self.assertEqual(len(a), len(b))
            for x, y in zip(a, b):
                for key in keys:
                    self.assertAlmostEqual(2*x[key], y[key], places=11)
                if 'scaling_function' in x:
                    self.assertAlmostEqual(x['scaling_function'], y['scaling_function'], places=11)

    def test_free_fermion_gap_and_interacting_oracle(self):
        exact = 2*math.hypot(1, math.pi)
        self.assertAlmostEqual(self.cases['sg-free-nlie'][1]['gap'][0]['gap'], exact, places=11)
        self.assertAlmostEqual(self.cases['sg-free-by'][1]['levels'][0]['energy'], exact, places=11)
        # Independent p=2 Gauss-grid oracle documented in sine-gordon-excited.md.
        self.assertAlmostEqual(self.cases['sg-nlie-l1'][1]['levels'][1]['casimir_energy'], 5.11055662981005, places=8)
        gaps = [self.cases[f'sg-nlie-l{l}'][1]['gap'][0] for l in SG_LENGTHS]
        self.assertLess(abs(gaps[0]['scaled_gap']-.75), .006)
        difference = [g['gap']-self.cases[f'sg-by-l{l}'][1]['levels'][0]['energy']
                      for l, g in zip(SG_LENGTHS, gaps)]
        self.assertTrue(all(x > 0 for x in difference))
        self.assertGreater(difference[3], .13)
        self.assertLess(difference[-1], 1e-5)
        self.assertGreater(difference[-1], 1e-6)  # Detect accidentally plotting BY twice.

    def test_lee_yang_uv_ir_and_vacuum_oracle(self):
        values = [self.cases[f'ly-vacuum-l{l}'][1]['vacuum'][0] for l in LY_LENGTHS]
        charges = [r['effective_central_charge'] for r in values]
        self.assertLess(abs(charges[0]-.4), 3e-7)  # Finite r=0.001, not the exact UV limit.
        self.assertTrue(all(a > b > 0 for a, b in zip(charges, charges[1:])))
        self.assertAlmostEqual(self.cases['ly-vacuum-l1'][1]['vacuum'][0]['scaling_function'], -.15320688011013006, places=12)
        self.assertAlmostEqual(.4+24*(-.2), -4.4, places=14)
        gaps = [self.cases[f'ly-excited-l{l}'][1]['gap'][0] for l in LY_EXCITED_LENGTHS]
        self.assertTrue(all(a['gap'] > b['gap'] > 1 for a, b in zip(gaps, gaps[1:])))
        for g in gaps:
            self.assertGreater(g['gap']-1, 1000*g['gap_error'])
        self.assertLess(gaps[-1]['gap']-1, 1e-7)
        # Gap must not be replaced by the one-particle Casimir energy.
        a = self.cases['ly-excited-l5'][1]
        self.assertGreater(a['gap'][0]['gap']-a['levels'][1]['casimir_energy'], .001)

    def assert_corruptions_rejected(self, name, edits):
        texts = {table: (DATA/f'{name}-{table}.csv').read_text() for table in SCHEMAS[CASES[name][0]]}
        for table, old, new in edits:
            with self.subTest(name=name, old=old):
                self.assertIn(old, texts[table])
                corrupt = dict(texts)
                corrupt[table] = corrupt[table].replace(old, new, 1)
                with self.assertRaises(ValueError):
                    read_exports(corrupt, CASES[name])

    def test_failed_or_inconsistent_finite_volume_tables_rejected(self):
        self.assert_corruptions_rejected('sg-nlie-l1', [
            ('gap', '# Outcome: success', '# Outcome: partial'),
            ('gap', '5.4440185956186031', '5.110556629808821'),
            ('source', '0.5,2,', '1.5,2,'),
            ('levels', 'true,converged', 'false,converged'),
            ('levels', 'two_soliton,', 'vacuum,'),
            ('source', '# Circumference: 1', '# Circumference: 2')])
        self.assert_corruptions_rejected('ly-excited-l5', [
            ('gap', '1.0306379212700283', 'nan'),
            ('gap', '1.0177227524048455e-12', '0'),
            ('source', '0.024923220090329563', '0.1'),
            ('source', '# Status: converged', '# Status: failed')])

    def test_dispersion_and_asymptotic_corruption_rejected(self):
        self.assert_corruptions_rejected('sg-particles-p0.4', [
            ('dispersion', 'soliton,particle,1,', 'soliton,particle,0,'),
            ('dispersion', 'B2,particle,0,', 'B3,particle,0,'),
            ('dispersion', '# Points per branch: 65', '# Points per branch: 64')])
        self.assert_corruptions_rejected('sg-by-l10', [
            ('levels', '0.5,2,', '0.5,0,'),
            ('levels', 'converged,converged', 'quadrature_limit,converged')])


if __name__ == '__main__':
    unittest.main()
