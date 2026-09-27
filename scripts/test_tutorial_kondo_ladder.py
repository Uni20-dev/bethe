import csv
import io
import math
import unittest

from plot_kondo_ladder_tutorial import (A0, DATA, KONDO, KONDO_FIELDS, LADDER, LADDER_FIELDS,
                                       SCHEMAS, ladder_envelope, load_cases, read_kondo, read_ladder)


class KondoLadderTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.cases = load_cases()

    def row(self, name): return self.cases[name][1]

    def test_kondo_independent_high_precision_references(self):
        # Independent mpmath references also recorded in the model guide/native tests.
        for b, m, e in [(.1, .0241204367140108666, -.00120793395564716882),
                        (.5, .112566048096924548, -.029147209883951834),
                        (1, .192111646210937728, -.106792242303035304),
                        (2, .275163637910673047, -.346291700274578678),
                        (10, .387179325881091209, -3.17987143725043218),
                        (1000000, .480548715930107963, -478804.124266610745)]:
            row = self.row(f'kondo-b{b}')
            self.assertLess(abs(row['magnetization']-m), max(row['magnetization_error'], 2e-13))
            self.assertLess(abs(row['energy_change']-e)/b, max(row['scaled_energy_error'], 2e-13))

    def test_kondo_symmetries_scale_and_low_field(self):
        positive, negative, scaled = [self.row(name) for name in ('kondo-b2', 'kondo-negative', 'kondo-scaled')]
        self.assertEqual(positive['energy_change'], negative['energy_change'])
        self.assertEqual(positive['magnetization'], -negative['magnetization'])
        self.assertAlmostEqual(scaled['energy_change'], 3*positive['energy_change'], places=12)
        self.assertAlmostEqual(scaled['magnetization'], positive['magnetization'], places=12)
        self.assertAlmostEqual(scaled['zero_field_susceptibility'], A0/3, places=14)
        weak = self.row('kondo-b0.0001')
        self.assertAlmostEqual(weak['magnetization']/.0001, A0, delta=1e-9)
        self.assertAlmostEqual(weak['energy_change']/(.0001**2), -A0/2, delta=1e-9)
        self.assertEqual(self.row('kondo-b0')['energy_change'], 0)

    def test_kondo_energy_derivative_crossover_and_saturation(self):
        for b in (1, 2):
            derivative = -(self.row(f'kondo-deriv{b}-plus')['energy_change']-
                           self.row(f'kondo-deriv{b}-minus')['energy_change'])/.002
            self.assertAlmostEqual(derivative, self.row(f'kondo-b{b}')['magnetization'], delta=3e-8)
        values = [self.row(f'kondo-b{b}')['magnetization'] for b in KONDO_FIELDS]
        self.assertTrue(all(a < b for a, b in zip(values, values[1:])))
        self.assertTrue(.47 < values[-1] < .5)
        self.assertTrue(-.5 < self.row('kondo-b1000000')['energy_change']/1e6 < -.47)
        self.assertGreater(self.row('kondo-b1.01')['evaluations'], 0)

    def test_ladder_sector_envelope_and_crossing_ties(self):
        fixed = [self.row(f'ladder-m{m}')['energy'] for m in range(7)]
        crossings = [b-a for a, b in zip(fixed, fixed[1:])]
        self.assertEqual(crossings, sorted(crossings))
        self.assertEqual(crossings[0], 1)
        self.assertEqual(crossings[-1], 9)
        for h in LADDER_FIELDS:
            row = self.row(f'ladder-h{h}')
            best, _ = ladder_envelope(self.cases, h)
            self.assertAlmostEqual(row['energy'], best, places=12)
            self.assertAlmostEqual(fixed[row['magnetization']]-h*row['magnetization'], best, places=12)
        self.assertEqual(self.row('ladder-h8.5')['magnetization'], 5)
        self.assertEqual(self.row('ladder-h0')['energy'], -18)
        self.assertEqual(self.row('ladder-h10')['energy'], -48)
        # At h=9 the envelope and native analytic path may choose different tied states.
        self.assertAlmostEqual(fixed[5]-9*5, fixed[6]-9*6, places=12)

    def test_ladder_descendants_fields_and_product_limits(self):
        d = self.row('ladder-descendant')
        self.assertEqual(d['populations'], [4, 1, 1, 0])
        self.assertEqual(d['shape'], [4, 2, 0, 0])
        self.assertAlmostEqual(d['permutation_energy'], 1-math.sqrt(5), places=12)
        self.assertEqual(self.row('ladder-desc-plus')['magnetization'], 2)
        self.assertEqual(self.row('ladder-desc-minus')['magnetization'], -2)
        self.assertEqual(self.row('ladder-desc-plus')['energy'], self.row('ladder-desc-minus')['energy'])
        self.assertAlmostEqual(self.row('ladder-desc-rung')['energy'], d['energy']-2.5, places=12)
        self.assertEqual(self.row('ladder-desc-rung')['levels'], d['levels'])
        shifted, original = self.row('ladder-fixed-field'), self.row('ladder-m2')
        self.assertAlmostEqual(shifted['energy'], original['energy']-6, places=12)
        self.assertEqual(shifted['levels'], original['levels'])
        self.assertEqual(shifted['energy'], self.row('ladder-reversed')['energy'])
        self.assertAlmostEqual(self.row('ladder-two')['energy'], -2.5, places=12)
        self.assertAlmostEqual(self.row('ladder-four')['energy'], -5, places=12)
        self.assertEqual(list(map(len, self.row('ladder-four')['levels'])), [3, 2, 1])
        for name in ('ladder-triplet-even', 'ladder-triplet-odd'):
            n, rung, *_ = LADDER[name]
            product = 3*n*(1-rung)/4
            self.assertAlmostEqual(self.row(name)['energy']-product,
                                   rung-4*math.sin(math.pi*(n//2)/n)**2, places=12)

    def test_corrupt_units_roots_and_status_are_rejected(self):
        for name, edits in {
                'kondo-b2': [('response', 'magnetization', ''), ('response', 'scaled_energy_error', '1'),
                             ('response', 'zero_field_susceptibility', '1')],
                'ladder-descendant': [('states', 'descendant', 'false'), ('states', 'magnetization', '2'),
                                      ('roots', 'level', '2'), ('roots', 'rapidity', '42'),
                                      ('representations', 'highest_weight', '3')]} .items():
            tables = ['response'] if name in KONDO else SCHEMAS
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in tables}
            def check(bad):
                return read_kondo(bad['response'], KONDO[name]) if name in KONDO else read_ladder(bad, LADDER[name])
            for table, key, value in edits:
                lines = texts[table].splitlines()
                reader = csv.DictReader(line for line in lines if not line.startswith('#'))
                rows = list(reader); rows[0][key] = value
                out = io.StringIO(); writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                writer.writeheader(); writer.writerows(rows)
                bad = dict(texts)
                bad[table] = '\n'.join(line for line in lines if line.startswith('#'))+'\n'+out.getvalue()
                with self.subTest(name=name, column=key):
                    with self.assertRaises(ValueError): check(bad)
            table = next(iter(tables)); bad = dict(texts)
            bad[table] = texts[table].replace('# Outcome: success', '# Outcome: partial')
            with self.assertRaises(ValueError): check(bad)


if __name__ == '__main__': unittest.main()
