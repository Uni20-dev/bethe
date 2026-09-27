import itertools
import math
import unittest

from plot_q_boson_tutorial import CASES, COUNT, DATA, L, N, TABLES, load_cases, read_exports


class QBosonTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = load_cases()

    def test_free_and_phase_roots(self):
        for name in ('q-boson-free', 'q-boson-phase'):
            _, states, roots = self.cases[name]
            for state_id, rr in roots.items():
                lifted_p = 2*math.pi/L*sum(r['m'] for r in rr)
                for r in rr:
                    k = 2*math.pi*r['m']/L if name.endswith('free') else (2*math.pi*r['I']+lifted_p)/(L+N)
                    self.assertAlmostEqual(r['k'], k, places=13)
            self.assertEqual(len(states), COUNT)
        self.assertEqual(float(self.cases['q-boson-free'][0]['Ground energy']), 0)
        self.assertAlmostEqual(float(self.cases['q-boson-phase'][0]['Ground energy']), 4-2*math.sqrt(2), places=13)
        self.assertGreater(float(self.cases['q-boson-eta05'][0]['Ground energy']), 0)
        self.assertLess(float(self.cases['q-boson-eta05'][0]['Ground energy']), 4-2*math.sqrt(2))
        self.assertEqual(COUNT, 35)
        self.assertEqual(math.comb(L, N), 10)  # Hard-core Hilbert space is different.

    def test_occupation_basis_spectral_moments(self):
        # Trace(H) and Trace(H^2) directly from occupation-basis hopping,
        # independent of the Bethe roots and their label conventions.
        occupations = [x for x in itertools.product(range(N+1), repeat=L) if sum(x) == N]
        self.assertEqual(len(occupations), COUNT)
        for name, eta in CASES.items():
            def qnumber(n):
                if n == 0: return 0
                if eta is None: return 1
                if eta == '0': return n
                return math.expm1(-2*float(eta)*n)/math.expm1(-2*float(eta))
            trace2 = COUNT*(2*N)**2
            for occupation in occupations:
                for source in range(L):
                    for target in ((source-1) % L, (source+1) % L):
                        trace2 += qnumber(occupation[source])*qnumber(occupation[target]+1)
            energies = [r['energy'] for r in self.cases[name][1]]
            self.assertAlmostEqual(sum(energies), COUNT*2*N, places=10)
            self.assertAlmostEqual(sum(e*e for e in energies), trace2, places=10)

    def test_bad_counts_joins_and_reference_rejected(self):
        texts = {table: (DATA / f'q-boson-eta05-{table}.csv').read_text() for table in TABLES}
        bad = []
        for table, old, new in [('states', '# Outcome: success', '# Outcome: partial'),
                                ('states', '# Converged count: 35', '# Converged count: 34'),
                                ('states', ',true,converged', ',false,failed'),
                                ('roots', '# Energy reference: +2N', '# Energy reference: -2N')]:
            changed = texts.copy(); changed[table] = changed[table].replace(old, new, 1); bad.append(changed)
        for table, column, value in [('roots', 0, '999'), ('roots', 1, '2'), ('roots', 2, '0'),
                                     ('roots', 3, '5'), ('roots', 4, 'nan'), ('states', 3, '1')]:
            changed = texts.copy(); lines = changed[table].splitlines()
            i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
            row = lines[i].split(','); row[column] = value; lines[i] = ','.join(row)
            changed[table] = '\n'.join(lines); bad.append(changed)
        changed = texts.copy(); lines = changed['roots'].splitlines()
        i = next(i for i, line in enumerate(lines) if line.startswith('0,'))
        changed['roots'] = '\n'.join(lines[:i]+lines[i+1:]); bad.append(changed)
        for candidate in bad:
            with self.assertRaises(ValueError):
                read_exports(candidate, '0.5')


if __name__ == '__main__':
    unittest.main()
