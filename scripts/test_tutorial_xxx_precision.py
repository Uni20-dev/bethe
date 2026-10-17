import copy
from decimal import Decimal, localcontext
import json
import unittest

from plot_xxx_precision_tutorial import DATA, analyze, load_cases, reference_pi, reference_shift


class XXXPrecisionTutorialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {(r['sites'], r['precision']): r for r in load_cases()}
        cls.cases = [json.loads(line) for line in DATA.read_text().splitlines()]

    def test_independent_reference(self):
        with localcontext() as ctx:
            ctx.prec = 110
            # Analytic values independent of the large-N Bethe calculation.
            self.assertLess(abs(reference_pi()-Decimal(
                '3.1415926535897932384626433832795028841971693993751058209749445923078164062862')),
                Decimal('1e-76'))
            self.assertLess(abs(reference_shift(4)-1), Decimal('1e-105'))
            self.assertLess(abs(reference_shift(6)-Decimal('0.5')), Decimal('1e-105'))
            n = 10**16
            asymptotic = 2*(reference_pi()/n)**2
            self.assertLess(abs(reference_shift(n)/asymptotic-1), Decimal('4e-32'))

    def test_more_than_extra_printed_digits(self):
        # The checked-in binary64 run has genuinely lost the whole energy shift.
        low = self.results[10**8, 'fp64']
        self.assertEqual(low['residual'], 0)
        self.assertEqual(low['shift'], 0)
        self.assertEqual(low['energy_error'], 1)
        self.assertLess(low['root_error'], Decimal('1e-6'))
        self.assertEqual(self.results[10**16, 'fp128']['shift'], 0)
        high = self.results[10**16, 'bits-256']
        self.assertLess(high['energy_error'], Decimal('1e-28'))
        self.assertLess(high['root_error'], Decimal('1e-58'))
        self.assertGreater(high['shift'], 0)

    def test_reject_wrong_state_or_lost_precision(self):
        original = next(c for c in self.cases if c['sites'] == 10**16 and c['precision'] == 'bits-256')
        for mutate in (
            lambda c: c['root'].update(quantum_number=0),
            lambda c: c['state'].update(momentum_index='2'),
            lambda c: c['state'].update(converged=False),
            lambda c: c['state']['energy'].update(value='2500000000000000'),
            lambda c: c['root']['rapidity'].update(value='1'),
            lambda c: c['root']['rapidity'].update(precision=53),
            lambda c: c['metadata'].update({'Working precision (bits)': '53'}),
            lambda c: c['metadata'].pop('Uni20 revision'),
            lambda c: c['summary'].update(Outcome='partial'),
        ):
            with self.subTest(mutate=mutate):
                case = copy.deepcopy(original)
                mutate(case)
                with self.assertRaises(ValueError):
                    analyze(case)


if __name__ == '__main__':
    unittest.main()
