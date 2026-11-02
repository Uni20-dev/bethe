# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from spectral_tables import RealValue, group_rows, index_rows, read_tables


def column(name, kind='string', **extra):
    return dict(id=name, label=name, unit='', description='', type=kind, nullable=False, **extra)


def fixture():
    columns = [column('candidate_id'), column('run', 'uint64', encoding='decimal_string'),
               column('accepted', 'bool'), column('root_index', 'int32', encoding='number'),
               column('lambda_real', 'real', radix=2, precision='per_value', encoding='value_precision'),
               dict(column('deviation', 'real', radix=2, precision_bits=113,
                           encoding='decimal_string'), nullable=True),
               column('spin', 'half_int', storage_bits=64, encoding='number')]
    rows = [['state-a', '9007199254740993', True, 0,
             {'value': '0.123456789012345678901234567890123456789', 'precision': 160}, '1e-400', 0],
            ['state-a', '9007199254740993', True, 1, {'value': '1/3', 'precision': 'exact'}, '0', 0],
            ['failed-state', '9007199254740993', False, 0,
             {'value': 'nan', 'precision': 160}, None, 0],
            ['state-a', '9007199254740993', True, 2,
             {'value': '-0', 'precision': 256}, '2e-300', 0]]
    table = dict(title='Roots including failures', first_row=7, columns=columns, rows=rows,
                 metadata={'Family': 'selected roots', 'Channel': 'zz', 'Precision': 'mixed'},
                 summary={'Outcome': 'partial', 'Failed': '1'})
    return {'tables': {'roots': table}, 'status': 'complete'}


class Tables(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name)/'data.json'

    def read(self, document, **kwargs):
        self.path.write_text(json.dumps(document), encoding='utf-8')
        return read_tables(self.path, **kwargs)

    def test_precision_failures_and_identities_are_retained(self):
        result = self.read(fixture(), expected_metadata={'roots': {'Channel': 'zz'}})
        table = result.tables['roots']
        self.assertEqual(result.envelope['status'], 'complete')
        self.assertEqual(table.summary['Outcome'], 'partial')
        self.assertEqual(table.first_row, 7)
        self.assertEqual(len(table.rows), 4)
        a, rational, failed, negative_zero = table.rows
        self.assertEqual(a['run'], 9007199254740993)
        self.assertEqual(a['lambda_real'].precision, 160)
        with localcontext() as ctx:
            ctx.prec = 6
            self.assertEqual(str(a['lambda_real'].decimal()), '0.123456789012345678901234567890123456789')
        self.assertEqual(a['deviation'].precision, 113)
        self.assertEqual(a['deviation'].decimal(), Decimal('1e-400'))
        with self.assertRaises(ValueError):
            a['deviation'].float64()
        self.assertEqual(rational['lambda_real'].fraction(), Fraction(1, 3))
        with self.assertRaises(ValueError):
            rational['lambda_real'].decimal()
        self.assertEqual(rational['lambda_real'].float64(), 1/3)
        self.assertFalse(failed['accepted'])
        self.assertIsNone(failed['deviation'])
        self.assertTrue(failed['lambda_real'].decimal().is_nan())
        with self.assertRaises(ValueError):
            failed['lambda_real'].finite_decimal()
        self.assertTrue(negative_zero['lambda_real'].decimal().is_signed())
        self.assertEqual(negative_zero['lambda_real'].precision, 256)
        self.assertEqual(len(group_rows(table, ['candidate_id'])[('state-a',)]), 3)
        self.assertEqual(len(index_rows(table, ['candidate_id', 'root_index'])), 4)
        with self.assertRaises(ValueError):
            index_rows(table, ['candidate_id'])
        for keys in ([], ['absent'], ['candidate_id', 'candidate_id'], 'candidate_id', ['lambda_real']):
            with self.subTest(keys=keys), self.assertRaises(ValueError):
                group_rows(table, keys)
        with self.assertRaises(TypeError):
            a['run'] = 2
        with self.assertRaises(TypeError):
            table.metadata['Channel'] = 'raising'

    def test_checksum_required_tables_and_metadata(self):
        result = self.read(fixture())
        self.assertEqual(result.sha256, hashlib.sha256(self.path.read_bytes()).hexdigest())
        self.assertEqual(read_tables(self.path, sha256=result.sha256).source, str(self.path.resolve()))
        for kwargs in ({'sha256': '0'*64}, {'table_name': 'roots'},
                       {'expected_metadata': {'missing': {}}},
                       {'expected_metadata': {'roots': {'Channel': 'raising'}}}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                read_tables(self.path, **kwargs)

    def test_single_snapshot_and_exact_half_integers(self):
        single = fixture()['tables']['roots']
        single.pop('summary')
        payload = json.dumps(single).replace('"1e-400", 0]', '"1e-400", 4503599627370496.5]')
        self.path.write_text(payload)
        table = read_tables(self.path, table_name='roots').tables['roots']
        self.assertIsNone(table.summary)
        self.assertEqual(table.rows[0]['spin'], Decimal('4503599627370496.5'))
        with self.assertRaises(ValueError):
            read_tables(self.path)

    def test_bad_schema_and_cells_are_errors(self):
        mutations = [lambda t: t['columns'].append(t['columns'][0]),
                     lambda t: t['rows'][0].pop(),
                     lambda t: t.update(first_row=-1),
                     lambda t: t['metadata'].update(Channel=2),
                     lambda t: t['columns'][4].update(radix=10),
                     lambda t: t['columns'][4].update(precision_bits=160),
                     lambda t: t['columns'][4].update(encoding='number'),
                     lambda t: t['rows'][0].__setitem__(0, None),
                     lambda t: t['rows'][0].__setitem__(1, 9007199254740993),
                     lambda t: t['rows'][0].__setitem__(1, str(2**64)),
                     lambda t: t['rows'][0].__setitem__(2, 1),
                     lambda t: t['rows'][0].__setitem__(3, True),
                     lambda t: t['rows'][0].__setitem__(4, {'value': '1', 'precision': True}),
                     lambda t: t['rows'][0].__setitem__(4, {'value': '1/3', 'precision': 160}),
                     lambda t: t['rows'][0].__setitem__(4, {'value': '1/0', 'precision': 'exact'}),
                     lambda t: t['rows'][0].__setitem__(4, {'value': 'nan', 'precision': 'exact'}),
                     lambda t: t['rows'][0].__setitem__(4, {'value': '1', 'precision': 0}),
                     lambda t: t['rows'][0].__setitem__(5, 1.2),
                     lambda t: t['rows'][0].__setitem__(6, .25)]
        for i, mutate in enumerate(mutations):
            with self.subTest(mutation=i), self.assertRaises(ValueError):
                document = fixture()
                mutate(document['tables']['roots'])
                self.read(document)

    def test_duplicate_keys_nonstandard_json_and_conversion_range(self):
        for text in ('{"tables":{},"tables":{}}', '{"x":NaN}', '{"x":Infinity}'):
            self.path.write_text(text)
            with self.assertRaises(ValueError):
                read_tables(self.path)
        for text in ('1e400', 'inf', '-nan', '1e-400'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                RealValue(text, 160).float64()
        self.assertEqual(RealValue('-0', 160).float64(), 0)
        self.assertEqual(RealValue('5e-324', 160).float64(), 5e-324)


if __name__ == '__main__':
    unittest.main()
