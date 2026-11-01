"""Cross-sector source permission preserves physics, source records and gaps."""
from copy import deepcopy
from decimal import Decimal
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import unittest

import plot_xxx_spinon_slice_tutorial as tutorial
from spinon_tutorial_exports import read_spectrum
from test_tutorial_spinon_exports import fixture
from test_tutorial_spinon_provenance import NEW, OLD, compatibility_fixture
from test_tutorial_xxx_spinon_slice import compact_fixture, six_document
from xxx_six_spinon_tutorial_data import read_compact_slice


def different_sources():
    data = compact_fixture()
    for run in data['six_slice']['runs']:
        run['metadata'].update(NEW)
    data['source_compatibility'] = compatibility_fixture()
    refresh(data)
    return data


def inputs(data):
    four = [tutorial.checked_four_run(run, momentum='1.5', exchange='1') for run in data['four_runs']]
    six = read_compact_slice(data['six_slice'], momentum='1.5')['points']
    return four, six


def refresh(data):
    data['points'] = tutorial.join_points(*inputs(data), compatibility=data.get('source_compatibility'))


class CrossSectorSourceJoin(unittest.TestCase):
    def test_default_is_strict_and_explicit_record_survives_roundtrip(self):
        data = different_sources()
        self.assertEqual(tutorial.validate(json.loads(json.dumps(data))), data)
        self.assertEqual(Decimal(data['points'][0]['s2_plus_s4_plus_s6']), Decimal('.6875'))
        with self.assertRaisesRegex(ValueError, 'different source'):
            tutorial.join_points(*inputs(data))
        self.assertEqual(data['four_runs'][0]['metadata']['Bethe revision'], OLD['Bethe revision'])
        self.assertEqual(data['six_slice']['runs'][0]['metadata']['Bethe revision'], NEW['Bethe revision'])

    def test_collect_embeds_record_without_modifying_either_input(self):
        raw = fixture()
        raw['tables']['spectrum']['metadata']['Command'] = 'bethe-xxx-structure-factor-thermo --sectors 2,4'
        six = six_document()
        six['runs'][0]['metadata'].update(NEW)
        record = compatibility_fixture()
        before_six, before_record = deepcopy(six), deepcopy(record)
        with TemporaryDirectory() as directory:
            path = Path(directory)/'four.json';path.write_text(json.dumps(raw))
            result = tutorial.collect([path], six, momentum='1.5', compatibility=record)
        self.assertEqual(six, before_six)
        self.assertEqual(record, before_record)
        self.assertEqual(result['source_compatibility'], record)
        record['evidence'].clear()
        self.assertTrue(result['source_compatibility']['evidence'])
        self.assertEqual(result['six_slice'], six)
        self.assertEqual(tutorial.validate(result), result)

    def test_source_permission_never_waives_physical_identity(self):
        for key, changed in [('Hamiltonian','XXZ Delta=2'), ('Family','another family'),
                             ('Program','another solver'), ('Channel','raising'),
                             ('Exchange J','2'), ('Density convention','another measure')]:
            for value in (changed, None):
                data = different_sources()
                four, six = inputs(data)
                four[0]['metadata'][key] = value
                with self.subTest(key=key,value=value), self.assertRaisesRegex(ValueError, key):
                    tutorial.join_points(four, six, compatibility=data['source_compatibility'])
        four, six = inputs(different_sources())
        four[0]['metadata']['Density relative tolerance'] = '1e-8'
        four[0]['metadata']['Exchange J'] = '1.000'
        self.assertEqual(len(tutorial.join_points(four,six,compatibility=compatibility_fixture())),1)

    def test_unapproved_wrong_sector_and_dirty_sources_are_rejected(self):
        for revision in (NEW['Bethe revision'], 'f'*40, OLD['Bethe revision']+'-dirty'):
            data = different_sources()
            data['four_runs'][0]['metadata']['Bethe revision'] = revision
            with self.subTest(revision=revision), self.assertRaises(ValueError):
                tutorial.validate(data)
        data = different_sources()
        data['source_compatibility']['sector_sources']['2,4'] = [NEW]
        with self.assertRaises(ValueError): tutorial.validate(data)

    def test_incomplete_six_grid_keeps_failed_point_and_partial_sum_absent(self):
        data = different_sources()
        data['six_slice']['frequencies'].append('3')
        data['six_slice']['unavailable'].append({'omega':'3','status':'failed',
             'reason':'density recovery exhausted','raw_sha256':'d'*64})
        extra = deepcopy(data['four_runs'][0])
        for row in extra['spectrum']: row['omega'] = '3'
        data['four_runs'].append(extra)
        refresh(data)
        self.assertIs(tutorial.validate(data),data)
        point = data['points'][1]
        self.assertEqual(point['s6_status'],'failed')
        self.assertIsNone(point['s6']);self.assertIsNone(point['s2_plus_s4_plus_s6'])
        self.assertEqual(Decimal(point['s2_plus_s4']),Decimal('.625'))
        self.assertEqual(point['reason'],'density recovery exhausted')

    def test_uncomputed_future_source_is_not_implicitly_approved(self):
        data = different_sources()
        data['six_slice']['runs'] = []
        data['six_slice']['unavailable'] = [{'omega':'2.0','status':'unavailable',
                     'reason':'not calculated','raw_sha256':None}]
        with self.assertRaisesRegex(ValueError,'source'):
            refresh(data)

    def test_embedded_evidence_and_plot_values_are_revalidated(self):
        mutations = [lambda d: d['source_compatibility'].update(evidence=[]),
                     lambda d: d['source_compatibility']['evidence'][0].update(sha256='bad'),
                     lambda d: d['source_compatibility']['sources'].append({'Bethe revision':'f'*40,'Uni20 revision':'b'*40}),
                     lambda d: d['points'][0].update(six_run=1),
                     lambda d: d['points'][0].update(s2_plus_s4_plus_s6='1'),
                     lambda d: d['six_slice'].update(frequencies=['2.0','3'])]
        for mutation in mutations:
            data = different_sources();mutation(data)
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):tutorial.validate(data)

    def test_cli_check_imports_record_without_writing_data_or_figure(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);raw=fixture()
            raw['tables']['spectrum']['metadata']['Command']='bethe-xxx-structure-factor-thermo --sectors 2,4'
            six=six_document();six['runs'][0]['metadata'].update(NEW)
            for name,obj in [('four',raw),('six',six),('compatibility',compatibility_fixture())]:
                (root/(name+'.json')).write_text(json.dumps(obj))
            command=['plot','--four-json',str(root/'four.json'),'--six-data',str(root/'six.json'),
                     '--source-compatibility',str(root/'compatibility.json'),'--momentum','1.5',
                     '--data',str(root/'combined.json'),'--figure',str(root/'figure.svg'),'--check']
            with patch('sys.argv',command):tutorial.main()
            self.assertFalse((root/'combined.json').exists());self.assertFalse((root/'figure.svg').exists())
            with patch('sys.argv',['plot','--source-compatibility',str(root/'compatibility.json'),'--check']):
                with self.assertRaises(SystemExit):tutorial.main()


if __name__ == '__main__':
    unittest.main()
