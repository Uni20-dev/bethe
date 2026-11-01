"""Reviewed source cohorts cannot hide changed controls, gaps or replacements."""
from copy import deepcopy
from decimal import Decimal
from importlib.util import find_spec
import json
from pathlib import Path
import unittest

from spinon_tutorial_provenance import (PHYSICAL_KEYS, SCHEMA, compatibility_sources,
                                       physical_identity, validate_physics, validate_sources)
from test_tutorial_xxx_six_spinon import fixture


OLD = {'Bethe revision': 'a' * 40, 'Uni20 revision': 'b' * 40}
NEW = {'Bethe revision': 'c' * 40, 'Uni20 revision': 'b' * 40}


def compatibility_fixture():
    return {'schema': SCHEMA, 'sources': [dict(OLD), dict(NEW)],
            'reason': 'Reviewed unchanged two/four algorithms and checked six arithmetic; synthetic evidence',
            'evidence': [{'description': 'Matched observations and source audit', 'sha256': 'd' * 64}],
            'sector_sources': {'2,4': [dict(OLD)], '6': [dict(NEW)]}}


class PhysicalIdentity(unittest.TestCase):
    def test_every_identity_field_is_required_and_missing_is_not_an_agreement(self):
        original = fixture()['tables']['spectrum']['metadata']
        for key in PHYSICAL_KEYS:
            for value in (None, '', '   ', 1):
                with self.subTest(key=key, value=value):
                    metadata = dict(original)
                    metadata[key] = value
                    with self.assertRaisesRegex(ValueError, key):
                        validate_physics([metadata, metadata])
        with self.assertRaisesRegex(ValueError, 'No physical'):
            validate_physics([])

    def test_each_physical_field_is_compared_without_source_compatibility(self):
        original = fixture()['tables']['spectrum']['metadata']
        for key in PHYSICAL_KEYS:
            with self.subTest(key=key):
                changed = dict(original)
                changed[key] = '2' if key == 'Exchange J' else 'different ' + changed[key]
                with self.assertRaisesRegex(ValueError, key):
                    validate_physics([original, changed])
        changed = dict(original, **NEW)
        changed.update({'Exchange J': '1.000', 'Density guard bits': '320',
                        'Selected spinon sectors': '2,4'})
        self.assertEqual(validate_physics([original, changed]), physical_identity(original))
        for value in ('NaN', 'sNaN', 'Infinity', '-Infinity', '0', '-1', 'unknown'):
            with self.subTest(exchange=value), self.assertRaisesRegex(ValueError, 'Exchange J'):
                physical_identity(dict(original, **{'Exchange J': value}))


class SourceCompatibility(unittest.TestCase):
    def test_default_join_stays_strict_and_allowlist_is_sector_scoped(self):
        record = compatibility_fixture()
        metadata = [{**OLD, 'Selected spinon sectors': '2,4'}, {**NEW, 'Selected spinon sectors': '6'}]
        with self.assertRaisesRegex(ValueError, 'different source'):
            validate_sources(metadata)
        self.assertEqual(len(validate_sources(metadata, record)), 2)
        metadata[1]['Selected spinon sectors'] = '2,4'
        with self.assertRaisesRegex(ValueError, 'these sectors'):
            validate_sources(metadata, record)
        metadata[1]['Bethe revision'] = 'f' * 40
        with self.assertRaisesRegex(ValueError, 'differ'):
            validate_sources(metadata, record)

    def test_missing_malformed_or_duplicate_evidence_fails(self):
        changes = [('reason', ''), ('evidence', []), ('sources', [OLD, OLD]),
                   ('evidence', [{'description': 'claimed', 'sha256': 'not-a-hash'}]),
                   ('evidence', [{'description': '', 'sha256': 'd' * 64}]),
                   ('sector_sources', {'6': [OLD, OLD]}), ('sector_sources', {}),
                   ('sector_sources', {'2,4': [OLD], '6': [OLD]}), ('cohort_controls', {})]
        for key, value in changes:
            with self.subTest(key=key, value=value):
                record = compatibility_fixture()
                record[key] = value
                with self.assertRaises(ValueError):
                    compatibility_sources(record)
        record = compatibility_fixture()
        record['evidence'] *= 2
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            compatibility_sources(record)


if __name__ == '__main__':
    unittest.main()
