#!/usr/bin/env python3
"""Reconstruct the small tutorial's native/MPFR comparisons from compact records.

This checks arithmetic agreement and sample identities, not sampling convergence
or omitted rapidity tails. Raw hashes identify the separately retained exports;
the compact reader does not claim to have inspected those archives.
"""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import re

from xxx_six_spinon_tutorial_data import read_compact_slice

DATA = Path(__file__).resolve().parents[1] / 'docs/validation/six-spinon-native-comparison-20261031.json'
Q = '1.57079632679489661923132169163975'


def validate(data):
    if (data['schema'] != 'bethe-six-spinon-native-comparison-v1'
            or [case['seed'] for case in data['cases']] != [42, 2, 35]
            or not re.fullmatch('[0-9a-f]{64}', data['binary_sha256'])):
        raise ValueError('Missing case or native-binary provenance')
    results = {}
    for case in data['cases']:
        checked = []
        for name, low, high in [('native', 64, 113), ('reference', 128, 192)]:
            record = case[name]
            if (record['frequencies'] != ['1.8', '2', '2.2'] or record['rapidity_cutoff'] != '8'
                    or len(record['runs']) != 1 or record['unavailable']):
                raise ValueError('Changed grid or missing observations')
            meta = record['runs'][0]['metadata']
            expected = {'Geometry and output bits': str(high), 'Density comparison bits': str(low),
                        'Density guard bits': str(high), 'Density relative tolerance': '1e-5',
                        'Density absolute tolerance': '0', 'First density residue budget': '64',
                        'Second density residue budget': '80', 'First density tail order': '16',
                        'Second density tail order': '20', 'Six-spinon seed': str(case['seed']),
                        'Six-spinon replicates': '4', 'Samples per replicate': '32',
                        'Scheduler concurrency limit': '2', 'Uni20 revision': data['uni20']}
            if name == 'native':
                expected['Bethe revision'] = data['source']
            if any(meta[key] != value for key, value in expected.items()):
                raise ValueError('Changed arithmetic, sampling controls or source')
            checked.append(read_compact_slice(record, momentum=Q)['points'])
        native_meta = case['native']['runs'][0]['metadata']
        reference_meta = case['reference']['runs'][0]['metadata']
        if any(native_meta[key] != reference_meta[key] for key in native_meta
               if key.startswith(('Six-spinon', 'Recovery seconds'))):
            raise ValueError('Changed proposal or numerical budget')
        differences = {'mean': [], 'sampling_se': [], 'replicate_mean_or_maximum': []}
        with localcontext() as context:
            context.prec = 200

            def change(x, y):
                x, y = Decimal(x), Decimal(y)
                return abs(x-y)/max(abs(x), abs(y)) if x or y else Decimal(0)

            for native, reference in zip(*checked):
                for field in ('omega', 'supported_roots', 'density_attempts', 'status'):
                    if native[field] != reference[field]:
                        raise ValueError('Changed frequency or work counts')
                differences['mean'].append(change(native['s6'], reference['s6']))
                differences['sampling_se'].append(change(native['sampling_se'], reference['sampling_se']))
                a, b = native['replicate_records'], reference['replicate_records']
                if len(a) != len(b) or len(a) != 4:
                    raise ValueError('Missing replicate')
                for x, y in zip(a, b):
                    for field in ('replicate', 'completed_samples', 'supported_roots', 'density_attempts', 'failed_sample'):
                        if x[field] != y[field]:
                            raise ValueError('Changed sample identities or work counts')
                    for field in ('density', 'maximum_absolute_density'):
                        differences['replicate_mean_or_maximum'].append(change(x[field], y[field]))
            maxima = {key: str(max(values)) for key, values in differences.items()}
        if (maxima != case['maximum_relative_changes']
                or max(Decimal(value) for value in maxima.values()) > Decimal('5e-30')):
            raise ValueError('Arithmetic comparison differs from the retained result')
        if case['all_work_counts_and_stream_identities_match'] is not True:
            raise ValueError('Missing identity-check result')
        results[str(case['seed'])] = maxima
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DATA)
    args = parser.parse_args()
    print(json.dumps(validate(json.loads(args.data.read_text())), indent=2))
