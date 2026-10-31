#!/usr/bin/env python3
"""Check the retained binary128 grid against matched MPFR128 S4 references.

This verifies the compact record, not the separately archived raw files. The
same finite quadrature is compared; no sampling, quadrature or tail bound follows.
"""
import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import re

from check_xxx_four_spinon_neighbourhood import DATA as ORIGINAL, validate as validate_original
from check_xxx_four_spinon_precision import number, close

DATA = Path(__file__).resolve().parents[1] / 'docs/validation/four-spinon-native-neighbourhood-20261031.json'


def validate(data, original):
    validate_original(original)
    if data['schema'] != 'bethe.four-spinon-native-neighbourhood.v1':
        raise ValueError('Unknown native neighbourhood schema')
    initial = [c for c in original['cases'] if c['precision'] == 'fp64'
               and c['relative_tolerance'] == '1e-8']
    points = [(c['q'], c['omega']) for c in initial if c['status'] != 'physical_zero']
    excluded = [[c['q'], c['omega']] for c in initial if c['status'] == 'physical_zero']
    if data['excluded_physical_zero_points'] != excluded:
        raise ValueError('The original support audit was changed')
    refs = data['references']
    if [(r['q'], r['omega']) for r in refs] != points:
        raise ValueError('Missing or duplicate matched reference')
    native_meta, ref_meta = data['native_metadata'], data['reference_metadata']
    for meta, pair, sectors in [(native_meta, ('64', '113'), '4'), (ref_meta, ('64', '128'), '2,4')]:
        expected = {'Density comparison bits':pair[0], 'Density guard bits':pair[1],
                    'Geometry and output bits':pair[1], 'Selected spinon sectors':sectors,
                    'Working precision (bits)':pair[1],
                    'Channel':'zz', 'Exchange J':'1', 'Four-spinon order':'4',
                    'Four-spinon log-radius cutoff':'8', 'Density absolute tolerance':'1e-30',
                    'First density residue budget':'64', 'Second density residue budget':'80',
                    'First density tail order':'16', 'Second density tail order':'20'}
        expected['Scheduler concurrency limit'] = '2' if sectors == '4' else '24'
        if any(meta[k] != v for k,v in expected.items()):
            raise ValueError('Changed numerical controls or precision pair')
    if (not re.fullmatch('[0-9a-f]{40}', native_meta['Bethe revision'])
            or native_meta['Precision'] != 'fp128' or ref_meta['Precision'] != 'bits-128'
            or ref_meta['Density relative tolerance'] != '1e-8'):
        raise ValueError('Missing format/source identity')
    if native_meta['Uni20 revision'] != ref_meta['Uni20 revision']:
        raise ValueError('Changed numerical dependency')

    def record(row, bits):
        for k in ('computed_q', 'computed_omega', 'density', 'contour_diagnostic', 'arithmetic_diagnostic'):
            if row[k] is not None and not isinstance(row[k], str):
                raise ValueError('Physical values must remain decimal strings')
        if (not close(row['computed_q'], row['q'], Decimal(2) ** (8-bits))
                or not close(row['computed_omega'], row['omega'], Decimal(2) ** (8-bits))):
            raise ValueError('Wrong computed coordinates')
        if not re.fullmatch('[0-9a-f]{64}', row['raw_sha256']):
            raise ValueError('Missing raw identity')
        if int(row['evaluations']) < 0:
            raise ValueError('Invalid work count')
        if row['status'] == 'estimated':
            if number(row['density']) <= 0 or any(number(row[k]) < 0 for k in ('contour_diagnostic','arithmetic_diagnostic')):
                raise ValueError('Missing successful estimate or diagnostics')
        elif row['status'] == 'numerical_failure':
            if row['density'] is not None or not row['reason']:
                raise ValueError('Failed density must remain unavailable')
        else:
            raise ValueError('Unexpected status in the supported grid')

    by_point = {}
    for ref in refs:
        record(ref, 128)
        if ref['status'] != 'estimated':
            raise ValueError('Missing complete reference')
        by_point[ref['q'], ref['omega']] = ref
    cases = data['cases']
    expected = [(tol,q,w) for tol in ('1e-8','1e-6') for q,w in points]
    if [(c['relative_tolerance'],c['q'],c['omega']) for c in cases] != expected:
        raise ValueError('Missing, duplicate or reordered native case')
    summaries = []
    with localcontext() as context:
        context.prec = 160
        for tol in ('1e-8','1e-6'):
            differences, failures, cpu = [], [], Decimal(0)
            for row in (c for c in cases if c['relative_tolerance'] == tol):
                record(row, 113)
                for field in ('application_cpu_seconds','process_cpu_seconds','application_wall_seconds','external_wall_seconds'):
                    if number(row[field]) <= 0:
                        raise ValueError('Missing measured cost')
                cpu += number(row['application_cpu_seconds'])
                ref = by_point[row['q'],row['omega']]
                if row['reference_sha256'] != ref['raw_sha256']:
                    raise ValueError('Wrong reference identity')
                if row['status'] == 'estimated':
                    if row['exit_code'] != 0:
                        raise ValueError('Successful record has a failure exit code')
                    difference = abs(number(row['density'])-number(ref['density'])) / abs(number(ref['density']))
                    if difference > Decimal('1.3e-30'):
                        raise ValueError('Native/wide value disagreement')
                    differences.append(difference)
                else:
                    if row['exit_code'] != 2:
                        raise ValueError('Failed record has a successful exit code')
                    failures.append([row['q'],row['omega']])
            summaries.append({'relative_tolerance':tol,'complete':len(differences),'failed_points':failures,
                              'maximum_relative_change_vs_mpfr128':str(max(differences)),
                              'total_application_cpu_seconds':str(cpu)})
    if summaries != data['summaries']:
        raise ValueError('Summary does not reconstruct from complete cases')
    return summaries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DATA)
    args = parser.parse_args()
    print(json.dumps(validate(json.loads(args.data.read_text()), json.loads(ORIGINAL.read_text())), indent=2))
