"""Whole-observation precision stages preserve physical results and retry work."""
from decimal import Decimal, localcontext
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr_flag, fp128_flag = sys.argv[1:]
mpfr = mpfr_flag.upper() in ('ON', 'TRUE', '1')
fp128 = fp128_flag.upper() in ('ON', 'TRUE', '1')
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')


def invoke(arguments, expected=0):
    result = subprocess.run([program, *map(str, arguments)], env=env,
                            capture_output=True, text=True, timeout=180)
    assert result.returncode == expected, (arguments, result.returncode, result.stdout[-2000:], result.stderr)
    return result


def rows(document, name):
    table = document['tables'][name]
    names = [c['id'] for c in table['columns']]
    return [dict(zip(names, [x['value'] if isinstance(x, dict) else x for x in row]))
            for row in table['rows']]


def close(actual, expected):
    with localcontext() as context:
        context.prec = 90
        a, b = Decimal(actual), Decimal(expected)
        assert abs(a - b) <= Decimal('1e-12') * max(abs(a), abs(b), Decimal('1e-100')), (a, b)


point = ['--sectors', '6', '--momentum', '1.5', '--omega', '2',
         '--replicates', '2', '--samples', '4', '--seed', '42', '--threads', '2', '--samples-per-task', '3']
record_options = ['--diagnostics', '--replicate-output', '--format', 'json']


def query(arguments):
    # CLI11 intentionally rejects duplicate scalar options; replace base inputs
    # rather than letting an unrelated parser error stand in for validation.
    base = [value for index in range(0, len(point), 2) if point[index] not in arguments
            for value in point[index:index + 2]]
    return [*base, *arguments]


def document(arguments, expected=0):
    return json.loads(invoke([*query(arguments), *record_options], expected).stdout)


with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    protected, named, absent = [directory / name for name in ('keep.json', 'keep.csv', 'absent.tsv')]

    def rejected(arguments):
        protected.write_text('keep JSON\n')
        named.write_text('keep CSV\n')
        result = invoke([*query(arguments), '--json', protected, '--csv-table', f'spectrum={named}',
                         '--tsv-table', f'spectrum={absent}', '--force'], 1)
        assert result.stderr
        assert protected.read_text() == 'keep JSON\n'
        assert named.read_text() == 'keep CSV\n'
        assert not absent.exists()

    for arguments in [
            ['--six-precision-stages', 'fp80,reference', '--precision', 'fp64'],
            ['--six-precision-stages', 'reference', '--precision', 'bits-192', '--sectors', '2'],
            ['--six-precision-stages', 'fp128,fp80,reference', '--precision', 'bits-192'],
            ['--six-precision-stages', 'reference,reference', '--precision', 'bits-192'],
            ['--six-precision-stages', 'reference,', '--precision', 'bits-192']]:
        rejected(arguments)
    if not mpfr:
        rejected(['--six-precision-stages', 'reference', '--precision', 'bits-192'])
        print('stage capability and preflight checks passed; MPFR is disabled')
        sys.exit(0)

    native = json.loads(invoke(['--momentum', '1.5', '--omega', '2', '--precision', 'long-double',
                               '--format', 'json']).stdout)
    fp80 = native['tables']['spectrum']['metadata']['Arithmetic format'] == 'extended80'
    policies = ['reference']
    if fp80:
        policies.append('fp80,reference')
    if fp128:
        policies.append('fp128,reference')
    if fp80 and fp128:
        policies.append('fp80,fp128,reference')
    fixed = document(['--precision', 'bits-192'])
    assert 'recovery' not in fixed['tables']
    assert fixed['tables']['spectrum']['metadata']['Geometry and output bits'] == '192'
    reference_replicates = rows(fixed, 'replicates')

    for policy in policies:
        result = document(['--precision', 'bits-192', '--six-precision-stages', policy])
        metadata = result['tables']['spectrum']['metadata']
        assert 'Geometry and output bits' not in metadata
        assert 'Density comparison bits' not in metadata and 'Density guard bits' not in metadata
        assert metadata['Accumulation and output bits'] == '192'
        assert metadata['Reference density comparison bits'] == '128'
        assert metadata['Reference density guard bits'] == '192'
        assert metadata['Six-spinon precision stages'] == policy
        stages = policy.split(',')
        for stage in stages:
            prefix = 'Six-spinon ' + stage + ' '
            bits = {'fp80': '64', 'fp128': '113', 'reference': '192'}[stage]
            assert metadata[prefix + 'geometry bits'] == bits
            assert metadata[prefix + 'density guard bits'] == bits
        for table in result['tables'].values():
            assert table['metadata'] == metadata
            for row in table['rows']:
                for cell in row:
                    if isinstance(cell, dict):
                        assert cell['precision'] == 192
        close(rows(result, 'spectrum')[0]['density'], rows(fixed, 'spectrum')[0]['density'])
        if policy == 'reference':
            for table in ('spectrum', 'diagnostics', 'replicates'):
                assert result['tables'][table]['rows'] == fixed['tables'][table]['rows'], table
        accounting = rows(result, 'recovery')
        for old, physical in zip(reference_replicates, rows(result, 'replicates')):
            assert old['replicate'] == physical['replicate']
            assert old['supported_roots'] == physical['supported_roots']
            close(physical['density'], old['density'])
            selected = [r for r in accounting if r['replicate'] == physical['replicate']]
            totals = {r['stage']: r for r in selected if r['kind'] == 'stage'}
            assert list(totals) == stages
            assert int(totals[stages[0]]['attempts']) == int(physical['completed_samples'])
            accepted_roots = accepted_density = accepted_samples = 0
            for i, stage in enumerate(stages):
                total = totals[stage]
                failures = [r for r in selected if r['kind'] == 'failure' and r['stage'] == stage]
                assert all(total[key] == '' for key in ('failure', 'density_status', 'density_failure', 'reason'))
                assert all(r['accepted'] is None and r['retried'] is None and r['failure'] for r in failures)
                assert int(total['accepted']) + sum(int(r['attempts']) for r in failures) == int(total['attempts'])
                assert int(total['retried']) == (int(totals[stages[i + 1]]['attempts']) if i + 1 < len(stages) else 0)
                accepted_samples += int(total['accepted'])
                accepted_roots += int(total['supported_roots']) - sum(int(r['supported_roots']) for r in failures)
                accepted_density += int(total['density_attempts']) - sum(int(r['density_attempts']) for r in failures)
            assert accepted_samples == int(physical['completed_samples'])
            assert accepted_roots == int(physical['supported_roots'])
            assert accepted_density == int(physical['density_attempts'])

    policy = policies[-1]
    mixed = ['--precision', 'bits-192', '--six-precision-stages', policy]
    boundary = document([*mixed, '--momentum', '0'])
    assert rows(boundary, 'spectrum')[0]['status'] == 'physical_zero'
    assert rows(boundary, 'spectrum')[0]['density'] == '0'
    assert rows(boundary, 'replicates') == rows(boundary, 'recovery') == []

    failed = document([*mixed, '--density-relative-tolerance', '1e-100'], 2)
    assert rows(failed, 'spectrum')[0]['density'] is None
    assert failed['tables']['spectrum']['summary']['Outcome'] == 'partial'
    assert failed['tables']['spectrum']['summary']['Status'] == 'incomplete'
    failures = [r for r in rows(failed, 'recovery') if r['kind'] == 'failure']
    assert any(r['stage'] == 'reference' and r['density_status'] == 'exhausted' for r in failures)
    assert all(r['reason'] and r['accepted'] is None and r['retried'] is None for r in failures)
    assert any(r['failed_sample'] is not None for r in rows(failed, 'replicates'))

    timed = document([*mixed, '--density-seconds', '1e-9'], 2)
    counters = rows(timed, 'recovery')
    assert any(r['kind'] == 'failure' and r['density_status'] == 'time_limit' for r in counters)
    assert all(int(r['attempts']) == 0 for r in counters if r['kind'] == 'stage' and r['stage'] != policy.split(',')[0])

    rejected([*mixed, '--replicate-output', '--detail-row-limit', '1'])
    if fp80:
        rejected(['--precision', 'bits-192', '--six-precision-stages', 'fp80,reference', '--density-guard-bits', '63'])
        rejected(['--precision', 'bits-192', '--six-precision-stages', 'fp80,reference', '--exchange', '1e4932', '--omega', '6e4932'])
    if fp128:
        rejected(['--precision', 'bits-192', '--six-precision-stages', 'fp128,reference', '--density-guard-bits', '112'])
    else:
        rejected(['--precision', 'bits-192', '--six-precision-stages', 'fp128,reference'])

    csv_path, tsv_path = directory/'recovery.csv', directory/'recovery.tsv'
    displayed = invoke([*point, *mixed, '--csv-table', f'recovery={csv_path}',
                        '--tsv-table', f'recovery={tsv_path}', '--no-preamble'])
    assert 'Six-spinon precision recovery' not in displayed.stdout
    with csv_path.open() as stream:
        csv_rows = list(csv.DictReader(stream))
    with tsv_path.open() as stream:
        tsv_rows = list(csv.DictReader(stream, delimiter='\t'))
    assert csv_rows == tsv_rows and len(csv_rows) >= 2 * len(policy.split(','))
    assert {r['stage'] for r in csv_rows} == set(policy.split(','))

print('precision stage exports, physical comparisons, failure accounting and file protection passed')
