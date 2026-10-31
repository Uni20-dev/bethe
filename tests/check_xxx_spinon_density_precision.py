"""Independent MPFR density precision preserves geometry, results and exports."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

spectrum, moments, mpfr, fp128 = sys.argv[1:]
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')


def invoke(program, args, expected=0, error=None):
    result = subprocess.run([program, *map(str, args)], env=env, capture_output=True, text=True, timeout=180)
    assert result.returncode == expected, (args, result.returncode, result.stdout[-3000:], result.stderr)
    if error is not None:
        assert error in result.stderr, (args, result.stderr)
    return result.stdout


def rows(document, name):
    table = document['tables'][name]
    names = [column['id'] for column in table['columns']]
    return [dict(zip(names, [cell['value'] if isinstance(cell, dict) else cell for cell in row]))
            for row in table['rows']]


def close(a, b):
    with localcontext() as context:
        context.prec = 160
        a, b = Decimal(a), Decimal(b)
        assert abs(a - b) <= Decimal('1e-20') * max(abs(a), abs(b), Decimal('1e-100')), (a, b)


def compare(program, arguments, table, value):
    reference = json.loads(invoke(program, [*arguments, '--format', 'json']))
    reduced = json.loads(invoke(program, [*arguments, '--density-guard-bits', '192', '--format', 'json']))
    meta = reduced['tables'][table]['metadata']
    assert meta['Precision'] == 'bits-384'
    if table == 'spectrum':
        assert meta['Geometry and output bits'] == '384'
    assert meta['Density comparison bits'] == '128'
    assert meta['Density guard bits'] == '192'
    default = reference['tables'][table]['metadata']
    assert default['Density comparison bits'] == '320' and default['Density guard bits'] == '384'
    for document in (reference, reduced):
        assert document['tables'][table]['summary']['Outcome'] == 'success'
        for row in document['tables'][table]['rows']:
            for cell in row:
                if isinstance(cell, dict):
                    assert cell['precision'] == 384, cell
    a, b = rows(reference, table), rows(reduced, table)
    assert len(a) == len(b) and a
    for old, new in zip(a, b):
        assert old['sector'] == new['sector'] and old['status'] == new['status']
        close(old[value], new[value])
    if 'replicates' in reference['tables']:
        a, b = rows(reference, 'replicates'), rows(reduced, 'replicates')
        assert len(a) == len(b) and a
        for old, new in zip(a, b):
            for key in ('replicate', 'completed_samples', 'supported_roots', 'failed_sample'):
                assert old[key] == new[key], key
            close(old['density'], new['density'])
            close(old['maximum_absolute_density'], new['maximum_absolute_density'])
    return reference, reduced


for program, point, table in ((spectrum, ['--momentum', '1.5', '--omega', '1.8'], 'spectrum'),
                               (moments, ['--momentum', '1.5'], 'moments')):
    assert '--density-guard-bits' in invoke(program, ['--help'])
    with tempfile.TemporaryDirectory() as directory:
        output, named, absent = [Path(directory) / name for name in ('keep.json', 'keep.csv', 'absent.tsv')]
        cases = [(['--sectors', '4', '--precision', 'fp64', '--density-guard-bits', '128'],
                  'requires --precision bits-N'),
                 (['--sectors', '2', '--density-guard-bits', '128'], 'requires sector 4 or 6')]
        if mpfr.upper() in ('ON', 'TRUE', '1'):
            cases += [(['--sectors', '4', '--precision', 'bits-192', '--density-guard-bits', bits],
                       'must be in 25..geometry bits') for bits in ('0', '24', '193', '2147483648')]
        for controls, error in cases:
            output.write_text('preserve JSON\n')
            named.write_text('preserve CSV\n')
            invoke(program, [*point, *controls, '--json', output, '--csv-table', f'{table}={named}',
                             '--tsv-table', f'{table}={absent}', '--force'], 1, error)
            assert output.read_text() == 'preserve JSON\n'
            assert named.read_text() == 'preserve CSV\n'
            assert not absent.exists()

if fp128.upper() in ('ON', 'TRUE', '1'):
    native = json.loads(invoke(spectrum, ['--momentum', '1.5', '--omega', '1.8',
                                         '--precision', 'long-double', '--format', 'json']))
    native_bits = int(native['tables']['spectrum']['metadata']['Working precision (bits)'])
    if 53 < native_bits < 113:
        # The old 53/113 pair rejects this complete quadrature. The intermediate
        # native pair must pass without relaxing either contour or arithmetic checks.
        document = json.loads(invoke(spectrum, ['--sectors', '4', '--momentum', '1.5', '--omega', '1.8',
            '--precision', 'fp128', '--four-order', '8', '--four-cutoff', '8', '--threads', '2',
            '--density-relative-tolerance', '1e-8', '--density-absolute-tolerance', '1e-30',
            '--density-terms', '64', '--density-tail-order', '16', '--format', 'json']))
        meta = document['tables']['spectrum']['metadata']
        assert meta['Density comparison bits'] == str(native_bits)
        assert meta['Density guard bits'] == '113'
        point, = rows(document, 'spectrum')
        assert point['status'] == 'estimated' and point['evaluations'] == '2304'
        # Independently retained MPFR192 quadrature, identical physical controls.
        close(point['density'], '0.18753256731755974112333306196814292064896682968245603409552')

if mpfr.upper() in ('ON', 'TRUE', '1'):
    common = ['--precision', 'bits-384', '--threads', '2', '--density-relative-tolerance', '1e-5',
              '--density-absolute-tolerance', '1e-30', '--diagnostics']
    four = ['--sectors', '4', '--four-order', '2', '--four-cutoff', '2']
    compare(spectrum, [*common, *four, '--momentum', '1.5', '--omega', '1.8'], 'spectrum', 'density')
    six = ['--sectors', '6', '--momentum', '1.5', '--omega', '2', '--replicates', '2', '--samples', '4',
           '--seed', '42', '--six-cutoff', '8', '--replicate-output']
    compare(spectrum, [*common, *six], 'spectrum', 'density')
    compare(moments, [*common, *four, '--momentum', '1.5'], 'moments', 'value')
    compare(moments, [*common, '--sectors', '4', '--global', '--four-order', '2',
                       '--four-gap-cutoffs', '1,1,1'], 'moments', 'value')
print('Density precision overrides preserve geometry, replicate aggregates, values and protected exports')
