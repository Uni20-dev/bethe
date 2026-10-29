"""Explicit Lee-Yang arithmetic, independent solves and complete exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

vacuum, excited, available = sys.argv[1:]
enabled = available.upper() in ('ON', 'TRUE', '1')
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')


def run(program, args, status=0):
    done = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=180)
    assert done.returncode == status, (program, args, done.returncode, done.stdout, done.stderr)
    return done.stdout


def value(cell):
    return cell['value'] if isinstance(cell, dict) else cell


def result(program, args, bits, status=0):
    doc = json.loads(run(program, [*args, '--format', 'json'], status))
    assert doc['status'] == 'complete'
    for table in doc['tables'].values():
        assert int(table['metadata']['Working precision (bits)']) == bits
        assert table['summary']['Outcome'] == ('success' if status == 0 else 'partial')
        for row in table['rows']:
            for cell in row:
                if isinstance(cell, dict):
                    assert cell['precision'] == bits
    return doc['tables']


for program in (vacuum, excited):
    help_text = run(program, ['--help'])
    assert 'bits-N' in help_text and 'digits-N' in help_text
    base = ['--mass', '1', '--length', '30', '--tolerance', '1e-10']
    if enabled:
        for bits in (80, 160):
            args = [*base, '--precision', f'bits-{bits}']
            direct = result(program, args, bits)
            if program == excited:
                concurrent = result(program, [*args, '--threads', '2'], bits)
                assert {n: t['rows'] for n, t in direct.items()} == {n: t['rows'] for n, t in concurrent.items()}
            failed = result(program, [*args, '--max-kernel-products', '0'], bits, 2)
            if program == excited:
                concurrent = result(program, [*args, '--threads', '2', '--max-kernel-products', '0'], bits, 2)
                assert {n: t['rows'] for n, t in failed.items()} == {n: t['rows'] for n, t in concurrent.items()}
            with tempfile.TemporaryDirectory() as folder:
                target = Path(folder)
                for sep, flag, extension in ((',', '--csv-table', 'csv'), ('\t', '--tsv-table', 'tsv')):
                    export = []
                    for name in direct:
                        export += [flag, f'{name}={target / (name + "." + extension)}']
                    tables = result(program, [*args, *export], bits)
                    for name, table in tables.items():
                        lines = (target / f'{name}.{extension}').read_text().splitlines()
                        rows = list(csv.reader((line for line in lines if not line.startswith('#')), delimiter=sep))
                        expected = [[str(value(c)).lower() if isinstance(value(c), bool) else '' if value(c) is None else str(value(c)) for c in row] for row in table['rows']]
                        assert rows[0] == [c['id'] for c in table['columns']]
                        assert rows[1:] == expected, (name, rows, expected)
        parsed = result(program, ['--length', '30', '--precision', 'digits-100', '--tolerance', '1e-80', '--max-kernel-products', '0'], 333, 2)
        for table in parsed.values():
            with localcontext() as context:
                context.prec = 110
                assert abs(Decimal(table['metadata']['Absolute Y tolerance' if program == vacuum else 'Absolute Y tolerance per state']) / Decimal('1e-80') - 1) < Decimal('1e-95')
    else:
        rejected = run(program, [*base, '--precision', 'bits-160'], 1)
        assert not rejected
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder)/'keep.json'
        for args in (['--length', '-1'], ['--length', '30', '--max-intervals', '1'], ['--length', '30', '--tolerance', 'nan']):
            output.write_text('keep\n')
            run(program, [*args, '--json', str(output), '--force'], 1)
            assert output.read_text() == 'keep\n'
print('Lee-Yang explicit precision and export checks passed')
