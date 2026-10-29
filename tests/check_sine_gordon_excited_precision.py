"""Exact two-hole precision, independent gaps, Async completion and exports."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, available = sys.argv[1:]
enabled = available.upper() in ('ON', 'TRUE', '1')
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def run(args, code=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=300)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout

def result(args, bits, code=0):
    doc = json.loads(run([*args, '--format', 'json'], code))
    assert doc['status'] == 'complete'
    assert set(doc['tables']) == {'levels', 'source', 'gap'}
    for t in doc['tables'].values():
        assert int(t['metadata']['Working precision (bits)']) == bits
        assert t['summary']['Outcome'] == ('success' if code == 0 else 'partial')
        for row in t['rows']:
            for col, cell in zip(t['columns'], row):
                if col['type'] == 'real' and cell is not None:
                    assert isinstance(cell, dict) and cell['precision'] == bits, (col, cell)
    return doc['tables']

def value(cell):
    return cell['value'] if isinstance(cell, dict) else cell

assert 'bits-N' in run(['--help']) and 'digits-N' in run(['--help'])
base = ['--length', '1', '--p', '1']
if enabled:
    pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647')
    previous = Decimal(1)
    for bits in (80, 160, 256, 333):
        args = [*base, '--precision', f'bits-{bits}', '--max-intervals', '8192', '--max-root-iterations', '8192']
        tables = result(args, bits)
        gap = tables['gap']['rows'][0][0]
        with localcontext() as c:
            c.prec = 120
            expected = 2*(1+pi*pi).sqrt()
            error = abs(Decimal(value(gap))-expected)
            assert error < previous*Decimal('1e-15'), (bits, error, previous)
            previous = error
        assert tables['source']['rows'][0][-2:] == [True, 'converged']
    for bits in (80, 160):
        args = [*base, '--precision', f'bits-{bits}']
        direct = result(args, bits)
        concurrent = result([*args, '--threads', '2'], bits)
        assert all(direct[k]['rows'] == concurrent[k]['rows'] for k in direct)
        failed = result([*args, '--max-root-iterations', '0', '--threads', '2'], bits, 2)
        assert failed['levels']['rows'][0][-2:] == [True, 'converged']
        assert failed['levels']['rows'][1][1:3] == [None, None]
        assert failed['source']['rows'][0][2] is None
        assert failed['gap']['rows'][0][:2] == [None, None]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for sep, flag, suffix in ((',', '--csv-table', 'csv'), ('\t', '--tsv-table', 'tsv')):
                paths = {name: root/f'{name}.{suffix}' for name in direct}
                options = [part for name, path in paths.items() for part in (flag, f'{name}={path}')]
                actual = result([*args, *options], bits)
                for name, path in paths.items():
                    lines = path.read_text().splitlines()
                    assert any(line.startswith('#') for line in lines)
                    rows = list(csv.reader((l for l in lines if not l.startswith('#')), delimiter=sep))
                    t = actual[name]
                    assert rows[0] == [c['id'] for c in t['columns']]
                    def exported(c):
                        v = value(c)
                        return str(v).lower() if isinstance(v, bool) else '' if v is None else str(v)
                    assert rows[1:] == [[exported(c) for c in row] for row in t['rows']]
            path = root/'levels.json'
            run([*args, '--threads', '2', '--json', str(path), '--no-retain', '--quiet'])
            actual = json.loads(path.read_text())['tables']
            assert all(direct[k]['rows'] == actual[k]['rows'] for k in direct)
    parsed = result(['--mass', '1.23456789012345678901234567890123456789', '--length', '1', '--p', '2',
                     '--precision', 'digits-100', '--tolerance', '1e-80', '--max-root-iterations', '0',
                     '--max-cutoffs', '0'], 333, 2)
    with localcontext() as c:
        c.prec = 110
        meta = parsed['levels']['metadata']
        assert abs(Decimal(meta['Soliton mass'])-Decimal('1.23456789012345678901234567890123456789')) < Decimal('1e-95')
        assert abs(Decimal(meta['Absolute Y tolerance per state'])/Decimal('1e-80')-1) < Decimal('1e-95')
else:
    assert not run([*base, '--precision', 'bits-160'], 1)
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder)/'keep.json'
    for args in (['--length', '-1', '--p', '1'], ['--length', '1', '--p', '0.5'],
                 [*base, '--number', '1'], [*base, '--tolerance', 'nan'],
                 [*base, '--max-intervals', '1'], [*base, '--contour-shift', '2']):
        for precision in (['fp64', 'bits-160'] if enabled else ['fp64']):
            path.write_text('keep\n')
            run([*args, '--precision', precision, '--threads', '2', '--json', str(path), '--force'], 1)
            assert path.read_text() == 'keep\n'
print('Sine-Gordon excited MPFR accuracy, complete precision, Async and exports passed')
