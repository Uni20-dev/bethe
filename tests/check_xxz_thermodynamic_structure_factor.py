"""Massive XXZ spectral CLI: precision, branch metadata, exports and failures."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128, mpfr = sys.argv[1:]
precisions = ['fp64', 'long-double'] + (['fp128'] if fp128.upper() in ('ON', 'TRUE', '1') else []) + (['bits-160'] if mpfr.upper() in ('ON', 'TRUE', '1') else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')


def run(args, code=0):
    p = subprocess.run([program, *map(str, args)], capture_output=True, text=True, env=env, timeout=90)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def tables(args, code=0):
    document = json.loads(run([*args, '--format', 'json'], code))
    assert document['status'] == 'complete'
    for table in document['tables'].values():
        assert table['summary']['Outcome'] == ('success' if code == 0 else 'partial')
    return document['tables']


def rows(table):
    return [dict(zip([column['id'] for column in table['columns']], [cell['value'] if isinstance(cell, dict) else cell for cell in row])) for row in table['rows']]


base = ['--delta', '1.5', '--momentum', '1.5', '--omega', '2']
for precision in precisions:
    result = tables([*base, '--precision', precision])
    assert set(result) == {'spectrum', 'continuum'}
    point = rows(result['spectrum'])[0]
    assert point['status'] == 'converged'
    with localcontext() as context:
        context.prec = 80
        expected = Decimal('4.159749286581976836416246520552805055874094165803454070291')
        assert abs(Decimal(str(point['density']))-expected) < Decimal(str(point['estimated_error']))
        tolerance = {'fp64': '1e-10', 'long-double': '1e-13', 'fp128': '1e-28', 'bits-160': '1e-42'}[precision]
        assert abs(Decimal(str(point['density']))-expected) < Decimal(tolerance)
    bands = rows(result['continuum'])
    assert {r['sheet'] for r in bands} == {'C_plus', 'C_minus'}
    assert all(r['present'] for r in bands)
    failure = rows(tables([*base, '--precision', precision, '--max-terms', '0'], 2)['spectrum'])[0]
    assert failure['density'] is None and failure['status'] == 'numerical_failure'
    grid = ['--delta', '1.5', '--momentum-points', '7', '--points', '21', '--precision', precision]
    serial, parallel = tables(grid), tables([*grid, '--threads', '2', '--no-retain'])
    assert serial['spectrum']['rows'] == parallel['spectrum']['rows']
    assert any(r['status'] == 'threshold' for r in rows(serial['spectrum']))
    k0 = rows(tables(['--delta', '1.5', '--momentum', '0', '--omega', '2', '--precision', precision])['spectrum'])[0]
    assert float(k0['density']) > 0

with tempfile.TemporaryDirectory() as directory:
    p = Path(directory)
    target = p/'result.json'
    run([*base, '--json', target, '--quiet'])
    previous = target.read_text()
    for bad in [['--delta', '1'], ['--delta', 'nan'], ['--delta', 'inf'], ['--exchange', '0'],
                ['--momentum', '-1'], ['--omega', 'inf'], ['--points', '1'], ['--tolerance', '0'],
                ['--threads', '0'], ['--csv-table', f'failed={p/"bad.csv"}'], ['--channel', 'zz'],
                ['--momentum-points', '1000000', '--points', '1000000']]:
        delta = [] if '--delta' in bad else ['--delta', '1.5']
        run([*delta, *bad, '--json', target, '--force'], 1)
        assert target.read_text() == previous
    assert not (p/'bad.csv').exists()
    for mode in ('csv', 'tsv'):
        run([*base, f'--{mode}', p/mode, '--stream', '--no-retain', '--quiet'])
        assert 'density' in (p/mode).read_text()
    screen = run([*base, '--format', 'csv', '--table', 'spectrum', '--csv-table', f'continuum={p/"bands.csv"}'])
    assert 'Two-spinon spectral sheets' not in screen
    assert 'C_plus' in (p/'bands.csv').read_text()
    help_text = run(['--help'])
    assert 'spectrum' in help_text and 'continuum' in help_text and 'Used for:' not in help_text
    assert 'Caux' in run(['--references'])


    if 'bits-160' in precisions:
        high = tables(['--delta', '1.5', '--momentum', '1.500000000000000000000000000000000000000000000000001',
                       '--omega', '2', '--precision', 'digits-100'])
        for t in high.values():
            assert int(t['metadata']['Working precision (bits)']) == 333
            for row in t['rows']:
                for cell in row:
                    if isinstance(cell, dict): assert cell['precision'] == 333
        with localcontext() as context:
            context.prec = 120
            k = Decimal(rows(high['spectrum'])[0]['k'])
            assert abs(k - Decimal('1.500000000000000000000000000000000000000000000000001')) < Decimal('1e-95')
print('Massive XXZ spectral CLI contracts passed')
