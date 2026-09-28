"""Thermodynamic XXX density: precision, statuses, scaling and table contracts."""
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, fp128 = sys.argv[1:]
precisions = ['fp64', 'long-double'] + (['fp128'] if fp128.upper() in ('ON', 'TRUE', '1') else [])
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def run(args, code=0):
    p = subprocess.run([program, *map(str, args)], capture_output=True, text=True, env=env, timeout=90)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout

def tables(args, code=0):
    d = json.loads(run([*args, '--format', 'json'], code))
    assert d['status'] == 'complete'
    for t in d['tables'].values():
        assert t['summary']['Outcome'] == ('success' if code == 0 else 'partial')
        assert Decimal(t['summary']['Run CPU seconds']) >= 0
    return d['tables']

def rows(table):
    return [dict(zip([c['id'] for c in table['columns']], row)) for row in table['rows']]

for precision in precisions:
    base = ['--momentum', '1.5', '--omega', '1.8', '--precision', precision]
    z = tables(base)
    assert set(z) == {'spectrum', 'continuum'}
    point = rows(z['spectrum'])[0]
    assert point['status'] == 'converged'
    assert abs(float(point['density']) - .87313045511634) < 1e-12
    raised = rows(tables([*base, '--channel', 'raising'])['spectrum'])[0]
    with localcontext() as ctx:
        ctx.prec = 80
        assert abs(Decimal(raised['density']) - 2*Decimal(point['density'])) < Decimal('1e-14')
    parallel = tables([*base, '--threads', '4'])
    assert parallel['spectrum']['rows'] == z['spectrum']['rows']
    failure = rows(tables([*base, '--max-evaluations', '0'], 2)['spectrum'])[0]
    assert failure['density'] is None and failure['status'] == 'numerical_failure'
    threshold = tables(['--momentum-points', '3', '--points', '3', '--precision', precision])
    singular = [r for r in rows(threshold['spectrum']) if r['status'] == 'lower_threshold']
    assert len(singular) == 1 and singular[0]['density'] is None
    grid = ['--momentum-points', '9', '--points', '17', '--precision', precision]
    serial, parallel = tables(grid), tables([*grid, '--threads', '4'])
    assert serial['spectrum']['rows'] == parallel['spectrum']['rows']

with tempfile.TemporaryDirectory() as directory:
    p = Path(directory)
    destination = p/'result.json'
    base = ['--momentum', '1.5', '--omega', '1.8']
    run([*base, '--json', destination, '--quiet'])
    before = destination.read_text()
    for bad in [['--exchange', '0'], ['--exchange', 'nan'], ['--momentum', '-1'], ['--momentum', '7'],
                ['--omega', 'inf'], ['--points', '1'], ['--momentum-points', '1'], ['--tolerance', '0'],
                ['--tolerance', 'nan'], ['--omega-min', '4'], ['--threads', '0'],
                ['--channel', 'xx'], ['--table', 'failed'], ['--csv-table', f'failed={p/"bad.csv"}'],
                ['--momentum-points', '1000000', '--points', '1000000']]:
        run([*bad, '--json', destination, '--force'], 1)
        assert destination.read_text() == before
    assert not (p/'bad.csv').exists()
    for mode in ['csv', 'tsv']:
        run([*base, f'--{mode}', p/mode, '--stream', '--no-retain', '--quiet'])
        assert 'density' in (p/mode).read_text()
    screen = run([*base, '--format', 'csv', '--table', 'spectrum', '--csv-table', f'continuum={p/"edges.csv"}'])
    assert 'Two-spinon continuum' not in screen
    assert (p/'edges.csv').exists()
    help_text = run(['--help'])
    assert 'spectrum' in help_text and 'continuum' in help_text and 'Used for:' not in help_text
    assert 'Caux' in run(['--references'])
print('Thermodynamic XXX CLI contracts passed')
