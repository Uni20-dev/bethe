"""Potts selection/counting, native precision, failures and export contracts."""
import csv
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


def run(args, status=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, timeout=30, env=env)
    assert p.returncode == status, (args, p.returncode, p.stdout, p.stderr)
    return p.stdout


def records(t):
    return [dict(zip([c['id'] for c in t['columns']], row)) for row in t['rows']]


def table(args, status=0):
    doc = json.loads(run([*args, '--format', 'json'], status))
    assert doc['status'] == 'complete'
    t = doc['tables']['levels']
    assert t['summary']['Outcome'] == ('success' if status == 0 else 'partial')
    assert Decimal(t['summary']['Run CPU seconds']) >= 0
    assert 'PBC' in t['metadata']['Hamiltonian']
    assert 'not a complete spectrum' in t['metadata']['State selection']
    return t


for precision in precisions:
    base = ['3', '--precision', precision]
    t = table(base)
    rows = records(t)
    assert len(rows) == 7
    assert [(r['charge'], int(r['k'])) for r in rows] == [(0, 0)] + [(q, k) for k in range(3) for q in (1, -1)]
    assert all(r['status'] == 'converged' for r in rows)
    assert [r['branch'] for r in rows] == ['ground']+['charged-one-hole']*6
    assert rows[0]['gap'] == rows[0]['x_scaled'] == rows[0]['momentum'] == '0'
    with localcontext() as ctx:
        ctx.prec = 80
        tolerance = {'fp64': Decimal('2e-13'), 'long-double': Decimal('1e-16'), 'fp128': Decimal('2e-31')}[precision]
        e0 = Decimal('-7.6846584384264908247321147839611155377207988380604306515979503596')
        eq = Decimal('-6.9243439920202489313622513170042223878990963349288496338373686912')
        e1 = Decimal('-2.3775516419230208751905199448242467572796214788739660597038311918')
        for r, expected in zip(rows, [e0, eq, eq, e1, e1, e1, e1]):
            assert abs(Decimal(r['energy'])-expected) < tolerance, (precision, r)
            assert abs(Decimal(r['gap'])-(expected-e0)) < 2*tolerance
            assert Decimal(r['residual']) <= Decimal(t['metadata']['Residual tolerance'])
        for i in (1, 3, 5):
            assert rows[i]['energy'] == rows[i+1]['energy']
    assert records(table([*base, '--no-retain'])) == rows
    assert records(table([*base, '--branch', 'ground'])) == rows[:1]
    assert records(table([*base, '--branch', 'charged'])) == rows[1:]
    assert records(table([*base, '--momentum-index', '1'])) == rows[:1]+rows[3:5]
    assert records(table([*base, '--charge', '-1'])) == rows[:1]+rows[2::2]
    failed = table([*base, '--max-iterations', '0'], 2)
    assert 'Vacuum energy' not in failed['metadata']
    assert all(r['energy'] is None and r['gap'] is None and r['x_scaled'] is None and
               r['status'] == 'iteration_limit' and r['iterations'] == '0' for r in records(failed))

partial = records(table(['3', '--max-iterations', '28'], 2))
assert partial[0]['energy'] is not None and partial[0]['gap'] == '0'
assert any(r['energy'] is None for r in partial[1:])
missing_reference = records(table(['2', '--momentum-index', '1', '--max-iterations', '20'], 2))
assert missing_reference[0]['energy'] is None
assert all(r['energy'] is not None and r['gap'] is None and r['x_scaled'] is None and
           r['status'] == 'reference_unavailable' for r in missing_reference[1:])
even = records(table(['4']))
assert len(even) == 9 and sum(int(r['k']) == 2 for r in even) == 2  # no duplicated zone boundary

invalid = [['1'], ['-1'], ['3', '--momentum-index', '3'], ['3', '--momentum-index', '-1'],
           ['3', '--charge', '0'], ['3', '--charge', '2'], ['3', '--branch', 'complete'],
           ['3', '--branch', 'ground', '--momentum-index', '0'], ['3', '--branch', 'ground', '--charge', '-1'],
           ['3', '--tolerance', '0'], ['3', '--tolerance', '-1'], ['3', '--tolerance', 'nan'],
           ['3', '--tolerance', 'inf'], ['3', '--max-iterations', '-1'], ['3', '--max-sites', '2'],
           ['513'], ['18446744073709551615', '--max-sites', '18446744073709551615']]
for args in invalid:
    run(args, 1)

with tempfile.TemporaryDirectory() as directory:
    js, cs, ts = [Path(directory)/('potts.'+suffix) for suffix in ('json', 'csv', 'tsv')]
    args = ['4', '--json', str(js), '--csv', str(cs), '--tsv', str(ts), '--format', 'plain']
    screen = run(args)
    assert 'CPU time' in screen and 'Used for:' not in screen
    t = json.loads(js.read_text())['tables']['levels']
    expected = [{k: '' if v is None else str(v) for k, v in row.items()} for row in records(t)]
    for path, delimiter in ((cs, ','), (ts, '\t')):
        text = path.read_text()
        assert 'Hamiltonian' in text and 'Momentum convention' in text
        actual = list(csv.DictReader((line for line in text.splitlines() if not line.startswith('#')), delimiter=delimiter))
        assert actual == expected
    original = js.read_text()
    run(['4', '--json', str(js)], 1)
    for bad in invalid:
        run([*bad, '--json', str(js), '--force'], 1)
        assert js.read_text() == original
    run(['4', '--json', str(js), '--force', '--quiet', '--no-retain'])
    assert json.loads(js.read_text())['tables']['levels']['rows'] == t['rows']
    run(['4', '--stream', '--no-retain', '--format', 'plain'])
    run(['4', '--format', 'pretty'])
    failure = run(['3', '--max-iterations', '0', '--format', 'csv', '--no-preamble'], 2)
    assert all(r['energy'] == r['gap'] == r['x_scaled'] == '' for r in csv.DictReader(failure.splitlines()))

print('Potts CLI contracts passed')
