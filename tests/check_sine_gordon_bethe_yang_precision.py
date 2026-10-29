"""Explicit Bethe-Yang precision, full ordered output and independent energies."""
import csv
from decimal import Decimal, localcontext, InvalidOperation
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, available = sys.argv[1:3]
long_running = sys.argv[3:] == ["--long-running"]
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def run(args, code=0):
    result = subprocess.run([program, *args], env=env, text=True, capture_output=True, timeout=1200 if long_running else 240)
    assert result.returncode == code, (args, result.returncode, result.stdout, result.stderr)
    return result

def table(args, code=0):
    document = json.loads(run([*args, '--format', 'json'], code).stdout)
    assert document['status'] == 'complete'
    return document['tables']['levels']

def rows(t):
    return [dict(zip([c['id'] for c in t['columns']],
                     [v['value'] if isinstance(v, dict) else v for v in row])) for row in t['rows']]

base = ['--p', '2', '--length', '10']
if available.upper() not in ('ON', 'TRUE', '1'):
    assert 'UNI20_ENABLE_MPREAL=ON' in run([*base, '--precision', 'bits-128'], 1).stderr
    raise SystemExit(0)

if long_running:
    with localcontext() as context:
        context.prec = 130
        t = table([*base, '--precision', 'bits-333', '--tolerance', '1e-80', '--max-iterations', '512'])
        assert rows(t)[0]['status'] == 'converged'
        assert abs(Decimal(t['metadata']['Counting tolerance']) / Decimal('1e-80') - 1) < Decimal('1e-95')
    print('Bethe-Yang 333-bit interacting refinement passed')
    raise SystemExit(0)

with localcontext() as context:
    context.prec = 130
    reference = Decimal('2.0857356187851404880309779322375079570900382654030915686633279793820126035814552563203616916129861034106042386768200734652')
    errors = []
    for bits in (128, 256):
        t = table([*base, '--precision', f'bits-{bits}'])
        row = rows(t)[0]
        assert row['status'] == 'converged'
        assert int(t['metadata']['Working precision (bits)']) == bits
        error = abs(Decimal(row['energy']) - reference)
        errors.append(error)
        assert error < Decimal('1e-32' if bits == 128 else '1e-70'), row
    assert errors[1] < errors[0] * Decimal('1e-30')

    # Scheduling is checked at a fixed numerical budget, separate from bit growth.
    scan = [*base, '--precision', 'bits-160', '--levels', '3', '--tolerance', '1e-12']
    expected = table(scan)
    for extra in (['--threads', '2'], ['--threads', '2', '--no-retain']):
        got = table([*scan, *extra])
        assert got['columns'] == expected['columns'] and got['rows'] == expected['rows']
    for fmt, delimiter in [('csv', ','), ('tsv', '\t')]:
        content = run([*scan, '--format', fmt]).stdout
        actual = list(csv.DictReader(io.StringIO('\n'.join(line for line in content.splitlines()
                                                        if not line.startswith('#'))), delimiter=delimiter))
        assert len(actual) == len(expected['rows'])
        for a, b in zip(actual, rows(expected)):
            assert set(a) == set(b)
            for key, value in b.items():
                if value is None:
                    assert a[key] == ''
                else:
                    try:
                        assert Decimal(a[key]) == Decimal(str(value)), (key, a[key], value)
                    except InvalidOperation:
                        assert a[key] == str(value), (key, a[key], value)

    for flag, status in [('--max-iterations', 'iteration_limit'),
                         ('--max-phase-evaluations', 'phase_limit'),
                         ('--max-fourier-cutoffs', 'phase_limit')]:
        a = table([*scan, flag, '0'], 2)
        b = table([*scan, flag, '0', '--threads', '2'], 2)
        assert a['rows'] == b['rows']
        assert all(r['energy'] is None and r['rapidity'] is None and r['status'] == status for r in rows(a))
    # Parsing is independent of expensive interacting refinement (tested separately).
    t = table(['--p', '1', '--length', '10', '--precision', 'bits-333', '--tolerance', '1e-80'])
    assert rows(t)[0]['status'] == 'converged'
    assert abs(Decimal(t['metadata']['Counting tolerance']) / Decimal('1e-80') - 1) < Decimal('1e-95')
    free = table(['--p', '1', '--length', '10', '--precision', 'bits-160', '--levels', '3', '--threads', '2'])
    assert free['metadata']['Execution mode'] == 'direct serial'
    assert all(int(r['evaluations']) == 0 for r in rows(free))

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory) / 'preserve.json'
    for extra in (['--precision', 'bits-0'], ['--tolerance', 'nan'], ['--number', '1'],
                  ['--p', '0.5'], ['--max-phase-levels', '25']):
        args = [*base, '--precision', 'bits-128']
        for flag, value in zip(extra[::2], extra[1::2]):
            if flag in args:
                args[args.index(flag) + 1] = value
            else:
                args += [flag, value]
        path.write_text('sentinel')
        run([*args, '--json', str(path), '--force'], 1)
        assert path.read_text() == 'sentinel'
print('Sine-Gordon Bethe-Yang MPFR accuracy and complete exports passed')
