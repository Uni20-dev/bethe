"""XXZ arbitrary precision, range failures and complete export contracts."""
import csv
from decimal import Decimal, localcontext
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, enabled = sys.argv[1:]
enabled = enabled.upper() in ('ON', 'TRUE', '1')
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')

def run(args, code=0):
    p = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=120)
    assert p.returncode == code, (args, p.returncode, p.stdout, p.stderr)
    return p

def table(args, code=0):
    d = json.loads(run([*args, '--format', 'json'], code).stdout)
    assert d['status'] == 'complete'
    return d['tables']['dispersion']

def rows(t):
    return [dict(zip([c['id'] for c in t['columns']],
                    [x['value'] if isinstance(x, dict) else x for x in row])) for row in t['rows']]

if not enabled:
    assert 'UNI20_ENABLE_MPREAL=ON' in run(['--delta', '2', '--precision', 'bits-128'], 1).stderr
    raise SystemExit(0)
with localcontext() as ctx:
    ctx.prec = 120
    reference = Decimal('0.19490113571001470727701163515916987988002921827177622392600310333179')
    previous = Decimal(1)
    for bits, limit in [(80, '1e-21'), (160, '1e-44'), (256, '1e-62')]:
        args = ['--delta', '2', '--points', '3', '--precision', f'bits-{bits}']
        t = table(args)
        assert int(t['metadata']['Working precision (bits)']) == bits
        values = rows(t)
        assert len(values) == 6 and all(r['status'] == 'converged' for r in values)
        error = abs(Decimal(values[0]['energy']) - reference)
        assert error < Decimal(limit) and error < previous*Decimal('1e-15')
        previous = error
        assert t['rows'] == table([*args, '--no-retain'])['rows']
        for fmt, separator in [('csv', ','), ('tsv', '\t')]:
            payload = run([*args, '--format', fmt]).stdout
            parsed = list(csv.DictReader(io.StringIO('\n'.join(line for line in payload.splitlines()
                                                               if not line.startswith('#'))), delimiter=separator))
            assert len(parsed) == len(values)
            for actual, expected in zip(parsed, values):
                assert actual.keys() == expected.keys()
                for key, value in expected.items():
                    if value is None:
                        assert actual[key] == ''
                    else:
                        try:
                            assert Decimal(actual[key]) == Decimal(str(value)), (key, actual[key], value)
                        except ArithmeticError:
                            assert actual[key] == str(value), (key, actual[key], value)
    args = ['--delta', '2', '--momentum', '0', '--precision', 'bits-160']
    failed = table([*args, '--max-evaluations', '0'], 2)
    assert failed['metadata']['Bulk status'] == 'evaluation_limit'
    assert 'Bulk energy/site' not in failed['metadata']
    assert all(r['status'] == 'converged' for r in rows(failed))
    underflow = table(['--delta', '1.000000000000000000000000000001', '--points', '2', '--precision', 'bits-160'], 2)
    assert 'Bulk energy/site' in underflow['metadata']
    assert all(r['status'] == 'precision_limit' and r['energy'] is None for r in rows(underflow))
    tiny = table(['--delta', '2', '--momentum', '0', '--precision', 'digits-100', '--tolerance', '1e-80'])
    assert int(tiny['metadata']['Working precision (bits)']) >= 333
    assert abs(Decimal(tiny['metadata']['Bulk tolerance'])/Decimal('1e-80')-1) < Decimal('1e-95')
with tempfile.TemporaryDirectory() as directory:
    output = Path(directory)/'existing.json'
    for args in [['--delta', '-1'], ['--delta', 'nan'], ['--delta', '2', '--tolerance', 'nan'],
                 ['--delta', '2', '--precision', 'bits-0']]:
        output.write_text('sentinel')
        if '--precision' not in args:
            args += ['--precision', 'bits-160']
        run([*args, '--json', str(output), '--force'], 1)
        assert output.read_text() == 'sentinel'
print('XXZ explicit precision and complete export contracts passed')
