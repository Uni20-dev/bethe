"""TASEP explicit precision, complete exports and failure diagnostics."""
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
base = ['5', '--particles', '2']


def run(args, status=0):
    done = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert done.returncode == status, (args, done.returncode, done.stdout, done.stderr)
    return done


def value(cell):
    return cell['value'] if isinstance(cell, dict) else cell


def result(args, bits, status=0):
    doc = json.loads(run([*args, '--format', 'json'], status).stdout)
    assert doc['status'] == 'complete'
    for table in doc['tables'].values():
        assert int(table['metadata']['Working precision (bits)']) == bits
        assert table['summary']['Outcome'] == ('success' if status == 0 else 'partial')
        for row in table['rows']:
            for column, cell in zip(table['columns'], row):
                if column['type'] == 'real' and cell is not None:
                    assert isinstance(cell, dict) and cell['precision'] == bits, (column, cell)
    return doc['tables']


def row(tables):
    table = tables['relaxation']
    return dict(zip([c['id'] for c in table['columns']], map(value, table['rows'][0])))


help_text = run(['--help']).stdout
assert 'bits-N' in help_text and 'digits-N' in help_text
if enabled:
    with localcontext() as ctx:
        ctx.prec = 110
        reference = Decimal('0.710290069123028672493762532365387684613821535715147285599991234743117333842877636225644555')
        previous = Decimal(1)
        for bits in (80, 160, 256):
            flags = ['--precision', f'bits-{bits}', '--roots']
            direct = result([*base, *flags], bits)
            r = row(direct)
            error = abs(Decimal(r['gap']) - reference)
            assert error < Decimal(8192) * Decimal(2) ** (1-bits)
            assert error < previous * Decimal('1e-15')
            previous = error
            if bits == 256:
                assert error < Decimal('1e-70')
            for rate in ('3.125', '1e-100'):
                scaled = row(result([*base, *flags, '--rate', rate], bits))
                assert abs(Decimal(scaled['gap'])/Decimal(rate) - reference) < Decimal(8192)*Decimal(2)**(1-bits)
            streamed = result([*base, *flags, '--no-retain'], bits)
            assert {k: t['rows'] for k, t in direct.items()} == {k: t['rows'] for k, t in streamed.items()}
            for budget, status in (('--max-iterations', 'iteration_limit'), ('--max-seed-iterations', 'seed_limit')):
                failed = result([*base, *flags, budget, '0'], bits, 2)
                fr = row(failed)
                assert fr['status'] == status and not fr['converged']
                assert all(fr[k] is None for k in ('gap', 'frequency', 'lambda_real', 'lambda_imag'))
                assert all(r[1:] == [None, None] for r in failed['roots']['rows'])
            for n in ('0', '6'):
                empty = result(['6', '--particles', n, *flags], bits)
                er = row(empty)
                assert er['status'] == 'stationary_only' and er['converged'] and not er['has_mode']
                assert er['gap'] is None and empty['roots']['rows'] == []
            with tempfile.TemporaryDirectory() as folder:
                target = Path(folder)
                for sep, extension in ((',', 'csv'), ('\t', 'tsv')):
                    export = []
                    for name in direct:
                        export += [f'--{extension}-table', f'{name}={target / (name + "." + extension)}']
                    tables = result([*base, *flags, *export, '--no-retain'], bits)
                    for name, table in tables.items():
                        lines = (target / f'{name}.{extension}').read_text().splitlines()
                        rows = list(csv.reader((line for line in lines if not line.startswith('#')), delimiter=sep))
                        expected = [[str(value(c)).lower() if isinstance(value(c), bool) else '' if value(c) is None else str(value(c)) for c in r] for r in table['rows']]
                        assert rows[0] == [c['id'] for c in table['columns']]
                        assert rows[1:] == expected
        analytic = result(['6', '--particles', '1', '--precision', 'digits-100', '--rate', '3.125', '--max-iterations', '0', '--max-seed-iterations', '0'], 333)
        assert abs(Decimal(row(analytic)['gap']) - Decimal('1.5625')) < Decimal('1e-98')
else:
    rejected = run([*base, '--precision', 'bits-160'], 1)
    assert not rejected.stdout and 'UNI20_ENABLE_MPREAL=ON' in rejected.stderr

with tempfile.TemporaryDirectory() as folder:
    output = Path(folder)/'keep.json'
    precision = 'bits-160' if enabled else 'bits-80'
    for args in (['--precision', precision, '--rate', 'nan'], ['--precision', precision, '--tolerance', '0'],
                 ['--precision', precision, '--max-sites', '4'], ['--precision', 'bits-0']):
        output.write_text('keep\n')
        run([*base, *args, '--json', str(output), '--force'], 1)
        assert output.read_text() == 'keep\n'
print('TASEP explicit precision and export checks passed')
