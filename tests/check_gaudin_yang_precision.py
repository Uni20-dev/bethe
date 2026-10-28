"""Gaudin-Yang MPFR input, endpoint/failure, provider and export contracts."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr, mplapack = sys.argv[1:]
enabled = all(v.upper() in ('ON', 'TRUE', '1') for v in (mpfr, mplapack))
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')
pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647')


def run(args, code=0):
    result = subprocess.run([program, *args], env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == code, (args, result.stdout, result.stderr)
    return result


def tables(args, code=0):
    return json.loads(run([*args, '--format', 'json'], code).stdout)['tables']


def rows(table):
    return [dict(zip([c['id'] for c in table['columns']],
                     [v['value'] if isinstance(v, dict) else v for v in row])) for row in table['rows']]


def check_precision(document, bits):
    for table in document.values():
        assert table['metadata']['Arithmetic format'] == 'arbitrary binary'
        assert int(table['metadata']['Working precision (bits)']) == bits
        for row in table['rows']:
            for value in row:
                if isinstance(value, dict):
                    assert value['precision'] == bits, value


if enabled:
    for precision, bits in (('bits-128', 128), ('bits-256', 256), ('digits-100', 333)):
        args = ['2', '--length', '1', '--c', str(pi), '--precision', precision, '--roots']
        result = tables(args)
        check_precision(result, bits)
        row = rows(result['states'])[0]
        assert row['converged']
        with localcontext() as ctx:
            ctx.prec = 140
            assert abs(Decimal(row['energy']) - pi*pi/2) < Decimal(2)**(14-bits)
        for n, sz, c in ((0, '0', '1'), (2, '0', '0'), (4, '2', '1'), (4, '-1', '1')):
            state = tables([str(n), '--sz', sz, '--length', '4', '--c', c, '--precision', precision, '--roots'])
            check_precision(state, bits)
            assert rows(state['states'])[0]['converged']
            if n == 0: assert Decimal(rows(state['states'])[0]['energy']) == 0
        failed = tables(['6', '--length', '6', '--c', '0.01', '--precision', precision,
                         '--max-iterations', '0', '--roots'], 2)
        check_precision(failed, bits)
        row = rows(failed['states'])[0]
        assert not row['converged'] and Decimal(row['root_c']) > Decimal('0.01')

    length = '1.00000000000000000000000000000000000000000000000000000000000000000000001'
    c = '1.000000000000000000000000000000000000000000000000000000000001'
    tolerance = '1.000000000000000000000000000000000000000000001e-90'
    a = tables(['2', '--length', length, '--c', c, '--tolerance', tolerance, '--precision', 'digits-100'])
    b = tables(['2', '--precision', 'bits-333', '--c', c, '--length', length, '--tolerance', tolerance])
    assert a['states']['rows'] == b['states']['rows']
    with localcontext() as ctx:
        ctx.prec = 140
        metadata = a['states']['metadata']
        for key, value in (('Length', length), ('c', c), ('Residual tolerance', tolerance)):
            assert abs(Decimal(metadata[key])/Decimal(value)-1) < Decimal('1e-98')

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        run(['6', '--length', '6', '--c', '1', '--precision', 'digits-100', '--no-preamble',
             '--csv', str(path/'states.csv'), '--json', str(path/'states.json'),
             '--tsv-table', f'charge_roots={path/"roots.tsv"}'])
        exported = json.loads((path/'states.json').read_text())
        if 'tables' in exported: exported = exported['tables']['states']
        expected = [{k: '' if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                     for k, v in row.items()} for row in rows(exported)]
        assert list(csv.DictReader((path/'states.csv').read_text().splitlines())) == expected
        root_rows = list(csv.DictReader((path/'roots.tsv').read_text().splitlines(), delimiter='\t'))
        document = tables(['6', '--length', '6', '--c', '1', '--precision', 'digits-100', '--roots'])
        assert [Decimal(r['k']) for r in root_rows] == [Decimal(r['k']) for r in rows(document['charge_roots'])]

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/'preserve.csv'
    path.write_text('preserve\n')
    invalid = [['--precision', p] for p in ('bits-0', 'digits--2')]
    if not enabled: invalid += [['--precision', p] for p in ('bits-128', 'digits-100')]
    else: invalid += [['--precision', 'digits-100', '--tolerance', t] for t in ('0', 'nan')]
    for opts in invalid:
        result = run(['2', '--length', '1', '--c', '1', *opts, '--csv', str(path), '--force'], 1)
        assert not result.stdout and path.read_text() == 'preserve\n'
    # Higher arithmetic precision does not add even-population interacting shells.
    opts = ['--precision', 'bits-128'] if enabled else []
    result = run(['4', '--length', '4', '--c', '1', *opts, '--csv', str(path), '--force'], 1)
    assert not result.stdout and path.read_text() == 'preserve\n'
print('Gaudin-Yang precision contracts passed')
