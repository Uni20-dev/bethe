"""Explicit Bose-Fermi precision, reductions, failures and exported zeros."""
import csv
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

program, mpfr = sys.argv[1:]
enabled = mpfr.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')
pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647')


def run(args, code=0):
    result = subprocess.run([program, *args], env=env, text=True, capture_output=True, timeout=60)
    assert result.returncode == code, (args, result.returncode, result.stdout, result.stderr)
    return result


def tables(args, code=0):
    doc = json.loads(run([*args, '--format', 'json'], code).stdout)
    assert doc['status'] == 'complete'
    for table in doc['tables'].values():
        assert table['summary']['Outcome'] == ('success' if code == 0 else 'partial')
    return doc['tables']


def rows(table):
    return [dict(zip([c['id'] for c in table['columns']],
                    [v['value'] if isinstance(v, dict) else v for v in row])) for row in table['rows']]


def precision(document, bits):
    for table in document.values():
        assert table['metadata']['Arithmetic format'] == 'arbitrary binary'
        assert int(table['metadata']['Working precision (bits)']) == bits
        for row in table['rows']:
            for value in row:
                if isinstance(value, dict): assert value['precision'] == bits, value


def args(b, f, length='1', c='1'):
    return ['--bosons', str(b), '--fermions', str(f), '--length', length, '--c', c]


if enabled:
    for name, bits in (('bits-128', 128), ('bits-256', 256), ('digits-100', 333)):
        flags = ['--precision', name, '--roots']
        for b, f in ((1, 1), (2, 0)):
            result = tables([*args(b, f, c=str(pi)), *flags])
            precision(result, bits)
            assert rows(result['states'])[0]['converged']
            with localcontext() as ctx:
                ctx.prec = 140
                assert abs(Decimal(rows(result['states'])[0]['energy'])-pi*pi/2) < Decimal(2)**(14-bits)
        for b, f, c in ((0, 0, '1'), (3, 2, '0'), (0, 4, '1'), (1, 0, '1')):
            result = tables([*args(b, f, c=c), *flags])
            precision(result, bits)
            assert rows(result['states'])[0]['converged']
            if c == '0':
                for row in rows(result['free_modes']):
                    if row['species'] == 'boson': assert Decimal(row['k']) == 0
        for b, f in ((2, 3), (3, 0)):
            failed = tables([*args(b, f), *flags, '--max-iterations', '0'], 2)
            precision(failed, bits)
            state = rows(failed['states'])[0]
            assert state['status'] == 'iteration_limit' and state['energy'] is None and state['momentum'] is None
            for table, key in (('charge_roots', 'k'), ('auxiliary_roots', 'lambda')):
                assert all(row[key] is None for row in rows(failed[table]))

    length = '1.000000000000000000000000000000000000000000000000000000000000000000001'
    c = '1.00000000000000000000000000000000000000000000000000000001'
    tolerance = '1.000000000000000000000000000000000000000000001e-90'
    base = args(1, 1, length=length, c=c)
    a = tables([*base, '--tolerance', tolerance, '--precision', 'digits-100'])
    b = tables([*base, '--precision', 'bits-333', '--tolerance', tolerance, '--no-retain'])
    assert a['states']['rows'] == b['states']['rows']
    with localcontext() as ctx:
        ctx.prec = 140
        metadata = a['states']['metadata']
        for key, value in (('Circumference', length), ('Contact parameter c', c), ('Residual tolerance', tolerance)):
            assert abs(Decimal(metadata[key])/Decimal(value)-1) < Decimal('1e-98')

    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)
        base = [*args(2, 3, length='5'), '--precision', 'digits-100', '--roots']
        flags = ['--json', str(path/'state.json'), '--no-preamble']
        for name in ('states', 'charge_roots', 'auxiliary_roots'):
            for suffix in ('csv', 'tsv'): flags += [f'--{suffix}-table', f'{name}={path/name}.{suffix}']
        run([*base, *flags])
        document = json.loads((path/'state.json').read_text())['tables']
        precision(document, 333)
        for name in document:
            expected = [{k: '' if v is None else str(v).lower() if isinstance(v, bool) else str(v)
                         for k, v in row.items()} for row in rows(document[name])]
            for suffix, delimiter in (('csv', ','), ('tsv', '\t')):
                actual = list(csv.DictReader((path/f'{name}.{suffix}').read_text().splitlines(), delimiter=delimiter))
                assert actual == expected

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/'preserve.json'
    path.write_text('sentinel')
    invalid = [['--precision', p] for p in ('bits-0', 'digits--1')]
    if not enabled: invalid += [['--precision', p] for p in ('bits-128', 'digits-100')]
    else: invalid += [['--precision', 'bits-128', '--tolerance', t] for t in ('0', 'nan')]
    for options in invalid:
        result = run([*args(2, 3), *options, '--json', str(path), '--force'], 1)
        assert not result.stdout and path.read_text() == 'sentinel'
        if not enabled and options[-1] in ('bits-128', 'digits-100'):
            assert 'arbitrary-precision linear solves require' in result.stderr
    # Increased precision retains the mixed-sector shell restriction.
    flags = ['--precision', 'bits-128'] if enabled else []
    result = run([*args(2, 2), *flags, '--json', str(path), '--force'], 1)
    assert not result.stdout and path.read_text() == 'sentinel'
print('Bose-Fermi MPFR precision and export contracts passed')
