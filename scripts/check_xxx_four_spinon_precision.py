#!/usr/bin/env python3
"""Compare arithmetic on one fixed four-spinon tutorial query, retaining failures."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess

from tutorial_common import solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data/xxx-four-spinon-precision.json'
FORMATS = ('fp64', 'long-double', 'fp128', 'bits-128', 'bits-192')
CONTROLS = {'Exchange J': '1', 'Channel': 'zz', 'Selected spinon sectors': '2,4',
            'Four-spinon order': '4', 'Four-spinon log-radius cutoff': '8',
            'Density relative tolerance': '1e-8', 'Density absolute tolerance': '1e-30',
            'First density residue budget': '64', 'Second density residue budget': '80',
            'First density tail order': '16', 'Second density tail order': '20',
            'Spectral points': '1', 'Scheduler concurrency limit': '1',
            'Density convention': 'S(q,w); moments use dw/(2*pi); raising=2*zz'}
ARGUMENTS = ['--sectors', '2,4', '--momentum', '1.5', '--omega', '1.8', '--exchange', '1',
             '--channel', 'zz', '--four-order', '4', '--four-cutoff', '8',
             '--density-relative-tolerance', '1e-8', '--density-absolute-tolerance', '1e-30',
             '--density-terms', '64', '--density-tail-order', '16', '--tolerance', '1e-12',
             '--threads', '1', '--diagnostics', '--format', 'json']


def number(value):
    if value is None or isinstance(value, bool):
        raise ValueError('Missing numeric value')
    result = Decimal(str(value))
    if not result.is_finite():
        raise ValueError('Nonfinite numeric value')
    return result


def close(a, b, tolerance):
    with localcontext() as context:
        context.prec = 120
        a, b = number(a), number(b)
        return abs(a-b) <= tolerance * max(abs(a), abs(b), Decimal('1e-100'))


def decode(table, bits):
    names = [c['id'] for c in table['columns']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate column')
    rows = []
    for raw in table['rows']:
        if len(raw) != len(names):
            raise ValueError('Ragged row')
        values = []
        for value in raw:
            if isinstance(value, dict):
                if set(value) != {'value', 'precision'} or value['precision'] != bits:
                    raise ValueError('Wrong scalar precision')
                value = value['value']
            values.append(value)
        rows.append(dict(zip(names, values)))
    return rows


def capture(raw, precision, exit_code):
    tables = json.loads(raw)['tables']
    table = tables['spectrum']
    meta, summary = table['metadata'], table['summary']
    if tables['diagnostics']['metadata'] != meta or tables['diagnostics']['summary'] != summary:
        raise ValueError('Inconsistent table provenance')
    bits = int(meta['Geometry and output bits'])
    command = shlex.split(meta['Command'])
    command[0] = Path(command[0]).name
    return {'precision': precision, 'exit_code': exit_code,
            'metadata': {k: v for k, v in meta.items() if k != 'Command'},
            'command': command, 'summary': summary,
            'spectrum': decode(table, bits), 'diagnostics': decode(tables['diagnostics'], bits),
            'raw_sha256': hashlib.sha256(raw).hexdigest()}


def validate(data):
    if data['schema'] != 'bethe.four-spinon-precision.v1' or not data['cpu']:
        raise ValueError('Wrong schema or missing timing hardware')
    if [c['precision'] for c in data['cases']] != list(FORMATS):
        raise ValueError('Missing, duplicate or reordered precision case')
    sources = set()
    for case in data['cases']:
        meta, summary = case['metadata'], case['summary']
        if case['command'] != ['bethe-xxx-structure-factor-thermo', *ARGUMENTS, '--precision', case['precision']]:
            raise ValueError('Regeneration command differs from the fixed comparison')
        if (meta.get('Program') != 'bethe-xxx-structure-factor-thermo'
                or meta.get('Precision') != case['precision']
                or any(meta.get(k) != v for k, v in CONTROLS.items())):
            raise ValueError('Physical or integration controls changed with precision')
        bits = int(meta['Geometry and output bits'])
        expected_bits = {'fp64': 53, 'fp128': 113, 'bits-128': 128, 'bits-192': 192}
        if bits != expected_bits.get(case['precision'], bits) or bits < 53:
            raise ValueError('Incorrect arithmetic format')
        low = max(24, bits-64) if case['precision'].startswith('bits-') else (24 if bits == 53 else 53)
        if int(meta['Density comparison bits']) != low or int(meta['Density guard bits']) != bits:
            raise ValueError('Density comparison was changed')
        epsilon = Decimal(2) ** (8-bits)
        if not close(meta['Two-spinon transition-integral tolerance'], '1e-12', epsilon):
            raise ValueError('Two-spinon tolerance changed')
        for key in ('Bethe revision', 'Uni20 revision'):
            if not re.fullmatch('[0-9a-f]{40}', meta[key]):
                raise ValueError('Require clean source provenance')
        if not re.fullmatch('[0-9a-f]{64}', case['raw_sha256']):
            raise ValueError('Missing raw-source hash')
        sources.add((meta['Bethe revision'], meta['Uni20 revision'], meta['Compiler'], meta['Platform']))
        rows, diagnostics = case['spectrum'], case['diagnostics']
        if [r['sector'] for r in rows] != ['2', '4', 'sum:2+4'] or [r['sector'] for r in diagnostics] != ['2', '4']:
            raise ValueError('Incomplete selected-sector records')
        if any(not close(r['q'], '1.5', epsilon) or not close(r['omega'], '1.8', epsilon)
               for r in rows + diagnostics):
            raise ValueError('Wrong spectral point')
        if rows[0]['status'] != 'converged' or not close(rows[0]['density'], '0.87313045511634', Decimal('1e-11')):
            raise ValueError('Two-spinon reference failed')
        if rows[1]['status'] == 'estimated':
            if (case['exit_code'] != 0 or summary['Outcome'] != 'success'
                    or summary['Failed points'] != '0' or rows[2]['status'] != 'partial_sum'):
                raise ValueError('Inconsistent successful result')
            with localcontext() as context:
                context.prec = 120
                if not close(rows[2]['density'], number(rows[0]['density'])+number(rows[1]['density']), epsilon):
                    raise ValueError('Partial sum changed')
            if number(rows[1]['density']) <= 0:
                raise ValueError('Missing positive density')
            for key in ('contour_diagnostic', 'arithmetic_diagnostic'):
                if number(diagnostics[1][key]) < 0:
                    raise ValueError('Negative density diagnostic')
        else:
            if (case['exit_code'] != 2 or summary['Outcome'] != 'partial'
                    or summary['Failed points'] != '1' or rows[1]['status'] != 'numerical_failure'
                    or rows[1]['density'] is not None or rows[2]['density'] is not None
                    or rows[2]['status'] != 'incomplete' or not diagnostics[1]['reason']):
                raise ValueError('Failure must retain its reason and missing contributions')
        for key in ('Run CPU seconds', 'Elapsed seconds'):
            if number(summary[key]) <= 0:
                raise ValueError('Missing measured cost')
    if len(sources) != 1:
        raise ValueError('Mixed sources or toolchains')
    reference = data['cases'][-1]['spectrum'][1]
    if reference['status'] != 'estimated':
        raise ValueError('Missing complete highest-precision reference')
    for case in data['cases']:
        row = case['spectrum'][1]
        if row['status'] == 'estimated' and not close(row['density'], reference['density'], Decimal('1e-6')):
            raise ValueError('Arithmetic disagreement exceeds the demonstration budget')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable, help='Regenerate all five formats; includes expected failures')
    parser.add_argument('--raw-dir', type=Path, help='New export directory, or existing matrix to import')
    parser.add_argument('--cpu', help='Original CPU model; required when importing')
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--check', action='store_true', help='Validate the committed data without running a solver')
    args = parser.parse_args()
    if args.check and (args.solver or args.raw_dir):
        parser.error('--check does not run or import calculations')
    if args.solver and not args.raw_dir:
        parser.error('--solver requires a new --raw-dir to retain successes and failures')
    if args.raw_dir:
        if not args.solver and not args.cpu:
            parser.error('Import requires the original --cpu model')
        if args.solver:
            args.raw_dir.mkdir(parents=True, exist_ok=False)
        cases = []
        for precision in FORMATS:
            path = args.raw_dir / (precision+'.json')
            if args.solver:
                process = subprocess.run([str(args.solver), *ARGUMENTS, '--precision', precision],
                                         capture_output=True, timeout=600,
                                         env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never'))
                path.write_bytes(process.stdout)
                path.with_suffix('.err').write_bytes(process.stderr)
                path.with_suffix('.receipt.json').write_text(json.dumps({'exit_code': process.returncode})+'\n')
            receipt = json.loads(path.with_suffix('.receipt.json').read_text())
            cases.append(capture(path.read_bytes(), precision, receipt['exit_code']))
        cpu = args.cpu
        if not cpu and Path('/proc/cpuinfo').exists():
            cpu = next(s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name'))
        data = validate({'schema': 'bethe.four-spinon-precision.v1', 'cpu': cpu or platform.processor(), 'cases': cases})
        args.data.write_text(json.dumps(data, indent=2)+'\n')
    else:
        data = validate(json.loads(args.data.read_text()))
    for case in data['cases']:
        row = case['spectrum'][1]
        print(case['precision'], row['status'], row['density'], case['summary']['Run CPU seconds'])


if __name__ == '__main__':
    main()
