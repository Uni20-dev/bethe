#!/usr/bin/env python3
"""Validate or reproduce the small API examples and native binary128 refinement grid."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

from check_xxx_four_spinon_precision import close, decode, number
from tutorial_common import solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/validation/spinon-native-precision-20261031.json'
REFERENCE = ROOT / 'docs/validation/four-spinon-refinement-precision-20261031.json'
FORMATS = ('fp64', 'long-double', 'fp128', 'bits-128', 'bits-192', 'bits-384')
Q = '1.57079632679489661923132169163975'


def queries():
    common = ['--momentum', Q, '--sectors', '6', '--six-cutoff', '8', '--seed', '42',
              '--threads', '2', '--diagnostics', '--replicate-output']
    for kind, extra in [('spectrum', ['--omega', '2', '--replicates', '3', '--samples', '8']),
                        ('moments', ['--replicates', '2', '--samples', '4', '--covariance'])]:
        for precision in FORMATS:
            yield f'{kind}-{precision}', kind, [*common, *extra, '--precision', precision]
    for cutoff in (8, 12):
        for order in (2, 4, 8, 16):
            yield f'four-fp128-n{order}-T{cutoff}', 'spectrum', [
                '--sectors', '2,4', '--momentum', '1.5', '--omega', '1.8', '--precision', 'fp128',
                '--four-order', str(order), '--four-cutoff', str(cutoff),
                '--density-relative-tolerance', '1e-8', '--density-absolute-tolerance', '1e-30',
                '--density-terms', '64', '--density-tail-order', '16', '--tolerance', '1e-12',
                '--threads', '2', '--diagnostics']


def capture(path, name, kind, exit_code):
    raw = path.read_bytes()
    tables = json.loads(raw)['tables']
    meta, summary = tables[kind]['metadata'], tables[kind]['summary']
    if any(t['metadata'] != meta or t['summary'] != summary for t in tables.values()):
        raise ValueError('Inconsistent table provenance')
    command = shlex.split(meta['Command'])
    command[0] = Path(command[0]).name
    # Retain reproducible numerical arguments, without a machine-specific export path.
    i = command.index('--json')
    del command[i:i+2]
    command.remove('--quiet')
    bits = int(meta['Working precision (bits)'])
    return {'name': name, 'exit_code': exit_code, 'command': command,
            'metadata': {k: v for k, v in meta.items() if k != 'Command'}, 'summary': summary,
            'tables': {k: decode(t, bits) for k, t in tables.items()},
            'raw_sha256': hashlib.sha256(raw).hexdigest()}


def validate(data):
    if data['schema'] != 'bethe.native-density-precision.v1' or not data['cpu']:
        raise ValueError('Missing schema or hardware')
    cases = data['cases']
    expected = list(queries())
    if [c['name'] for c in cases] != [q[0] for q in expected]:
        raise ValueError('Missing, duplicated or reordered query')
    native_bits = int(cases[1]['metadata']['Working precision (bits)'])
    sources = set()
    for case, (name, kind, arguments) in zip(cases, expected):
        program = 'bethe-xxx-structure-factor-thermo' if kind == 'spectrum' else 'bethe-xxx-moments-thermo'
        if case['command'] != [program, *arguments]:
            raise ValueError('Changed physical or sampling command')
        meta, summary, tables = case['metadata'], case['summary'], case['tables']
        precision = arguments[arguments.index('--precision')+1]
        bits = int(meta['Working precision (bits)'])
        wanted_bits = {'fp64': 53, 'long-double': native_bits, 'fp128': 113,
                       'bits-128': 128, 'bits-192': 192, 'bits-384': 384}[precision]
        low = max(24, bits-64) if precision.startswith('bits-') else 24 if bits <= 53 else 53
        if precision == 'fp128' and 53 < native_bits < 113:
            low = native_bits
        if (bits != wanted_bits or meta['Precision'] != precision
                or int(meta['Density comparison bits']) != low or int(meta['Density guard bits']) != bits):
            raise ValueError('Incorrect native or MPFR comparison pair')
        four = name.startswith('four-')
        controls = {'Channel': 'zz', 'Exchange J': '1', 'Scheduler concurrency limit': '2',
                    'Density relative tolerance': '1e-8' if four else '1e-6',
                    'Density absolute tolerance': '1e-30' if four else '0',
                    'First density tail order': '16', 'Second density tail order': '20'}
        term_name = 'density residue budget' if kind == 'spectrum' else 'density terms'
        controls.update({f'First {term_name}': '64', f'Second {term_name}': '80'})
        if four:
            controls.update({'Four-spinon order': arguments[arguments.index('--four-order')+1],
                             'Four-spinon log-radius cutoff': arguments[arguments.index('--four-cutoff')+1]})
        else:
            controls.update({'Six-spinon absolute rapidity cutoff': '8', 'Six-spinon seed': '42',
                             'Six-spinon replicates': '3' if kind == 'spectrum' else '2',
                             'Samples per replicate': '8' if kind == 'spectrum' else '4'})
        if any(meta.get(k) != v for k, v in controls.items()):
            raise ValueError('Changed density controls')
        for key in ('Bethe revision', 'Uni20 revision'):
            if not re.fullmatch('[0-9a-f]{40}', meta[key]):
                raise ValueError('Missing clean source revision')
        sources.add((meta['Bethe revision'], meta['Uni20 revision'], meta['Compiler'], meta['Platform']))
        if not re.fullmatch('[0-9a-f]{64}', case['raw_sha256']):
            raise ValueError('Missing raw hash')
        if any(number(summary[k]) <= 0 for k in ('Run CPU seconds', 'Elapsed seconds')):
            raise ValueError('Missing measured cost')
        value_key = 'density' if kind == 'spectrum' else 'value'
        success = case['exit_code'] == 0
        if (case['exit_code'] not in (0, 2)
                or summary['Outcome'] != ('success' if success else 'partial')):
            raise ValueError('Inconsistent outcome')
        values = tables[kind]
        if not values or not tables['diagnostics']:
            raise ValueError('Missing results or diagnostics')
        sectors = ['2', '4', 'sum:2+4'] if four else ['6'] if kind == 'spectrum' else ['6', '6']
        if [r['sector'] for r in values] != sectors:
            raise ValueError('Changed sector identity')
        if any(not close(r['q'], '1.5' if four else Q, Decimal('1e-15')) for r in values):
            raise ValueError('Wrong momentum')
        if kind == 'spectrum' and any(not close(r['omega'], '1.8' if four else '2', Decimal('1e-15'))
                                      for r in values):
            raise ValueError('Wrong frequency')
        if success:
            if any(number(r[value_key]) <= 0 for r in values):
                raise ValueError('Missing successful value')
        elif (any(r[value_key] is not None for r in values)
              or any(r['status'] != 'numerical_failure' for r in values)
              or not any(r['reason'] for r in tables['diagnostics'])):
            raise ValueError('Failed result cannot become a partial estimate')
        if not four:
            count, samples = (3, 8) if kind == 'spectrum' else (2, 4)
            if [int(r['replicate']) for r in tables['replicates']] != list(range(count)):
                raise ValueError('Missing replicate')
            if success and any(int(r['completed_samples']) != samples or r['failed_sample'] is not None
                               for r in tables['replicates']):
                raise ValueError('Incomplete successful replicate')
            if kind == 'moments' and [(r['moment_i'], r['moment_j']) for r in tables['covariance']] != [
                    (0, 0), (0, 1), (1, 0), (1, 1)]:
                raise ValueError('Missing moment covariance')
            if not success and kind == 'moments' and any(r['covariance'] is not None for r in tables['covariance']):
                raise ValueError('Failed moment covariance must be absent')
    if len(sources) != 1:
        raise ValueError('Mixed calculation sources')
    return data


def comparisons(data):
    validate(data)
    results = {}
    fields = {'spectrum': {'spectrum': ['density'], 'diagnostics': ['sampling_error'],
                          'replicates': ['density', 'maximum_absolute_density']},
              'moments': {'moments': ['value', 'fraction'], 'diagnostics': ['sampling_error'],
                          'replicates': ['m0', 'm1', 'maximum_m0', 'maximum_m1'],
                          'covariance': ['covariance']}}
    cases = {c['name']: c for c in data['cases']}
    with localcontext() as context:
        context.prec = 160
        for kind, spec in fields.items():
            reference = cases[f'{kind}-bits-384']
            if reference['exit_code'] != 0:
                raise ValueError('Missing wider reference')
            for precision in FORMATS:
                case = cases[f'{kind}-{precision}']
                if case['exit_code'] != 0:
                    continue
                differences = []
                for table, columns in spec.items():
                    a, b = case['tables'][table], reference['tables'][table]
                    if len(a) != len(b):
                        raise ValueError('Changed observations')
                    for r, s in zip(a, b):
                        for column in columns:
                            if r[column] is None or s[column] is None:
                                if r[column] is not s[column]:
                                    raise ValueError('Missing observable')
                                continue
                            x, y = number(r[column]), number(s[column])
                            differences.append(abs(x-y)/max(abs(x), abs(y), Decimal('1e-150')))
                for r, s in zip(case['tables']['replicates'], reference['tables']['replicates']):
                    for k in ('completed_samples', 'supported_roots', 'supported_configurations', 'density_attempts'):
                        if k in r and r[k] != s[k]:
                            raise ValueError('Sampling stream counts changed')
                maximum = max(differences)
                if maximum > Decimal('1e-12'):
                    raise ValueError('Disagreement with the wider reference')
                results[case['name']] = str(maximum)
        historical = json.loads(REFERENCE.read_text())
        for case in data['cases'][12:]:
            meta = case['metadata']
            order, cutoff = int(meta['Four-spinon order']), int(meta['Four-spinon log-radius cutoff'])
            reference, = [c for c in historical['reference_192']['cases']
                           if c['order'] == order and c['cutoff'] == cutoff]
            if case['exit_code'] != 0:
                raise ValueError('Incomplete native refinement grid')
            value = number(case['tables']['spectrum'][1]['density'])
            difference = abs(value-number(reference['s4']))/abs(number(reference['s4']))
            if difference > Decimal('1e-28'):
                raise ValueError('Native grid disagrees with MPFR192')
            results[case['name']] = str(difference)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, help='Import retained exports, or create fresh ones with both solvers')
    parser.add_argument('--solver', type=solver_executable)
    parser.add_argument('--moments', type=solver_executable)
    parser.add_argument('--cpu', help='Original CPU model, required for import')
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if bool(args.solver) != bool(args.moments) or (args.solver and not args.raw_dir):
        parser.error('Regeneration needs both solvers and a new --raw-dir')
    if args.check and (args.raw_dir or args.solver):
        parser.error('--check only reads the compact record')
    if args.raw_dir:
        if args.solver:
            args.raw_dir.mkdir(parents=True, exist_ok=False)
        elif not args.cpu:
            parser.error('Import needs the original --cpu')
        cases = []
        receipt = args.raw_dir/'receipt.json'
        exits = {r['case']: r['exit_code'] for r in json.loads(receipt.read_text())['runs']} if receipt.exists() else {}
        for name, kind, arguments in queries():
            path = args.raw_dir/(name+'.json')
            if args.solver:
                program = args.solver if kind == 'spectrum' else args.moments
                command = [str(program), *arguments, '--json', str(path), '--quiet']
                process = subprocess.run(command, capture_output=True, timeout=1200,
                    env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never'))
                path.with_suffix('.stderr').write_bytes(process.stderr)
                path.with_suffix('.receipt.json').write_text(json.dumps({'exit_code': process.returncode})+'\n')
                exits[name] = process.returncode
            if name not in exits:
                exits[name] = json.loads(path.with_suffix('.receipt.json').read_text())['exit_code']
            cases.append(capture(path, name, kind, exits[name]))
        cpu = args.cpu
        if not cpu and Path('/proc/cpuinfo').exists():
            cpu = next(s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                       if s.startswith('model name'))
        data = {'schema': 'bethe.native-density-precision.v1', 'cpu': cpu,
                'scope': 'Exact small API commands and fixed order/cutoff grid; no omitted-tail or total-error bound.',
                'cases': cases}
        data['maximum_relative_differences'] = comparisons(data)
        args.data.write_text(json.dumps(data, indent=2)+'\n')
    else:
        data = json.loads(args.data.read_text())
    result = comparisons(data)
    if result != data['maximum_relative_differences']:
        raise ValueError('Stored comparison summary differs')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
