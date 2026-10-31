#!/usr/bin/env python3
"""Reproduce a checked order/cutoff study of one physical four-spinon point."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import platform
import subprocess

from tutorial_common import save_svg, solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data/xxx-four-spinon-convergence.json'
FIGURE = ROOT / 'docs/tutorials/figures/xxx-four-spinon-convergence.svg'
ORDERS, CUTOFFS = (2, 4, 8, 16), (8, 12)
CONTROL = {'momentum': '1.5', 'omega': '1.8', 'exchange': '1', 'channel': 'zz',
           'precision': 'bits-128', 'relative_tolerance': '1e-8', 'absolute_tolerance': '1e-30',
           'two_kernel_tolerance': '1e-12'}
BITS = 128
ROUNDING_ALLOWANCE = Decimal(2) ** (8 - BITS)


def command(solver, order, cutoff, threads):
    return [str(solver), '--sectors', '2,4', '--momentum', CONTROL['momentum'],
            '--omega', CONTROL['omega'], '--exchange', CONTROL['exchange'],
            '--channel', CONTROL['channel'], '--precision', CONTROL['precision'],
            '--four-order', str(order), '--four-cutoff', str(cutoff),
            '--density-relative-tolerance', CONTROL['relative_tolerance'],
            '--density-absolute-tolerance', CONTROL['absolute_tolerance'],
            '--density-terms', '64', '--density-tail-order', '16',
            '--tolerance', CONTROL['two_kernel_tolerance'],
            '--threads', str(threads), '--diagnostics', '--format', 'json']


def number(value):
    result = Decimal(str(value))
    if not result.is_finite():
        raise ValueError('Nonfinite numerical value')
    return result


def rows(table):
    columns = [column['id'] for column in table['columns']]
    result = []
    for row in table['rows']:
        if len(row) != len(columns):
            raise ValueError('Ragged table')
        values = []
        for cell in row:
            if isinstance(cell, dict):
                if cell.get('precision') != BITS:
                    raise ValueError('Unexpected output precision')
                cell = cell['value']
            values.append(cell)
        result.append(dict(zip(columns, values)))
    return result


def read_case(document, order, cutoff):
    spectrum, diagnostics = (document['tables'][name] for name in ('spectrum', 'diagnostics'))
    meta, summary = spectrum['metadata'], spectrum['summary']
    expected = {'Program': 'bethe-xxx-structure-factor-thermo', 'Precision': CONTROL['precision'],
                'Exchange J': '1', 'Channel': 'zz', 'Selected spinon sectors': '2,4',
                'Spectral points': '1', 'Four-spinon order': str(order),
                'Four-spinon log-radius cutoff': str(cutoff), 'Geometry and output bits': str(BITS),
                'Density comparison bits': str(BITS-64), 'Density guard bits': str(BITS),
                'Density relative tolerance': '1e-8', 'Density absolute tolerance': '1e-30',
                'First density residue budget': '64', 'Second density residue budget': '80',
                'First density tail order': '16', 'Second density tail order': '20',
                'Density convention': 'S(q,w); moments use dw/(2*pi); raising=2*zz'}
    if any(meta.get(k) != v for k, v in expected.items()):
        raise ValueError('Unexpected model, normalization or numerical controls')
    if abs(number(meta['Two-spinon transition-integral tolerance']) / Decimal(CONTROL['two_kernel_tolerance']) - 1) > ROUNDING_ALLOWANCE:
        raise ValueError('Two-spinon tolerance changed')
    if meta != diagnostics['metadata'] or summary != diagnostics['summary']:
        raise ValueError('Mismatched table provenance')
    if summary.get('Outcome') != 'success' or summary.get('Failed points') != '0':
        raise ValueError('Incomplete calculation')
    values, checks = rows(spectrum), rows(diagnostics)
    if [r['sector'] for r in values] != ['2', '4', 'sum:2+4'] or [r['sector'] for r in checks] != ['2', '4']:
        raise ValueError('Missing or duplicate selected contribution')
    if [r['status'] for r in values] != ['converged', 'estimated', 'partial_sum']:
        raise ValueError('Failed or incorrectly labelled contribution')
    for row in values + checks:
        if abs(number(row['q']) - Decimal('1.5')) > ROUNDING_ALLOWANCE or abs(number(row['omega']) - Decimal('1.8')) > ROUNDING_ALLOWANCE:
            raise ValueError('Wrong spectral point')
    result = {'order': order, 'cutoff': cutoff, 's2': values[0]['density'], 's4': values[1]['density'],
              'sum': values[2]['density'], 'two_kernel_error': checks[0]['two_kernel_error'],
              'contour_diagnostic': checks[1]['contour_diagnostic'],
              'arithmetic_diagnostic': checks[1]['arithmetic_diagnostic'],
              'density_attempts': int(values[1]['evaluations']),
              'wall_seconds': summary['Elapsed seconds'], 'cpu_seconds': summary['Run CPU seconds'],
              'threads': int(meta['Scheduler concurrency limit']),
              'source_commit': meta['Bethe revision'], 'uni20_commit': meta['Uni20 revision'],
              'compiler': meta['Compiler'], 'platform': meta['Platform']}
    validate_case(result)
    return result


def validate_case(case):
    with localcontext() as context:
        context.prec = 100
        if any(number(case[key]) <= 0 for key in ('s2', 's4', 'sum', 'cpu_seconds', 'wall_seconds')):
            raise ValueError('Expected positive interior response and measured cost')
        if abs(number(case['sum']) - number(case['s2']) - number(case['s4'])) > ROUNDING_ALLOWANCE:
            raise ValueError('Partial sum differs from physical sector sum')
        if any(number(case[key]) < 0 for key in ('two_kernel_error', 'contour_diagnostic', 'arithmetic_diagnostic')):
            raise ValueError('Negative numerical diagnostic')
    if case['density_attempts'] <= 0 or case['threads'] <= 0:
        raise ValueError('Missing numerical work')
    for key in ('source_commit', 'uni20_commit'):
        value = case[key]
        if len(value) != 40 or any(c not in '0123456789abcdef' for c in value):
            raise ValueError('Require an identifiable clean source revision')


def validate(data):
    if data['schema'] != 'bethe.four-spinon-convergence.v1' or data['controls'] != CONTROL:
        raise ValueError('Unexpected tutorial data contract')
    expected = {(n, t) for n in ORDERS for t in CUTOFFS}
    keys = [(row['order'], row['cutoff']) for row in data['cases']]
    if set(keys) != expected or len(keys) != len(expected):
        raise ValueError('Incomplete or duplicate order/cutoff grid')
    for row in data['cases']:
        validate_case(row)
    for key in ('source_commit', 'uni20_commit', 'compiler', 'platform', 'threads'):
        if len({row[key] for row in data['cases']}) != 1:
            raise ValueError(f'Mixed {key} in the refinement study')
    if not data['cpu']:
        raise ValueError('Missing timing hardware')
    reference = number(data['cases'][0]['s2'])
    if abs(reference - Decimal('0.87313045511634')) > Decimal('1e-13'):
        raise ValueError('Two-spinon reference normalization changed')
    if any(abs(number(row['s2']) - reference) > ROUNDING_ALLOWANCE for row in data['cases']):
        raise ValueError('Two-spinon baseline changed during four-spinon refinement')
    return data


def filename(order, cutoff):
    return f'q1p5-w1p8-order{order}-cutoff{cutoff}-bits{BITS}.json'


def machine_name():
    if Path('/proc/cpuinfo').is_file():
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.startswith('model name'):
                return line.split(':', 1)[1].strip()
    return platform.processor() or platform.machine()


def collect(solver, raw, threads, machine):
    cases = []
    for cutoff in CUTOFFS:
        for order in ORDERS:
            path = raw / filename(order, cutoff) if raw else None
            if solver:
                if path and path.exists():
                    raise FileExistsError(f'Choose a new raw directory; refusing to overwrite {path}')
                result = subprocess.run(command(solver, order, cutoff, threads), capture_output=True,
                                        text=True, check=False, timeout=3600)
                if path:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(result.stdout)
                    path.with_suffix('.err').write_text(result.stderr)
                result.check_returncode()
                document = json.loads(result.stdout)
                case = read_case(document, order, cutoff)
                case['raw_sha256'] = hashlib.sha256(result.stdout.encode()).hexdigest()
            else:
                raw_text = path.read_text()
                document = json.loads(raw_text)
                case = read_case(document, order, cutoff)
                case['raw_sha256'] = hashlib.sha256(raw_text.encode()).hexdigest()
            cases.append(case)
    return validate({'schema': 'bethe.four-spinon-convergence.v1', 'controls': CONTROL,
                     'cpu': machine, 'timing_scope': 'One two-sector spectral query including geometry and output; shared host.',
                     'cases': cases})


def plot(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-four-spinon-convergence', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    for cutoff, color in zip(CUTOFFS, ('#0072B2', '#D55E00')):
        points = sorted((x for x in data['cases'] if x['cutoff'] == cutoff), key=lambda x: x['order'])
        x = [r['order'] for r in points]
        axes[0].plot(x, [float(r['s4']) for r in points], 'o-', color=color, label=f'Soft cutoff T={cutoff}')
        axes[1].plot(x, [float(r['cpu_seconds']) for r in points], 'o-', color=color, label=f'T={cutoff}')
    for axis in axes:
        axis.set_xscale('log', base=2)
        axis.set_xticks(ORDERS, [str(x) for x in ORDERS])
        axis.set_xlabel('Gauss order per panel')
        axis.grid(alpha=.2)
        axis.legend(frameon=False)
    axes[0].set(ylabel=r'$S^{zz}_4(q,\omega)$', title=r'Four spinons at $q=1.5$, $\omega=1.8J$')
    axes[1].set_yscale('log')
    workers = data['cases'][0]['threads']
    axes[1].set(ylabel='CPU seconds', title=f'Measured cost; {workers} workers per query')
    save_svg(fig, FIGURE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable, help='Regenerate all eight points; can take tens of minutes')
    parser.add_argument('--raw-dir', type=Path, help='Read existing JSON exports, or retain them when regenerating')
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--cpu', help='CPU model for imported exports; required with --raw-dir alone')
    parser.add_argument('--check', action='store_true', help='Validate data without plotting')
    args = parser.parse_args()
    if args.threads <= 0:
        parser.error('--threads must be positive')
    if args.raw_dir and not args.solver and not args.cpu:
        parser.error('imported timings need their original CPU model via --cpu')
    if args.solver or args.raw_dir:
        data = collect(args.solver, args.raw_dir, args.threads, args.cpu or machine_name())
        args.data.parent.mkdir(parents=True, exist_ok=True)
        args.data.write_text(json.dumps(data, indent=2) + '\n')
    else:
        data = validate(json.loads(args.data.read_text()))
    if not args.check:
        plot(data)
    print(f'Checked {len(data["cases"])} physical order/cutoff points')


if __name__ == '__main__':
    main()
