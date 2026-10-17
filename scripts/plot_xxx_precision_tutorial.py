#!/usr/bin/env python3
"""One-magnon precision study from XXX CLI exports, with an independent oracle."""

import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import json
import os
from pathlib import Path
import subprocess

from tutorial_common import save_svg, solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data/xxx-precision.jsonl'
FIGURE = ROOT / 'docs/tutorials/figures/xxx-precision.svg'
SIZES = tuple(10**power for power in range(2, 17))
PRECISIONS = ('fp64', 'long-double', 'fp128', 'bits-128', 'bits-256')
DECIMAL_DIGITS = 110  # More than the largest working precision (256 bits).


@lru_cache(maxsize=1)
def reference_pi():
    """Machin's formula, evaluated independently of the solver/math library."""
    with localcontext() as ctx:
        ctx.prec = DECIMAL_DIGITS

        def atan_inverse(n):
            x = Decimal(1) / n
            power, total, k = x, x, 1
            while True:
                power *= -x*x
                updated = total + power / (2*k+1)
                if updated == total:
                    return total
                total, k = updated, k+1

        return 16*atan_inverse(5) - 4*atan_inverse(239)


@lru_cache(maxsize=None)
def reference_shift(n):
    """2 sin^2(pi/N), using a Decimal Taylor series, not 1-cos(2*pi/N)."""
    with localcontext() as ctx:
        ctx.prec = DECIMAL_DIGITS
        x = reference_pi() / n
        term, total, k = x, x, 1
        while True:
            term *= -x*x / ((2*k)*(2*k+1))
            updated = total + term
            if updated == total:
                return 2*total*total
            total, k = updated, k+1


def value(cell):
    """Preserve arbitrary-precision decimal strings; never pass through float."""
    return cell['value'] if isinstance(cell, dict) else cell


def one_row(table):
    columns = [column['id'] for column in table['columns']]
    if len(table['rows']) != 1 or len(table['rows'][0]) != len(columns):
        raise ValueError('Expected one state and one root')
    return dict(zip(columns, table['rows'][0]))


def capture(solver, n, precision):
    command = [str(solver), str(n), '--quantum-numbers', str(n//2-1),
               '--precision', precision, '--roots', '--format', 'json']
    process = subprocess.run(command, check=True, capture_output=True, text=True,
                             timeout=30, env=dict(os.environ, OPENBLAS_NUM_THREADS='1',
                                                 OMP_NUM_THREADS='1', UNI20_COLOR='never'))
    document = json.loads(process.stdout)
    if document['status'] != 'complete':
        raise ValueError('Incomplete JSON export')
    states, roots = (document['tables'][name] for name in ('states', 'roots'))
    if states['metadata'] != roots['metadata'] or states['summary'] != roots['summary']:
        raise ValueError('Mismatched state/root provenance')
    # Compact the export, retaining numeric strings, precision tags and provenance.
    return dict(sites=n, precision=precision, metadata=states['metadata'],
                summary=states['summary'], state=one_row(states), root=one_row(roots))


def analyze(case):
    n, precision = case['sites'], case['precision']
    meta, state, root = case['metadata'], case['state'], case['root']
    expected = {'Program': 'bethe-xxx-pbc', 'Sites': str(n), 'Precision': precision,
                'Model': 'periodic spin-1/2, J=1, h=0', 'Momentum index': '1'}
    if any(meta.get(key) != val for key, val in expected.items()):
        raise ValueError('Wrong state or precision metadata')
    if any(not meta.get(key) for key in ('Command', 'Bethe revision', 'Uni20 revision',
                                        'Compiler', 'Platform', 'Arithmetic format')):
        raise ValueError('Missing provenance')
    if case['summary']['Outcome'] != 'success' or case['summary']['Status'] != 'converged':
        raise ValueError('Failed solver run')
    if (state['state_id'] != '0' or root['state_id'] != '0' or root['index'] != '0'
            or state['momentum_index'] != '1' or state['sz'] != n//2-1
            or root['quantum_number'] != n//2-1 or state['spin_reversed']
            or state['converged'] is not True or state['status'] != 'converged'):
        raise ValueError('Not the specified one-magnon state')
    bits = int(meta['Working precision (bits)'])
    requested = {'fp64': 53, 'fp128': 113, 'bits-128': 128, 'bits-256': 256}
    if bits != requested.get(precision, bits) or bits < 53:
        raise ValueError('Wrong working precision')
    for cell in (*state.values(), *root.values()):
        if isinstance(cell, dict) and cell.get('precision') != bits:
            raise ValueError('Mismatched per-value precision')
    with localcontext() as ctx:
        ctx.prec = DECIMAL_DIGITS
        energy = Decimal(value(state['energy']))
        z = Decimal(value(root['rapidity']))
        residual = Decimal(value(state['residual']))
        tolerance = Decimal(meta['Residual tolerance'])
        if (not all(x.is_finite() for x in (energy, z, residual, tolerance))
                or z <= 0 or tolerance <= 0 or not 0 <= residual <= tolerance):
            raise ValueError('Invalid numerical output')
        exact = reference_shift(n)
        subtracted = Decimal(n)/4 - energy
        from_root = 2/(1+z*z)
        errors = [abs(estimate/exact-1) for estimate in (subtracted, from_root)]
        # These are loose conditioning-aware sanity bounds, not accuracy claims.
        epsilon = Decimal(2)**(1-bits)
        for error, amplification in zip(errors, (n**3, n)):
            if error > 1000*epsilon*amplification:
                raise ValueError('Error exceeds the conditioning-aware check')
        return dict(sites=n, precision=precision, bits=bits, shift=subtracted,
                    from_root=from_root, exact=exact, residual=residual,
                    energy_error=errors[0], root_error=errors[1])


def load_cases(solver=None):
    cases = ([capture(solver, n, precision) for precision in PRECISIONS for n in SIZES]
             if solver else [json.loads(line) for line in DATA.read_text().splitlines()])
    keys = [(case['sites'], case['precision']) for case in cases]
    if len(keys) != len(set(keys)) or set(keys) != {(n, p) for n in SIZES for p in PRECISIONS}:
        raise ValueError('Incomplete or duplicated precision study')
    results = [analyze(case) for case in cases]
    if solver:
        # Validate the whole replacement before overwriting the checked-in data.
        DATA.write_text(''.join(json.dumps(case, ensure_ascii=False) + '\n' for case in cases))
    return results


def plot(results):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-xxx-precision', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout='constrained')
    styles = zip(PRECISIONS, ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#333333'),
                 ('o', 's', '^', 'D', 'v'))
    for precision, color, marker in styles:
        rows = sorted((row for row in results if row['precision'] == precision), key=lambda r: r['sites'])
        label = f"{precision} ({rows[0]['bits']} bits)"
        for ax, key in zip(axes, ('energy_error', 'root_error')):
            ax.loglog([row['sites'] for row in rows], [float(row[key]) for row in rows],
                      marker+'-', color=color, label=label, markersize=3.5, linewidth=1)
    axes[0].set_title(r'From exported total energy: $N/4-E$')
    axes[1].set_title(r'From exported root: $2/(1+z^2)$')
    for ax in axes:
        ax.axhline(1, color='0.6', linestyle=':', linewidth=1)
        ax.set(xlabel='Ring length N (one root at every size)', ylabel='Relative error in energy shift',
               xlim=(70, 2e16), ylim=(1e-80, 10), xticks=[1e2, 1e6, 1e10, 1e14],
               yticks=[1e-75, 1e-60, 1e-45, 1e-30, 1e-15, 1])
        ax.grid(alpha=.15)
    axes[0].legend(fontsize=8, loc='lower right')
    fig.suptitle('A zero Bethe residual does not certify the excitation energy', fontsize=12)
    save_svg(fig, FIGURE)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable,
                        help='Regenerate using all five formats (requires MPREAL and FLOAT128)')
    parser.add_argument('--check', action='store_true', help='Validate without plotting')
    args = parser.parse_args()
    results = load_cases(args.solver)
    if not args.check:
        plot(results)
    for row in results:
        if row['sites'] in (10**6, 10**12, 10**16):
            print(f"N={row['sites']} {row['precision']:12s} residual={row['residual']} "
                  f"subtraction_error={row['energy_error']:.4E} root_error={row['root_error']:.4E}")
    print(f'Validated {len(results)} one-magnon calculations against the analytic dispersion.')
