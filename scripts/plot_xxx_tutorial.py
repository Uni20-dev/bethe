#!/usr/bin/env python3
"""Finite XXX spinons and periodic/free-end real-root scans from CLI exports."""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
SIZES = (15, 31, 63)
STATE_COLUMNS = ['state_id', 'sz', 'spin_reversed', 'energy', 'gap', 'residual', 'iterations', 'converged', 'status']
MOMENTUM_COLUMNS = ['momentum_index', 'p']


def read_states(text, n, boundary, extra):
    expected = {'Program': f'bethe-xxx-{boundary}', 'Sites': str(n), 'Precision': 'fp64',
                'Model': 'periodic spin-1/2, J=1, h=0' if boundary == 'pbc' else 'open spin-1/2, free ends, J=1, h=0',
                **extra}
    columns = STATE_COLUMNS + (MOMENTUM_COLUMNS if boundary == 'pbc' else [])
    meta, rows = read_csv_export(text, columns, expected)
    tolerance = finite_number(meta, 'Residual tolerance')
    if tolerance <= 0:
        raise ValueError('Bad tolerance')
    parsed = []
    for i, row in enumerate(rows):
        numbers = {key: finite_number(row, key) for key in ('sz', 'energy', 'residual', 'iterations')}
        if (row['state_id'] != str(i) or row['converged'] != 'true' or row['spin_reversed'] != 'false'
                or not 0 <= numbers['residual'] <= tolerance or numbers['iterations'] < 0):
            raise ValueError('Invalid/unconverged state')
        if boundary == 'pbc':
            index = int(row['momentum_index']); p = finite_number(row, 'p')
            if not 0 <= index < n or not math.isclose(p, 2*math.pi*index/n, abs_tol=1e-14):
                raise ValueError('Bad lattice momentum')
            numbers.update(momentum_index=index, p=p)
        numbers['gap'] = None if row['gap'] == '' else finite_number(row, 'gap')
        parsed.append(numbers)
    return meta, parsed


def read_spinons(texts, n):
    expected = {'Bulk reference': 'e_inf = 1/4 - log(2)', 'Spinon momentum': 'k = pi/2 - 2*pi*hole/N',
                'Finite-size energy': 'E - N*e_inf retains finite-size corrections'}
    meta, states = read_states(texts['states'], n, 'pbc', expected)
    expected.update(Program='bethe-xxx-pbc', Sites=str(n), Precision='fp64')
    columns = ['state_id', 'hole', 'k', 'bulk_subtracted_energy', 'epsilon_inf']
    spin_meta, raw = read_csv_export(texts['spinons'], columns, expected, row_status=None)
    for key in ('Command', 'Bethe revision', 'Uni20 revision'):
        if spin_meta[key] != meta[key]:
            raise ValueError('Mismatched tables')
    if len(states) != (n+1)//2 or len(raw) != len(states):
        raise ValueError('Incomplete branch')
    parsed = []
    for i, (row, state) in enumerate(zip(raw, states)):
        numbers = {key: finite_number(row, key) for key in columns[1:]}
        if (row['state_id'] != str(i) or state['sz'] != .5 or state['gap'] is not None
                or numbers['hole'] != (n-1)/4-i
                or not math.isclose(numbers['k'], math.pi/2-2*math.pi*numbers['hole']/n, abs_tol=1e-14)
                or not math.isclose(numbers['bulk_subtracted_energy'], state['energy']-n*(.25-math.log(2)), abs_tol=1e-12)):
            raise ValueError('Invalid spinon coordinates/reference')
        parsed.append(dict(state, **numbers))
    return parsed


def read_scan(text, n, boundary):
    count = math.comb(n//2+1, n//2-1)
    expected = {'Family': 'restricted real-root highest-weight multiplets; NOT a complete spectrum',
                'S': '1', 'Multiplet size': '3', 'Candidates': str(count), 'Converged candidates': str(count),
                'Returned multiplets': str(count), 'Ordering': 'complete within supported family',
                'Ground converged': '1', 'Ground status': 'converged', 'Gap reference': 'E-E0; global ground state'}
    meta, rows = read_states(text, n, boundary, expected)
    ground = finite_number(meta, 'Ground energy')
    if len(rows) != count:
        raise ValueError('Incomplete scan')
    for row in rows:
        if (row['sz'] != 1 or row['gap'] is None or row['gap'] < 0
                or not math.isclose(row['gap'], row['energy']-ground, abs_tol=1e-13)):
            raise ValueError('Bad gap reference or spin')
    if any(a['energy'] > b['energy'] for a, b in zip(rows, rows[1:])):
        raise ValueError('Unsorted scan')
    return ground, rows


def load_cases(build_dir=None):
    exports, cases = {}, {}
    for n in SIZES:
        texts = (capture_csv_tables([build_dir / 'bethe-xxx-pbc', n, '--spinons'], ('states', 'spinons'))
                 if build_dir else {name: (DATA / f'xxx-n{n}-{name}.csv').read_text() for name in ('states', 'spinons')})
        cases[f'spinons-{n}'] = read_spinons(texts, n)
        exports.update({f'xxx-n{n}-{name}': text for name, text in texts.items()})
    for n in (4, 12):
        for boundary in ('pbc', 'obc'):
            name = f'xxx-n{n}-{boundary}'
            text = (subprocess.run([str(build_dir / f'bethe-xxx-{boundary}'), str(n), '--excitations', 'all',
                                    '--spin', '1', '--precision', 'fp64', '--format', 'csv'],
                                   check=True, capture_output=True, text=True).stdout
                    if build_dir else (DATA / f'{name}.csv').read_text())
            cases[f'{boundary}-{n}'] = read_scan(text, n, boundary)
            exports[name] = text
    if build_dir:
        for name, text in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-xxx-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), layout='constrained')
    rows = cases['spinons-63']
    axes[0].plot([r['k']/math.pi for r in rows], [r['epsilon_inf'] for r in rows], ':', color='0.3', label='Thermodynamic')
    for n, color, marker in zip(SIZES, ('#0072B2', '#D55E00', '#009E73'), ('o', 's', '^')):
        rows = cases[f'spinons-{n}']
        axes[0].plot([r['k']/math.pi for r in rows], [r['bulk_subtracted_energy'] for r in rows],
                     marker, color=color, markersize=4, label=f'N={n}')
        axes[1].plot([r['k']/math.pi for r in rows], [r['bulk_subtracted_energy']-r['epsilon_inf'] for r in rows],
                     marker+'-', color=color, markersize=3, linewidth=1, label=f'N={n}')
    axes[0].set(title='Odd-ring one-spinon branch', ylabel=r'$E_N-N e_\infty$')
    axes[1].set(title='Finite-size correction', ylabel=r'$E_N-N e_\infty-\epsilon_\infty(k)$')
    for ax in axes:
        ax.set(xlabel=r'Spinon momentum $k/\pi$', xlim=(0, 1))
        ax.grid(alpha=.15)
        ax.legend(fontsize=8)
    save_svg(fig, FIGURES / 'xxx-finite-spinons.svg')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    periodic = cases['pbc-12'][1]
    axes[0].scatter([r['p']/math.pi for r in periodic], [r['gap'] for r in periodic], color='#0072B2', s=30)
    axes[0].set(title='Periodic: momentum-resolved', xlabel=r'Lattice momentum $P/\pi$', ylabel=r'$E-E_0$',
                xlim=(0, 2), xticks=[0, .5, 1, 1.5, 2])
    for boundary, color, marker in [('pbc', '#0072B2', 'o'), ('obc', '#D55E00', 's')]:
        rows = cases[f'{boundary}-12'][1]
        axes[1].plot(range(1, len(rows)+1), [r['gap'] for r in rows], marker+'-', color=color,
                     markersize=4, linewidth=1, label='Periodic' if boundary == 'pbc' else 'Free ends')
    axes[1].set(title='Same family size, different boundaries', xlabel='Energy rank (not momentum)', ylabel=r'$E-E_0$',
                xticks=[1, 5, 10, 15, 21])
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.15)
    fig.suptitle('N=12, S=1: supported real-root multiplets only', fontsize=12)
    save_svg(fig, FIGURES / 'xxx-boundary-scans.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.build_dir.resolve() if args.build_dir else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} XXX calculations.')
