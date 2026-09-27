#!/usr/bin/env python3
"""Finite Sutherland label windows; use --solver to regenerate exact-rule exports."""

import argparse
import itertools
import math
from pathlib import Path
import subprocess

from tutorial_common import finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
CASES = {'sutherland-lambda0': (0, 4, 1), 'sutherland-lambda1': (1, 4, 1),
         'sutherland-lambda2': (2, 4, 1), 'sutherland-window2': (1, 4, 2),
         'sutherland-length8': (1, 8, 1)}


def read_export(text, exponent, length, window):
    count = math.comb(4+2*window, 4)
    columns = ['state_id', 'labels', 'energy', 'gap', 'momentum_index', 'p']
    expected = {'Program': 'bethe-sutherland-pbc', 'Lambda': str(exponent), 'Length': str(length),
                'Particles': '4', 'Precision': 'fp64', 'Calculation': 'label-window spectrum',
                'Label window': f'[-{window},{window}]', 'States enumerated': str(count),
                'Status': 'exact spectral rules',
                'Coverage': 'window only; no global low-energy completeness claimed',
                'Statistics/domain': 'periodic bosons; collision branch |x_i-x_j|^lambda',
                'Momentum': 'P=2*pi*momentum_index/L; not reduced modulo 2*pi'}
    # Exact spectral-rule outputs have no iterative per-row convergence status.
    metadata, rows = read_csv_export(text, columns, expected, row_status=None)
    ground = finite_number(metadata, 'Ground energy')
    if len(rows) != count or ground < 0:
        raise ValueError('Incomplete window or invalid ground energy')
    parsed = []
    for i, row in enumerate(rows):
        if row['state_id'] != str(i):
            raise ValueError('Bad state IDs')
        labels = tuple(map(int, row['labels'].split()))
        if len(labels) != 4 or tuple(sorted(labels)) != labels or any(abs(n) > window for n in labels):
            raise ValueError('Invalid labels')
        numbers = {key: finite_number(row, key) for key in ('energy', 'gap', 'p')}
        momentum = int(row['momentum_index'])
        if sum(labels) != momentum or min(numbers['energy'], numbers['gap']) < 0:
            raise ValueError('Invalid energy/momentum')
        if not math.isclose(numbers['p'], 2*math.pi/length*momentum, rel_tol=0, abs_tol=2e-14):
            raise ValueError('Wrong momentum units')
        if not math.isclose(numbers['energy'], ground+numbers['gap'], rel_tol=0, abs_tol=1e-12):
            raise ValueError('Energy and gap disagree')
        parsed.append(dict(numbers, labels=labels, momentum_index=momentum))
    wanted = set(itertools.combinations_with_replacement(range(-window, window+1), 4))
    if {r['labels'] for r in parsed} != wanted:
        raise ValueError('Duplicate or missing states in accepted window')
    return metadata, parsed


def load_cases(solver=None):
    exports = {}
    for name, case in CASES.items():
        exponent, length, window = case
        if solver:
            command = [str(solver), '4', '--length', str(length), '--lambda', str(exponent),
                       '--levels', 'all', '--window', str(window), '--precision', 'fp64', '--format', 'csv']
            text = subprocess.run(command, check=True, capture_output=True, text=True).stdout
        else:
            text = (DATA / f'{name}.csv').read_text()
        exports[name] = text, read_export(text, *case)
    if solver:
        for name, (text, _) in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return {name: case for name, (_, case) in exports.items()}


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-sutherland-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    q2 = (math.pi/2)**2

    def points(ax, rows, **kwargs):
        ax.scatter([r['momentum_index'] for r in rows], [r['gap']/q2 for r in rows], s=35, **kwargs)
        ax.set(xlabel=r'Momentum index $P L/(2\pi)$', ylabel=r'Gap / $(2\pi/L)^2$')
        ax.grid(alpha=.15)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True, layout='constrained')
    for ax, exponent, color in zip(axes, (0, 1, 2), ('#0072B2', '#D55E00', '#009E73')):
        points(ax, cases[f'sutherland-lambda{exponent}'][1], color=color)
        ax.set(title=rf'$\lambda={exponent}$, $N=L=4$', xticks=[-4, -2, 0, 2, 4], ylim=(-.5, 21))
    save_svg(fig, FIGURES / 'sutherland-collision-branches.svg')
    plt.close(fig)

    small = cases['sutherland-lambda1'][1]
    keys = {r['labels'] for r in small}
    extra = [r for r in cases['sutherland-window2'][1] if r['labels'] not in keys]
    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    points(ax, extra, facecolors='none', edgecolors='#D55E00', label='Additional W=2 states')
    points(ax, small, color='#0072B2', label='W=1 states')
    ax.set(title=r'Expanding the label window: $\lambda=1$, $N=L=4$', ylim=(-1, None))
    ax.legend(fontsize=9)
    save_svg(fig, FIGURES / 'sutherland-windows.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} Sutherland exports.')
