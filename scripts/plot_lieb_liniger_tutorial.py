#!/usr/bin/env python3
"""Lieb-Liniger examples from frontend exports; --solver regenerates all data."""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
POINTS = 129
CASES = {'ll-c1-n1': (1, 1), 'll-c10-n1': (10, 1), 'll-c100-n1': (100, 1),
         'll-c2-n2': (2, 2)}


def read_export(text, coupling, density):
    columns = ['branch', 'p', 'p_over_pi_n', 'energy', 'rapidity', 'edge_distance',
               'energy_error', 'momentum_error', 'iterations', 'status']
    expected = {'Program': 'bethe-lieb-liniger-dispersion', 'c': str(coupling),
                'Density N/L': str(density), 'Precision': 'fp64', 'Branches': 'all',
                'Background status': 'converged', 'Points per branch': str(POINTS),
                'Energy reference': 'fixed-N excitation gap above the bulk ground state',
                'Momentum convention': 'physical inverse length; not reduced modulo 2*pi'}
    metadata, rows = read_csv_export(text, columns, expected)
    for key in ('Ground energy/length', 'Chemical potential', 'Fermi rapidity Q'):
        if finite_number(metadata, key) <= 0:
            raise ValueError(f'Invalid background {key}')
    branches = {'type-i': [], 'type-ii': []}
    for row in rows:
        if row['branch'] not in branches:
            raise ValueError('Unknown branch')
        numbers = {key: finite_number(row, key) for key in columns[1:-1]}
        if min(numbers[key] for key in ('energy', 'edge_distance', 'energy_error', 'momentum_error', 'iterations')) < 0:
            raise ValueError('Negative gap/error/work')
        branches[row['branch']].append(numbers)
    for rows in branches.values():
        if len(rows) != POINTS:
            raise ValueError('Incomplete grid')
        for i, row in enumerate(rows):
            if not math.isclose(row['p_over_pi_n'], 2*i/(POINTS-1), rel_tol=0, abs_tol=1e-14):
                raise ValueError('Wrong momentum grid')
            if not math.isclose(row['p'], math.pi*density*row['p_over_pi_n'], rel_tol=0, abs_tol=2e-14):
                raise ValueError('Wrong momentum units')
    return metadata, branches


def load_cases(solver=None):
    exports = {}
    for name, (coupling, density) in CASES.items():
        if solver:
            command = [str(solver), '--c', str(coupling), '--density', str(density),
                       '--points', str(POINTS), '--precision', 'fp64', '--format', 'csv']
            text = subprocess.run(command, check=True, capture_output=True, text=True).stdout
        else:
            text = (DATA / f'{name}.csv').read_text()
        exports[name] = text, read_export(text, coupling, density)
    if solver:
        for name, (text, _) in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return {name: case for name, (_, case) in exports.items()}


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-ll-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    for branch, ax, title in zip(('type-i', 'type-ii'), axes, ('Type I · particle-like', 'Type II · hole-like')):
        for name, color in zip(('ll-c1-n1', 'll-c10-n1', 'll-c100-n1'), ('#0072B2', '#D55E00', '#009E73')):
            _, rows = cases[name]
            ax.plot([r['p_over_pi_n'] for r in rows[branch]], [r['energy'] for r in rows[branch]],
                    color=color, label=rf'$\gamma={CASES[name][0]}$', linewidth=2)
        # Explicitly labelled analytic limiting curve, not replacement data.
        grid = [r['p_over_pi_n'] for r in rows[branch]]
        ax.plot(grid, [math.pi**2*x*(2+x if branch == 'type-i' else 2-x) for x in grid],
                ':', color='0.3', label='Tonks limit', linewidth=1.5)
        ax.set(title=title, xlabel=r'$p/(\pi n)$', ylabel=r'$E/n^2$', xlim=(0, 2), ylim=(0, None))
        ax.grid(alpha=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=4, fontsize=9)
    save_svg(fig, FIGURES / 'll-branches.svg')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    for name, color in [('ll-c1-n1', '#0072B2'), ('ll-c2-n2', '#D55E00')]:
        n = CASES[name][1]
        for branch, style in [('type-i', '-'), ('type-ii', '--')]:
            rows = cases[name][1][branch]
            for ax, scaled in zip(axes, (False, True)):
                ax.plot([r['p']/(math.pi*n if scaled else math.pi) for r in rows],
                        [r['energy']/(n*n if scaled else 1) for r in rows],
                        linestyle=style, color=color, linewidth=2 if n == 1 else 1.2,
                        label=f'n={n}, {branch}')
    axes[0].set(title=r'Physical units, fixed $\gamma=1$', xlabel=r'$p/\pi$', ylabel='Excitation energy')
    axes[1].set(title='Density-scaled curves coincide', xlabel=r'$p/(\pi n)$', ylabel=r'$E/n^2$')
    for ax in axes:
        ax.set_ylim(0, None)
        ax.grid(alpha=.15)
        ax.legend(fontsize=8)
    save_svg(fig, FIGURES / 'll-density-scaling.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} Lieb-Liniger exports.')
