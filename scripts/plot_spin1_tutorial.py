#!/usr/bin/env python3
"""Plot actual SU(3)/ULS and TB exports; --bin-dir regenerates the four CSVs."""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
POINTS = 193  # Include thirds, halves, and every folded momentum image.
BRANCHES = {
    'su3': {'3': (4/3, '3'), 'bar3': (2/3, 'bar3'),
            'two-soliton': (2, '1+8'), 'four-soliton': (2, '(1+8)x(1+8)')},
    'tb': {'spinon': (1, 'S=1/2'), 'two-spinon': (2, 'S=0,1'),
           'four-spinon': (2, 'S=0,1,2')},
}


def read_export(text, model, folded):
    labels = ['representations'] if model == 'su3' else ['spin', 'sectors']
    columns = ['branch', *labels, 'p', 'p_over_pi', 'cell_momentum',
               'energy', 'lower', 'upper', 'status']
    expected = {'Program': f'bethe-{model}-dispersion', 'Exchange J': '1',
                'Precision': 'fp64', 'Branches': 'all', 'Points per branch': str(POINTS),
                'Continuum momentum': ('three-site folded envelope' if model == 'su3'
                                       else 'two-site folded envelope') if folded else 'one-site unfolded'}
    _, rows = read_csv_export(text, columns, expected)
    result = {b: [] for b in BRANCHES[model]}
    for row in rows:
        branch = row['branch']
        if branch not in result:
            raise ValueError(f'Unexpected branch {branch}')
        elementary = branch in ('3', 'bar3', 'spinon')
        if row[labels[-1]] != BRANCHES[model][branch][1]:
            raise ValueError('Incorrect representation/spin sector')
        if model == 'tb' and row['spin'] != ('0.5' if elementary else ''):
            raise ValueError('Incorrect elementary spin label')
        required = ('energy',) if elementary else ('lower', 'upper')
        absent = ('lower', 'upper') if elementary else ('energy',)
        if any(row[key] != '' for key in absent):
            raise ValueError('Non-applicable energy columns must be empty')
        numbers = {key: finite_number(row, key) for key in
                   ('p', 'p_over_pi', 'cell_momentum', *required)}
        if any(numbers[key] < 0 for key in required):
            raise ValueError('Negative excitation energy')
        if not elementary and numbers['lower'] > numbers['upper']:
            raise ValueError('Reversed continuum bounds')
        result[branch].append(numbers)
    cell = 3 if model == 'su3' else 2
    for branch, rows in result.items():
        if len(rows) != POINTS:
            raise ValueError('Incomplete branch')
        end = BRANCHES[model][branch][0]
        for i, row in enumerate(rows):
            if not math.isclose(row['p_over_pi'], end*i/(POINTS-1), rel_tol=0, abs_tol=1e-14):
                raise ValueError('Incorrect momentum grid')
            if not math.isclose(row['p'], math.pi*row['p_over_pi'], rel_tol=0, abs_tol=1e-14):
                raise ValueError('Incorrect momentum units')
            phase_error = math.remainder(row['cell_momentum'] - cell*row['p'], 2*math.pi)
            if abs(phase_error) > 1e-13:
                raise ValueError('Incorrect unit-cell momentum')
    return result


def load_cases(bin_dir=None):
    exports = {}
    for model in BRANCHES:
        for folded in (False, True):
            name = model + ('-folded' if folded else '-unfolded')
            if bin_dir:
                command = [str(bin_dir / f'bethe-{model}-dispersion'), '--exchange', '1',
                           '--points', str(POINTS), '--branch', 'all', '--precision', 'fp64', '--format', 'csv']
                if folded:
                    command.append('--folded')
                text = subprocess.run(command, check=True, capture_output=True, text=True).stdout
            else:
                text = (DATA / f'{name}.csv').read_text()
            exports[name] = text, read_export(text, model, folded)
    if bin_dir:
        for name, (text, _) in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return {name: rows for name, (_, rows) in exports.items()}


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-spin1-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})

    def bands(ax, rows, model):
        suffix = 'soliton' if model == 'su3' else 'spinon'
        for count, color, alpha in [('four', '#009E73', .12), ('two', '#0072B2', .25)]:
            branch = rows[f'{count}-{suffix}']
            x = [r['p_over_pi'] for r in branch]
            low, high = ([r[k] for r in branch] for k in ('lower', 'upper'))
            ax.fill_between(x, low, high, color=color, alpha=alpha)
            ax.plot(x, low, color=color, linewidth=1.5)
            ax.plot(x, high, color=color, linewidth=1.5, label=f'{count.title()} {suffix}s')
        ax.set(xlim=(0, 2), ylim=(0, 11 if model == 'su3' else 29),
               xlabel=r'Total one-site momentum $Q/\pi$')
        ax.legend(loc='upper center', fontsize=9, ncols=2)

    def save(fig, name):
        for ax in fig.axes:
            ax.set_ylabel(r'Excitation energy / $J$')
            ax.grid(alpha=.15)
        save_svg(fig, FIGURES / name)
        plt.close(fig)

    for model in BRANCHES:
        rows = cases[f'{model}-unfolded']
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
        particles = ('3', 'bar3') if model == 'su3' else ('spinon',)
        for branch, color, label in zip(particles, ('#D55E00', '#CC79A7'),
                                        (r'$\mathbf{3}$', r'$\overline{\mathbf{3}}$') if model == 'su3'
                                        else ('Spin 1/2 spinon',)):
            axes[0].plot([r['p_over_pi'] for r in rows[branch]], [r['energy'] for r in rows[branch]],
                         color=color, label=label, linewidth=2)
        axes[0].set(title='ULS elementary branches' if model == 'su3' else 'TB elementary spinon',
                    xlabel=r'One-site momentum $p/\pi$', ylim=(0, 4.4 if model == 'su3' else 7.5))
        axes[0].legend(loc='upper center', ncols=2)
        bands(axes[1], rows, model)
        axes[1].set_title('Specified-particle continua · unfolded')
        save(fig, f'{model}-spectrum.svg')

        fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
        for ax, mode in zip(axes, ('unfolded', 'folded')):
            bands(ax, cases[f'{model}-{mode}'], model)
            ax.set_title('One-site momentum' if mode == 'unfolded' else
                         ('Three-site envelope' if model == 'su3' else 'Two-site envelope'))
        save(fig, f'{model}-folding.svg')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, help='Directory containing both built frontends')
    parser.add_argument('--check', action='store_true', help='Validate without plotting dependencies')
    args = parser.parse_args()
    cases = load_cases(args.bin_dir.resolve() if args.bin_dir else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} spin-1 tutorial exports.')
