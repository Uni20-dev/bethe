#!/usr/bin/env python3
"""Validate and plot Hubbard tutorial exports, optionally regenerating with --solver.

All plotted energies are exported columns, not a Python dispersion calculation.
Data are replaced only after every requested solver run and validation succeeds.
"""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import solver_executable, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
POINTS = 129
# name: (density, interaction convention, energy reference)
CASES = {
    'hubbard-half-symmetric': (1, 'symmetric', 'hamiltonian'),
    'hubbard-half-unshifted': (1, 'unshifted', 'hamiltonian'),
    'hubbard-half-fermi': (1, 'unshifted', 'fermi'),
    'hubbard-doped-symmetric': (.75, 'symmetric', 'fermi'),
    'hubbard-doped-unshifted': (.75, 'unshifted', 'fermi'),
}


def momentum_ranges(density):
    if density == 1:
        return {'spinon': (0, 1), 'holon': (-1, 1), 'antiholon': (-1, 1)}
    return {'spinon': (0, density), 'holon': (-density / 2, 3 * density / 2),
            'charge-particle': (density / 2, 2 - 3 * density / 2)}


def read_export(text, density, convention, reference):
    error_column = 'energy_quad_error' if density == 1 else 'energy_mesh_error'
    columns = ['branch', 'p', 'p_over_pi', 'energy', 'symmetric_energy', 'delta_n', 'spin',
               'parameter', error_column, 'momentum_error', 'evaluations', 'iterations',
               'status', 'fermi_energy']
    expected = {'Program': 'bethe-hubbard-dispersion', 'U (t=1)': '4',
                'Density N/L': str(density), 'Energy convention': convention,
                'Energy reference': reference, 'Precision': 'fp64', 'Branch selection': 'all',
                'Points per branch': str(POINTS), 'Rows': str(3 * POINTS),
                'Momentum': ('one-site radians; antiholon p = holon p - pi at the same bare k'
                             if density == 1 else
                             'unwrapped one-site radians; hole offset pi*n/2; particle offset -pi*n/2')}
    if density < 1:
        expected['Background status'] = 'converged'
    metadata, rows = read_csv_export(text, columns, expected)
    for key in ('Mu symmetric', 'Mu unshifted'):
        finite_number(metadata, key)
    branches = {branch: [] for branch in momentum_ranges(density)}
    for row in rows:
        branch = row['branch']
        if branch not in branches:
            raise ValueError(f'Unexpected branch: {branch}')
        numbers = {key: finite_number(row, key) for key in (
            'p', 'p_over_pi', 'energy', 'symmetric_energy', 'fermi_energy',
            'delta_n', 'spin', error_column, 'momentum_error', 'iterations')}
        delta_n, spin = (0, .5) if branch == 'spinon' else (-1 if branch == 'holon' else 1, 0)
        if numbers['delta_n'] != delta_n or numbers['spin'] != spin:
            raise ValueError(f'Wrong quantum numbers for {branch}')
        if min(numbers[error_column], numbers['momentum_error'], numbers['iterations']) < 0:
            raise ValueError('Negative error estimate/work count')
        if numbers['fermi_energy'] < -1e-12:
            raise ValueError('Negative Fermi-referenced energy')
        # Infinite rapidities are physical at spinon endpoints, never energies.
        parameter = float(row['parameter'])
        if not math.isfinite(parameter) and not (
                branch == 'spinon' and numbers['p_over_pi'] in (0, density)
                and parameter == (math.inf if numbers['p_over_pi'] == 0 else -math.inf)):
            raise ValueError('Invalid Bethe parameter')
        if density == 1:
            if finite_number(row, 'evaluations') < 0:
                raise ValueError('Negative evaluation count')
        elif row['evaluations'] != '':
            raise ValueError('Doped rows do not supply per-point evaluation counts')
        branches[branch].append(numbers)
    for branch, (start, end) in momentum_ranges(density).items():
        rows = branches[branch]
        if len(rows) != POINTS:
            raise ValueError(f'{branch}: expected {POINTS} points, got {len(rows)}')
        for i, row in enumerate(rows):
            expected_p = start + (end - start) * i / (POINTS - 1)
            if not math.isclose(row['p_over_pi'], expected_p, rel_tol=0, abs_tol=1e-14):
                raise ValueError('Incomplete/out-of-order momentum grid')
            if not math.isclose(row['p'], math.pi * row['p_over_pi'], rel_tol=0, abs_tol=2e-14):
                raise ValueError('Momentum units disagree')
    return metadata, branches


def load_cases(solver=None):
    exports = {}
    for name, case in CASES.items():
        density, convention, reference = case
        if solver:
            command = [str(solver), '--u', '4', '--density', str(density), '--branch', 'all',
                       '--convention', convention, '--reference', reference,
                       '--precision', 'fp64', '--points', str(POINTS), '--format', 'csv']
            text = subprocess.run(command, check=True, text=True, capture_output=True).stdout
        else:
            text = (DATA / f'{name}.csv').read_text()
        exports[name] = (text, read_export(text, *case))
    if solver:
        DATA.mkdir(parents=True, exist_ok=True)
        for name, (text, _) in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return {name: result for name, (_, result) in exports.items()}


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'svg.hashsalt': 'bethe-hubbard-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = {'spinon': '#0072B2', 'holon': '#D55E00',
              'antiholon': '#009E73', 'charge-particle': '#009E73'}
    labels = {'spinon': 'Spinon', 'holon': r'Holon ($\Delta N=-1$)',
              'antiholon': r'Antiholon ($\Delta N=+1$)',
              'charge-particle': r'Charge particle ($\Delta N=+1$)'}

    def line(ax, branches, branch, column='energy'):
        rows = branches[branch]
        ax.plot([r['p_over_pi'] for r in rows], [r[column] for r in rows],
                color=colors[branch], linewidth=2, label=labels[branch],
                linestyle='--' if branch in ('antiholon', 'charge-particle') else '-')

    def finish(fig, name):
        for ax in fig.axes:
            ax.set_xlabel(r'One-site momentum $p/\pi$')
            ax.set_ylabel(r'Excitation energy / $t$')
            ax.axhline(0, color='0.4', linewidth=.8)
            ax.grid(alpha=.15)
        save_svg(fig, FIGURES / name)
        plt.close(fig)

    half = cases['hubbard-half-symmetric'][1]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    line(axes[0], half, 'spinon')
    axes[0].set(title=r'$U/t=4$, $n=1$ · gapless spin', xlim=(0, 1), ylim=(0, 1.4))
    for branch in ('holon', 'antiholon'):
        line(axes[1], half, branch)
    axes[1].set(title='Gapped charge · symmetric Hamiltonian', xlim=(-1, 1), ylim=(0, 5.4))
    axes[1].legend(loc='upper center', fontsize=9, ncols=2)
    finish(fig, 'hubbard-half-lines.svg')

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True, layout='constrained')
    for ax, name, title in zip(axes,
                              ('hubbard-half-symmetric', 'hubbard-half-unshifted', 'hubbard-half-fermi'),
                              ('Symmetric Hamiltonian', 'Unshifted Hamiltonian', 'Unshifted, Fermi reference')):
        for branch in ('holon', 'antiholon'):
            line(ax, cases[name][1], branch)
        ax.set(title=title, xlim=(-1, 1), ylim=(-1.8, 7.4), xticks=[-1, -.5, 0, .5, 1])
    axes[0].legend(loc='upper center', fontsize=9)
    finish(fig, 'hubbard-half-conventions.svg')

    doped = cases['hubbard-doped-symmetric'][1]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    line(axes[0], doped, 'spinon')
    axes[0].set(title=r'$U/t=4$, $n=0.75$ · spinon', xlim=(0, .75), ylim=(0, 1.05),
                xticks=[0, .25, .5, .75])
    for branch in ('holon', 'charge-particle'):
        line(axes[1], doped, branch)
    axes[1].set(title='Charge holes and particles · Fermi reference', xlim=(-.45, 1.2), ylim=(0, 3.5))
    axes[1].legend(loc='upper center', fontsize=9, ncols=2)
    finish(fig, 'hubbard-doped-lines.svg')

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True, layout='constrained')
    for ax, column, title in zip(axes, ('symmetric_energy', 'fermi_energy'),
                                 ('Symmetric Hamiltonian difference', r'Subtract $\mu_{\mathrm{sym}}\Delta N$')):
        for branch in ('holon', 'charge-particle'):
            line(ax, doped, branch, column)
        ax.set(title=title, xlim=(-.45, 1.2), ylim=(-2, 5.4))
    axes[0].legend(loc='upper center', fontsize=9, ncols=2)
    finish(fig, 'hubbard-doped-reference.svg')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable, help='Regenerate data with this executable')
    parser.add_argument('--check', action='store_true', help='Validate without plotting')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} exports ({3 * POINTS} rows each).')


if __name__ == '__main__':
    main()
