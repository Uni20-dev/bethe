#!/usr/bin/env python3
"""Deform a small bosonic ring from free hopping to the phase model."""

import argparse
import itertools
import math
from pathlib import Path

from tutorial_common import solver_executable, capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
CASES = {'q-boson-free': '0', 'q-boson-eta05': '0.5', 'q-boson-phase': None}
TABLES = ('states', 'reference', 'roots')
L, N = 5, 3
COUNT = math.comb(L+N-1, N)


def read_exports(texts, eta):
    expected = {'Program': 'bethe-q-boson-pbc', 'Sites L': str(L), 'Particles N': str(N), 'Precision': 'fp64',
                'Status': 'converged canonical family', 'Calculation': 'canonical fixed-N excitation scan',
                'Deformation': 'phase limit' if eta is None else 'free bosons' if eta == '0' else 'finite eta',
                'Candidate count': str(COUNT), 'Converged count': str(COUNT), 'Retained count': str(COUNT),
                'Energy reference': '+2N diagonal shift included; subtract 2N for pure hopping',
                'Momentum convention': '(-pi,pi]; index=sum(m) modulo L; lifted roots'}
    if eta is not None:
        expected['eta=log(q)'] = eta
    columns = ['state_id', 'energy', 'energy_per_site', 'gap', 'momentum_index', 'momentum',
               'residual', 'iterations', 'converged', 'status']
    meta, raw = read_csv_export(texts['states'], columns, expected)
    ref_meta, reference = read_csv_export(texts['reference'], columns, expected)
    root_columns = ['state_id', 'root_index', 'I', 'm', 'k', 'deviation']
    root_meta, raw_roots = read_csv_export(texts['roots'], root_columns, expected, row_status=None)
    for other in (ref_meta, root_meta):
        for key in ('Command', 'Bethe revision', 'Uni20 revision', 'Ground energy'):
            if other[key] != meta[key]:
                raise ValueError('Tables from different calculations')
    ground = finite_number(meta, 'Ground energy')
    tolerance = finite_number(meta, 'Residual tolerance')
    if len(raw) != COUNT or len(reference) != 1 or tolerance <= 0:
        raise ValueError('Incomplete scan/reference')
    states = []
    for i, row in enumerate(raw+reference):
        numbers = {key: finite_number(row, key) for key in ('energy', 'energy_per_site', 'momentum', 'residual', 'iterations')}
        index = int(row['momentum_index'])
        momentum = 2*math.pi/L*(index if index <= L/2 else index-L)
        if (row['state_id'] != str(i) or row['converged'] != 'true' or not 0 <= index < L
                or not 0 <= numbers['residual'] <= tolerance or numbers['iterations'] < 0
                or not math.isclose(numbers['momentum'], momentum, abs_tol=1e-14)
                or not math.isclose(numbers['energy_per_site'], numbers['energy']/L, abs_tol=1e-14)):
            raise ValueError('Invalid state')
        gap = None if i == COUNT else finite_number(row, 'gap')
        if i == COUNT:
            if row['gap'] != '' or not math.isclose(numbers['energy'], ground, abs_tol=1e-13):
                raise ValueError('Bad reference row')
        elif not math.isclose(gap, numbers['energy']-ground, abs_tol=1e-13):
            raise ValueError('Wrong gap reference')
        states.append(dict(numbers, gap=gap, momentum_index=index))
    roots = {i: [] for i in range(COUNT+1)}
    for row in raw_roots:
        i, j, mode = (int(row[k]) for k in ('state_id', 'root_index', 'm'))
        numbers = {key: finite_number(row, key) for key in ('I', 'k', 'deviation')}
        if (i not in roots or j != len(roots[i]) or j >= N or not 0 <= mode < L
                or numbers['I'] != mode+j-(N-1)/2
                or not math.isclose(numbers['k'], 2*math.pi*mode/L+numbers['deviation'], abs_tol=1e-13)):
            raise ValueError('Bad root coordinates or table join')
        roots[i].append(dict(numbers, m=mode))
    patterns = set()
    for i, row in enumerate(states):
        rr = roots[i]
        modes = tuple(r['m'] for r in rr)
        if len(rr) != N or modes != tuple(sorted(modes)) or sum(modes) % L != row['momentum_index']:
            raise ValueError('Missing/unordered roots or momentum mismatch')
        if not math.isclose(row['energy'], sum(4*math.sin(r['k']/2)**2 for r in rr), abs_tol=2e-13):
            raise ValueError('Energy/root mismatch')
        if i < COUNT:
            patterns.add(modes)
        elif modes != (0,)*N:
            raise ValueError('Incorrect reference labels')
    if patterns != set(itertools.combinations_with_replacement(range(L), N)):
        raise ValueError('Missing/duplicate canonical modes')
    return meta, states[:COUNT], roots


def load_cases(solver=None):
    exports, cases = {}, {}
    for name, eta in CASES.items():
        command = [solver, L, '--particles', N, '--excitations', 'all']
        command += ['--phase'] if eta is None else ['--eta', eta]
        texts = (capture_csv_tables(command, TABLES) if solver else
                 {table: (DATA / f'{name}-{table}.csv').read_text() for table in TABLES})
        cases[name] = read_exports(texts, eta)
        exports[name] = texts
    if solver:
        for name, texts in exports.items():
            for table, text in texts.items():
                (DATA / f'{name}-{table}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-q-boson-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    titles = ('Free bosons', r'$\eta=0.5$', 'Phase limit')
    colors = ('#0072B2', '#D55E00', '#009E73')
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True, layout='constrained')
    for ax, (name, _), title, color in zip(axes, CASES.items(), titles, colors):
        rows = cases[name][1]
        ax.scatter([r['momentum']/math.pi for r in rows], [r['energy'] for r in rows], s=35, color=color)
        ax.axhline(float(cases[name][0]['Ground energy']), color=color, linestyle=':', linewidth=1)
        ax.set(title=title, xlabel=r'$P/\pi$', xticks=[-.8, -.4, 0, .4, .8])
        ax.grid(alpha=.15)
    axes[0].set_ylabel(r'Energy (including $+2N$)')
    fig.suptitle('q-boson ring: L=5, N=3, 35 canonical states per panel', fontsize=12)
    save_svg(fig, FIGURES / 'q-boson-spectra.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout='constrained')
    for (name, _), title, color, marker in zip(CASES.items(), titles, colors, ('o', 's', '^')):
        ground_roots = cases[name][2][COUNT]
        axes[0].plot(range(N), [r['k']/math.pi for r in ground_roots], marker+'-', color=color, label=title)
    axes[0].set(title='Ground-state rapidities spread', xlabel='Root index j', ylabel=r'$k_j/\pi$', xticks=range(N))
    axes[0].legend(fontsize=9); axes[0].grid(alpha=.15)
    energies = [float(cases[name][0]['Ground energy']) for name in CASES]
    axes[1].bar(range(3), energies, color=colors)
    for i, energy in enumerate(energies):
        axes[1].text(i, energy+.025, f'{energy:.6f}', ha='center', fontsize=9)
    axes[1].set(title='Shifted ground energy', ylabel=r'$E_0$', xticks=range(3), xticklabels=titles,
                ylim=(0, max(energies)*1.2))
    save_svg(fig, FIGURES / 'q-boson-ground.svg'); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} q-boson scans and auxiliary tables.')
