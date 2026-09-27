#!/usr/bin/env python3
"""XYZ spinons and parity-labelled bound lines from native frontend exports."""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
POINTS = 129
CASES = {'xyz-eta04': (.4, 1), 'xyz-xy': (.5, 1), 'xyz-eta075': (.75, 1), 'xyz-exchange2': (.75, 2)}


def read_export(text, eta, exchange):
    indices = [s for s in range(1, 4) if s*(1-eta) < eta]
    expected = {'Program': 'bethe-xyz-dispersion', 'Theta t': '1', 'Exchange J': str(exchange),
                'Precision': 'fp64', 'Branches': 'all', 'Spinon status': 'converged',
                'Points per branch/copy': str(POINTS), 'Selected bound branches': str(len(indices)),
                'Energy reference': 'excitation energy above the infinite-chain ground state',
                'Two-spinon sectors': 'envelope over pi-shifted translation copies; no symmetry resolution',
                'Bound sectors': 'delta_rx=(-1)^s; delta_rz=+1,-1 relative to one vacuum',
                'Bound momentum': 'reduced_q=p-pi*(delta_rz==-1) modulo 2*pi',
                'Two-site translation phase': '2*p modulo 2*pi; radians/cell'}
    columns = ['branch', 's', 'delta_rx', 'delta_rz', 'p', 'p_over_pi', 'cell_momentum', 'reduced_q',
               'energy', 'lower', 'upper', 'status']
    meta, raw = read_csv_export(text, columns, expected)
    if finite_number(meta, 'Eta') != eta:
        raise ValueError('Wrong eta')
    for key in ('Jx', 'Jy', 'Jz', 'Single-spinon gap', 'Spinon maximum energy'):
        finite_number(meta, key)
    branches = {('spinon', 0, 0): [], ('two-spinon', 0, 0): []}
    branches.update({('bound', s, rz): [] for s in indices for rz in (1, -1)})
    for row in raw:
        branch = row['branch']
        s, rz = (int(row['s']), int(row['delta_rz'])) if branch == 'bound' else (0, 0)
        key = (branch, s, rz)
        if key not in branches:
            raise ValueError('Unexpected branch')
        used = ['p', 'p_over_pi', 'cell_momentum']
        used += ['lower', 'upper'] if branch == 'two-spinon' else ['energy']
        if branch == 'bound':
            used.append('reduced_q')
            if int(row['delta_rx']) != (-1)**s:
                raise ValueError('Wrong bound parity')
        elif any(row[k] != '' for k in ('s', 'delta_rx', 'delta_rz')):
            raise ValueError('Unexpected ring parity')
        if any(row[k] != '' for k in ('reduced_q', 'energy', 'lower', 'upper') if k not in used):
            raise ValueError('Unexpected energy/coordinate column')
        numbers = {k: finite_number(row, k) for k in used}
        if min(numbers.values()) < 0 or (branch == 'two-spinon' and numbers['lower'] > numbers['upper']):
            raise ValueError('Invalid energy/range')
        branches[key].append(numbers)
    for (branch, s, rz), rows in branches.items():
        if len(rows) != POINTS:
            raise ValueError('Incomplete branch/copy')
        end = math.pi if branch == 'spinon' else 2*math.pi
        for i, row in enumerate(rows):
            p = end*i/(POINTS-1)
            cell = 2*p % (2*math.pi)
            for k, wanted in [('p', p), ('p_over_pi', p/math.pi), ('cell_momentum', cell)]:
                if not math.isclose(row[k], wanted, rel_tol=0, abs_tol=2e-14):
                    raise ValueError(f'Wrong {k} grid/convention')
            if s:
                q = p if rz == 1 else p+math.pi if p < math.pi else p-math.pi
                if not math.isclose(row['reduced_q'], q, rel_tol=0, abs_tol=2e-14):
                    raise ValueError('Wrong reduced bound momentum')
    return meta, branches


def load_cases(solver=None):
    exports, cases = {}, {}
    for name, (eta, exchange) in CASES.items():
        text = (subprocess.run([str(solver), '--eta', str(eta), '--t', '1', '--exchange', str(exchange),
                                '--points', str(POINTS), '--precision', 'fp64', '--format', 'csv'],
                               check=True, capture_output=True, text=True).stdout
                if solver else (DATA / f'{name}.csv').read_text())
        cases[name] = read_export(text, eta, exchange)
        exports[name] = text
    if solver:
        for name, text in exports.items():
            (DATA / f'{name}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-xyz-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(7, 3.8), layout='constrained')
    for name, color in zip(('xyz-eta04', 'xyz-xy', 'xyz-eta075'), ('#0072B2', '#009E73', '#D55E00')):
        rows = cases[name][1]['spinon', 0, 0]
        ax.plot([r['p_over_pi'] for r in rows], [r['energy'] for r in rows], color=color,
                label=rf'$\eta={CASES[name][0]}$', linewidth=2)
    ax.set(title='XYZ single-spinon bands, t=1 and J=1', xlabel=r'$p/\pi$', ylabel='Excitation energy',
           xlim=(0, 1), ylim=(0, None))
    ax.legend(); ax.grid(alpha=.15)
    save_svg(fig, FIGURES / 'xyz-spinon-bands.svg'); plt.close(fig)

    rows = cases['xyz-eta075'][1]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout='constrained')
    continuum = rows['two-spinon', 0, 0]
    axes[0].fill_between([r['p_over_pi'] for r in continuum], [r['lower'] for r in continuum],
                         [r['upper'] for r in continuum], color='0.85', label='All-sector two-spinon envelope')
    for s, color in [(1, '#D55E00'), (2, '#0072B2')]:
        for rz, style in [(1, '-'), (-1, '--')]:
            branch = rows['bound', s, rz]
            label = rf'$s={s},\ \delta r_z={rz:+d}$'
            axes[0].plot([r['p_over_pi'] for r in branch], [r['energy'] for r in branch],
                         style, color=color, linewidth=1.6, label=label)
            ordered = sorted(branch, key=lambda r: r['reduced_q'])
            axes[1].plot([r['reduced_q']/math.pi for r in ordered], [r['energy'] for r in ordered],
                         style, color=color, linewidth=2 if rz == 1 else 1,
                         marker='o' if rz == -1 else None, markevery=16, markersize=3)
    axes[0].set(title='Physical momentum: two parity copies', xlabel=r'$p/\pi$')
    axes[1].set(title='Reduced momentum: copies coincide', xlabel=r'$Q/\pi$')
    for ax in axes:
        ax.set(ylabel='Excitation energy', xlim=(0, 2), ylim=(0, 2.05))
        ax.grid(alpha=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=3, fontsize=8)
    fig.suptitle(r'XYZ bound branches: $\eta=0.75$, t=1, J=1', fontsize=12)
    save_svg(fig, FIGURES / 'xyz-bound-copies.svg'); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} XYZ exports.')
