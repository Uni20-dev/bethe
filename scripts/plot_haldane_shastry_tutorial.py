#!/usr/bin/env python3
"""Complete small-ring motif spectra and their SU(2) content, from exports."""

import argparse
from collections import Counter
import itertools
import math
from pathlib import Path

from tutorial_common import solver_executable, capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
SIZES = (8, 9)


def read_exports(texts, n):
    expected = {'Program': 'bethe-haldane-shastry-pbc', 'Sites': str(n), 'Precision': 'fp64',
                'Calculation': 'energy-ordered motif spectrum', 'Status': 'exact spectral rules',
                'Degeneracy': 'whole Yangian multiplet; S_max is not a unique total spin',
                'Spin content': 'whole Yangian multiplet; multiplicity counts SU(2) irreps',
                'Momentum': 'P=2*pi*momentum_index/N modulo 2*pi'}
    columns = ['state_id', 'motif', 'energy', 'gap', 'momentum_index', 'p', 'spinons', 's_max', 'degeneracy']
    meta, raw = read_csv_export(texts['levels'], columns, expected, row_status=None)
    wanted = {tuple(i+1 for i, bit in enumerate(bits) if bit)
              for bits in itertools.product((0, 1), repeat=n-1)
              if not any(a and b for a, b in zip(bits, bits[1:]))}
    if len(raw) != len(wanted) or int(meta['Motifs enumerated']) != len(wanted):
        raise ValueError('Incomplete motif spectrum')
    levels = []
    for i, row in enumerate(raw):
        motif = () if row['motif'] == 'empty' else tuple(map(int, row['motif'].split()))
        numbers = {key: finite_number(row, key) for key in ('energy', 'gap', 'p', 's_max')}
        index, spinons, dimension = (int(row[key]) for key in ('momentum_index', 'spinons', 'degeneracy'))
        if (row['state_id'] != str(i) or motif not in wanted or index != sum(motif) % n
                or spinons != n-2*len(motif) or numbers['s_max'] != spinons/2 or dimension <= 0
                or numbers['gap'] < 0 or not math.isclose(numbers['p'], 2*math.pi*index/n, abs_tol=1e-14)):
            raise ValueError('Invalid motif or quantum numbers')
        levels.append(dict(numbers, motif=motif, momentum_index=index, spinons=spinons, degeneracy=dimension))
    if {r['motif'] for r in levels} != wanted:
        raise ValueError('Missing/duplicate motif')
    ground = min(r['energy'] for r in levels)
    for row in levels:
        if not math.isclose(row['gap'], row['energy']-ground, abs_tol=2e-13):
            raise ValueError('Wrong gap reference')
    columns = ['state_id', 'spin', 'multiplicity', 'updates', 'complete', 'status']
    spin_meta, raw = read_csv_export(texts['spin_content'], columns, expected, row_status=None)
    for key in ('Command', 'Bethe revision', 'Uni20 revision', 'Motifs enumerated'):
        if spin_meta[key] != meta[key]:
            raise ValueError('Tables from different calculations')
    content = {i: Counter() for i in range(len(levels))}
    for row in raw:
        state_id, multiplicity = int(row['state_id']), int(row['multiplicity'])
        spin = finite_number(row, 'spin')
        if (state_id not in content or row['complete'] != 'true' or row['status'] != 'complete'
                or multiplicity <= 0 or int(row['updates']) < 0 or spin < 0 or 2*spin != int(2*spin)
                or int(2*spin) % 2 != n % 2 or spin > levels[state_id]['s_max']
                or spin in content[state_id]):
            raise ValueError('Invalid/incomplete spin decomposition')
        content[state_id][spin] = multiplicity
    for i, terms in content.items():
        if sum(int(2*s+1)*m for s, m in terms.items()) != levels[i]['degeneracy']:
            raise ValueError('Spin content does not reproduce motif dimension')
    if sum(r['degeneracy'] for r in levels) != 2**n:
        raise ValueError('Hilbert-space count mismatch')
    return levels, content


def load_cases(solver=None):
    exports, cases = {}, {}
    for n in SIZES:
        texts = (capture_csv_tables([solver, n, '--levels', 'all'], ('levels', 'spin_content')) if solver else
                 {name: (DATA / f'hs-n{n}-{name}.csv').read_text() for name in ('levels', 'spin_content')})
        cases[n] = read_exports(texts, n)
        exports[n] = texts
    if solver:
        for n, texts in exports.items():
            for name, text in texts.items():
                (DATA / f'hs-n{n}-{name}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-hs-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    colors = ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#777777')
    for ax, n in zip(axes, SIZES):
        levels, _ = cases[n]
        for spinons, color in zip(range(n % 2, n+1, 2), colors):
            rows = [r for r in levels if r['spinons'] == spinons]
            ax.scatter([r['p']/math.pi for r in rows], [r['gap'] for r in rows],
                       s=35, color=color, label=f'{spinons} spinons')
        ax.set(title=f'Complete motif spectrum, N={n}', xlabel=r'$P/\pi$', ylabel=r'$E-E_0$', xlim=(-.05, 2.05))
        ax.grid(alpha=.15)
        ax.legend(fontsize=8)
    save_svg(fig, FIGURES / 'hs-motif-spectra.svg')
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout='constrained')
    for ax, n in zip(axes, SIZES):
        _, content = cases[n]
        totals = sum(content.values(), Counter())
        spins = sorted(totals)
        ax.bar(spins, [totals[s] for s in spins], width=.65, color='#0072B2')
        for s in spins:
            ax.text(s, totals[s]+.5, str(totals[s]), ha='center', fontsize=9)
        ax.set(title=f'All SU(2) irreps, N={n}', xlabel='Total spin S', ylabel='Number of multiplets',
               xticks=spins, ylim=(0, max(totals.values())*1.18))
        ax.grid(axis='y', alpha=.15)
    save_svg(fig, FIGURES / 'hs-spin-counts.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} complete Haldane-Shastry spectra and spin tables.')
