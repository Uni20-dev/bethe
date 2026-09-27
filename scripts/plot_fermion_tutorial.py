#!/usr/bin/env python3
"""Finite Gaudin-Yang and supersymmetric t-J examples from checked native roots."""

import argparse
import cmath
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
COUPLINGS = (0, .001, .01, .1, 1, 10, 100, 1000)
FILLINGS = (0, 2, 6, 10, 14, 16)
# model, circumference/sites, N_up, N_down, c (unused for t-J)
CASES = {f'gy-c{c}': ('gy', 6, 3, 3, c) for c in COUPLINGS}
CASES.update({'gy-polarized': ('gy', 6, 6, 0, 1), 'gy-scaled': ('gy', 12, 3, 3, .5),
              'gy-imbalanced': ('gy', 8, 5, 3, 1), 'gy-reversed': ('gy', 8, 3, 5, 1),
              'gy-two-body': ('gy', 1, 1, 1, math.pi), 'gy-free-shell': ('gy', 1, 2, 2, 0)})
for n in FILLINGS:
    CASES[f'tj-balanced-n{n}'] = ('tj', 16, n//2, n//2, None)
    CASES[f'tj-polarized-n{n}'] = ('tj', 16, n, 0, None)
CASES.update({'tj-four-noholes': ('tj', 4, 2, 2, None), 'tj-six-noholes': ('tj', 6, 3, 3, None),
              'tj-imbalanced': ('tj', 5, 3, 1, None), 'tj-reversed': ('tj', 5, 1, 3, None)})
STATE = 'state_id sz energy momentum_index p'.split()
END = 'iterations converged status'.split()
ROOT_COLUMNS = 'state_id index quantum_number'.split()


def schemas(case):
    model, length, up, down, c = case
    free = min(up, down) == 0 or (model == 'gy' and c == 0)
    result = {'states': STATE + (['charge_residual', 'spin_residual', 'residual', 'root_c'] if model == 'gy' else ['residual']) + END}
    if model == 'gy':
        result.update({'free_up': ['state_id', 'mode', 'k'], 'free_down': ['state_id', 'mode', 'k']} if free else
                      {'charge_roots': ROOT_COLUMNS+['k'], 'spin_roots': ROOT_COLUMNS+['rapidity']})
    else:
        result.update({'free_modes': ['state_id', 'index', 'mode']} if free else
                      {name: ROOT_COLUMNS+['rapidity'] for name in ('first_roots', 'second_roots')})
    return free, result


def close(actual, expected, tol=3e-11):
    if abs(actual-expected) > tol*max(1, abs(expected)):
        raise ValueError(f'Inconsistent observable/equation: {actual} != {expected}')


def phase(x, width):
    return complex(x, width/2)/complex(x, -width/2)


def read_exports(texts, case):
    model, length, up, down, c = case
    free, columns = schemas(case)
    n, minority = up+down, min(up, down)
    if model == 'gy':
        branch = 'exact free-fermion ground state' if free else 'odd-population sector ground state'
    else:
        branch = 'exact polarized free fermions' if free else 'no-hole XXX reduction' if n == length else 'Sutherland real-root sector'
    expected = {'Program': 'bethe-gaudin-yang-pbc' if model == 'gy' else 'bethe-tj-pbc',
                'Calculation': branch, 'Precision': 'fp64', 'Particles': str(n), 'N_up': str(up), 'N_down': str(down)}
    if model == 'gy': expected['Units'] = 'hbar^2/(2m)=1; interaction 2c delta'
    else: expected['Hamiltonian'] = 't=1, J=2; projected hopping + 2*(S.S-nn/4)'
    if model == 'gy' and not free: expected['Reference spin'] = 'down (spin reversed)' if up < down else 'up'
    tables, meta = {}, None
    for table, schema in columns.items():
        current, raw = read_csv_export(texts[table], schema, expected, row_status='status' if table == 'states' else None)
        if meta is not None and current != meta: raise ValueError('Mismatched table provenance')
        meta = current
        rows = []
        for row in raw:
            if row['state_id'] != '0': raise ValueError('Wrong state join')
            numbers = {key: finite_number(row, key) for key in schema if key not in ('converged', 'status')}
            if 'converged' in row and row['converged'] != 'true': raise ValueError('Unconverged state')
            rows.append(numbers)
        tables[table] = rows
    close(finite_number(meta, 'Length' if model == 'gy' else 'Sites'), length)
    if len(tables['states']) != 1: raise ValueError('Missing/duplicated state')
    state = tables['states'][0]
    close(state['sz'], (up-down)/2)
    close(state['energy'], finite_number(meta, 'Total energy'))
    close(state['p'], 2*math.pi*state['momentum_index']/length)
    if not state['momentum_index'].is_integer() or (model == 'tj' and not 0 <= state['momentum_index'] < length):
        raise ValueError('Invalid momentum index/convention')
    tol = finite_number(meta, 'Residual tolerance')
    if tol <= 0 or not 0 <= state['residual'] <= tol or state['iterations'] < 0:
        raise ValueError('Invalid convergence diagnostics')
    if model == 'gy':
        close(finite_number(meta, 'c'), c)
        close(state['root_c'], c)
        if not all(0 <= state[key] <= tol for key in ('charge_residual', 'spin_residual')):
            raise ValueError('Unresolved nested equation')
        counts = {'free_up': up, 'free_down': down} if free else {'charge_roots': n, 'spin_roots': minority}
    else:
        holes = length-n
        if meta['Holes'] != str(holes): raise ValueError('Wrong holes')
        counts = {'free_modes': n} if free else {'first_roots': holes+minority, 'second_roots': holes}
        for label, table in [('First-level roots', 'first_roots'), ('Second-level roots', 'second_roots')]:
            if int(meta[label]) != counts.get(table, 0): raise ValueError('Wrong nested root count')
    for table, count in counts.items():
        if len(tables[table]) != count: raise ValueError('Incomplete roots/occupations')
        for i, row in enumerate(tables[table]):
            if 'index' in row and row['index'] != i: raise ValueError('Wrong root index')
            if 'quantum_number' in row: close(row['quantum_number'], i-(count-1)/2)
    if free:
        energy, momentum_index = 0, 0
        for table, count in counts.items():
            for i, row in enumerate(tables[table]):
                mode = i-(count-1)//2
                if row['mode'] != mode: raise ValueError('Wrong free periodic shell')
                k = 2*math.pi*mode/length
                if model == 'gy': close(row['k'], k)
                energy += k*k if model == 'gy' else -2*math.cos(k)
                momentum_index += mode
        close(state['energy'], energy)
        close(state['momentum_index'], momentum_index if model == 'gy' else momentum_index % length)
    elif model == 'gy':
        charge = [r['k'] for r in tables['charge_roots']]
        spin = [r['rapidity'] for r in tables['spin_roots']]
        close(state['energy'], sum(k*k for k in charge))
        close(state['p'], sum(charge))
        for k in charge: close(cmath.exp(1j*k*length), math.prod(phase(k-a, c) for a in spin))
        for i, a in enumerate(spin):
            close(math.prod(phase(a-k, c) for k in charge), math.prod(phase(a-b, 2*c) for j, b in enumerate(spin) if i != j))
    else:
        first = [r['rapidity'] for r in tables['first_roots']]
        second = [r['rapidity'] for r in tables['second_roots']]
        close(state['energy'], 2*holes-sum(1/(a*a+.25) for a in first))
        for i, a in enumerate(first):
            close(phase(a, 1)**length, math.prod(phase(a-b, 2) for j, b in enumerate(first) if i != j)/
                  math.prod(phase(a-b, 1) for b in second))
        for a in second: close(math.prod(phase(a-b, 1) for b in first), 1)
        if holes: close(state['momentum_index'], 0)
        else: close(cmath.exp(1j*state['p']), (-1)**(length-1)*math.prod(phase(a, 1) for a in first))
    for table in counts:
        key = 'k' if table == 'charge_roots' else 'rapidity' if 'roots' in table else None
        if key and any(a[key] >= b[key] for a, b in zip(tables[table], tables[table][1:])):
            raise ValueError('Unordered or duplicate roots')
    return meta, tables


def load_cases(solver_dir=None):
    cases, exports = {}, {}
    for name, case in CASES.items():
        model, length, up, down, c = case
        _, columns = schemas(case)
        if solver_dir:
            program = 'bethe-gaudin-yang-pbc' if model == 'gy' else 'bethe-tj-pbc'
            args = [Path(solver_dir)/program, str(up+down if model == 'gy' else length), '--sz', str((up-down)/2)]
            args += ['--length', str(length), '--c', str(c)] if model == 'gy' else ['--particles', str(up+down)]
            texts = capture_csv_tables(args, columns)
        else:
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in columns}
        cases[name] = read_exports(texts, case)
        exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items(): (DATA/f'{name}-{table}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-fermion-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    einfty = math.pi**2*6*(6**2-1)/(3*6**2)
    energy = lambda name: cases[name][1]['states'][0]['energy']
    fig, ax = plt.subplots(figsize=(7, 4.2), layout='constrained')
    ax.semilogx(COUPLINGS[1:], [energy(f'gy-c{c}')/einfty for c in COUPLINGS[1:]], 'o-', color='#0072B2', label='Balanced interacting branch')
    ax.axhline(energy('gy-c0')/einfty, color='#009E73', linestyle=':', label='Free, 3 up + 3 down')
    ax.axhline(1, color='0.4', linestyle='--', label='Strong-coupling branch limit')
    ax.axhline(energy('gy-polarized')/einfty, color='#D55E00', linestyle=':', label='Polarized periodic free shell')
    ax.set(xlabel=r'$\gamma=c/(N/\ell)$', ylabel=r'$E/E_\infty$', title=r'Gaudin-Yang: N=6, $\ell=6$', ylim=(0, 1.16))
    ax.legend(fontsize=9, loc='center left'); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'gy-coupling.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4), layout='constrained')
    for c, color in [(.01, '#009E73'), (1, '#0072B2'), (1000, '#D55E00')]:
        tables = cases[f'gy-c{c}'][1]
        axes[0].plot([r['quantum_number'] for r in tables['charge_roots']], [r['k'] for r in tables['charge_roots']], 'o-', color=color, label=f'c={c}')
        axes[1].plot([r['quantum_number'] for r in tables['spin_roots']], [r['rapidity']/c for r in tables['spin_roots']], 'o-', color=color, label=f'c={c}')
    axes[0].set(xlabel='Charge label I', ylabel=r'$k$', title='Physical momenta determine E')
    axes[1].set(xlabel='Spin label J', ylabel=r'$\lambda/c$', title='Spin rapidities are auxiliary')
    axes[0].set_xticks([i-.5 for i in range(-2, 4)])
    axes[1].set_xticks([-1, 0, 1])
    axes[1].set_yscale('symlog', linthresh=.5)
    axes[1].set_yticks([-100, -10, -1, 0, 1, 10, 100])
    for ax in axes: ax.legend(fontsize=9); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'gy-roots.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    for family, label, color in [('balanced', 'Balanced supported sectors', '#0072B2'), ('polarized', 'Fully polarized free fermions', '#D55E00')]:
        axes[0].plot([n/16 for n in FILLINGS], [energy(f'tj-{family}-n{n}')/16 for n in FILLINGS], 'o--', color=color, label=label)
    axes[0].set(xlabel=r'$N_e/L$', ylabel='Energy per site (t=1)', title='t-J: selected populations on L=16')
    axes[0].legend(fontsize=8)
    tables = cases['tj-balanced-n14'][1]
    for table, label, color in [('first_roots', r'$\lambda$: 9 first-level roots', '#0072B2'), ('second_roots', r'$\mu$: 2 hole roots', '#D55E00')]:
        rows = tables[table]
        axes[1].plot([r['quantum_number'] for r in rows], [r['rapidity'] for r in rows], 'o-', color=color, label=label)
    axes[1].set(xlabel='Bethe label at the corresponding level', ylabel='Rapidity', title='Two holes: L=16, 7 up + 7 down')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'tj-filling.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path, help='regenerate with both native frontends')
    parser.add_argument('--check', action='store_true', help='validate exports without plotting packages')
    args = parser.parse_args()
    cases = load_cases(args.solver_dir)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} finite fermion calculations.')


if __name__ == '__main__': main()
