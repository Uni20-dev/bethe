#!/usr/bin/env python3
"""SU(kappa) fermion and Bose-Fermi tutorials, with native reduction benchmarks."""

import argparse
import cmath
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg
from plot_fermion_tutorial import read_exports as read_gy, schemas as gy_schemas

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'docs/tutorials/data'
FIGURES = ROOT/'docs/tutorials/figures'
COUPLINGS = (.001, .01, .1, 1, 10, 100, 1000)
# model, populations, circumference, c
CASES = {f'sun-c{c}': ('sun', (3, 3, 3), 9, c) for c in COUPLINGS}
CASES.update({
    'sun-free': ('sun', (3, 3, 3), 9, 0), 'sun-polarized': ('sun', (9, 0, 0), 9, 1),
    'sun-imbalanced': ('sun', (5, 3, 1), 9, 1), 'sun-permuted': ('sun', (1, 0, 5, 3), 9, 1),
    'sun-scaled': ('sun', (3, 3, 3), 18, .5), 'sun-two': ('sun', (3, 3), 6, 1),
    'sun-singlet': ('sun', (1, 1, 1), 3, 1), 'sun-four': ('sun', (3, 1, 1, 1), 6, 1),
    'sun-free-shell': ('sun', (2, 0, 2), 1, 0), 'sun-vacuum': ('sun', (0, 0, 0), 1, 1)})
CASES.update({f'bf-c{c}': ('bf', (2, 3), 5, c) for c in COUPLINGS})
CASES.update({
    'bf-free': ('bf', (2, 3), 5, 0), 'bf-fermions': ('bf', (0, 5), 5, 1),
    'bf-bosons': ('bf', (5, 0), 5, 1), 'bf-one-fermion': ('bf', (4, 1), 5, 1),
    'bf-scaled': ('bf', (2, 3), 10, .5), 'bf-free-shell': ('bf', (2, 2), 1, 0),
    'bf-vacuum': ('bf', (0, 0), 1, 1), 'bf-contact': ('bf', (1, 1), 1, math.pi),
    'ref-ll3': ('ll', (3,), 3, 1), 'ref-ll5': ('ll', (5,), 5, 1),
    'ref-gy6': ('gy', (3, 3), 6, 1)})
PROGRAMS = {'sun': 'sun-fermions', 'bf': 'bose-fermi', 'll': 'lieb-liniger', 'gy': 'gaudin-yang'}
END = ['iterations', 'converged', 'status']


def schemas(case):
    model, pop, length, c = case
    if model == 'gy': return gy_schemas(('gy', length, *pop, c))
    free = (c == 0 or (sum(n > 0 for n in pop) <= 1 if model == 'sun' else pop[0] == 0)) if model != 'll' else False
    if model == 'sun':
        result = {'states': 'state_id energy requested_c reached_c momentum_index p residual target_residual iterations stages converged status'.split(),
                  'components': 'state_id component particles nesting_rank'.split()}
        result.update({'free_modes': 'state_id component mode k'.split()} if free else
                      {'roots': 'state_id level index quantum_number rapidity'.split()})
    elif model == 'bf':
        result = {'states': 'state_id energy momentum_index momentum residual'.split()+END}
        result.update({'free_modes': 'state_id species index mode k'.split()} if free else
                      {'charge_roots': 'state_id index I k'.split(), 'auxiliary_roots': 'state_id index J lambda'.split()})
    else:
        result = {'states': 'state_id energy gap momentum_index p residual'.split()+END,
                  'roots': 'state_id index quantum_number k'.split()}
    return free, result


def close(actual, expected, tol=5e-11):
    if abs(actual-expected) > tol*max(1, abs(expected)): raise ValueError('Inconsistent root, label or observable')


def phase(x, c):
    return complex(x, c/2)/complex(x, -c/2)


def read_exports(texts, case):
    model, pop, length, c = case
    if model == 'gy': return read_gy(texts, ('gy', length, *pop, c))
    free, columns = schemas(case)
    n = sum(pop)
    branch = ('exact free-fermion ground state' if free else 'odd-population sector ground state') if model == 'sun' else (
        'exact free ground state' if free else 'odd-fermion mixed ground state' if pop[1] else 'Lieb-Liniger reduction') if model == 'bf' else 'ground state'
    expected = {'Program': f'bethe-{PROGRAMS[model]}-pbc', 'Precision': 'fp64', 'Calculation': branch}
    if model == 'bf': expected.update({'Bosons N_b': str(pop[0]), 'Fermions N_f': str(pop[1]),
                                      'Hamiltonian': '-sum d_j^2+2c sum delta; equal masses, equal BB/BF couplings'})
    else: expected.update({'Particles': str(n), 'Units': 'hbar^2/(2m)=1; interaction 2c delta'})
    meta, tables = None, {}
    for table, schema in columns.items():
        current, rows = read_csv_export(texts[table], schema, expected, row_status='status' if table == 'states' else None)
        if meta is not None and current != meta: raise ValueError('Mismatched run provenance')
        meta = current
        parsed = []
        for row in rows:
            if row['state_id'] != '0' or ('converged' in row and row['converged'] != 'true'):
                raise ValueError('Wrong state join/status')
            item = dict(row)
            for key in schema:
                if key in ('species', 'status', 'converged'): continue
                if key in ('nesting_rank', 'gap') and row[key] == '': item[key] = None
                else: item[key] = finite_number(row, key)
            parsed.append(item)
        tables[table] = parsed
    close(finite_number(meta, 'Circumference' if model == 'bf' else 'Length'), length)
    close(finite_number(meta, 'Requested c' if model == 'sun' else 'Contact parameter c' if model == 'bf' else 'c'), c)
    if len(tables['states']) != 1: raise ValueError('Missing/duplicate state')
    state = tables['states'][0]
    close(state['energy'], finite_number(meta, 'Energy' if model == 'bf' else 'Total energy'))
    momentum = state['momentum' if model == 'bf' else 'p']
    close(momentum, 2*math.pi*state['momentum_index']/length)
    if not state['momentum_index'].is_integer() or state['energy'] < 0: raise ValueError('Invalid observable')
    tolerance = finite_number(meta, 'Residual tolerance')
    if tolerance <= 0 or not 0 <= state['residual'] <= tolerance or state['iterations'] < 0:
        raise ValueError('Invalid convergence diagnostics')
    if model == 'sun':
        close(state['requested_c'], c); close(state['reached_c'], c)
        close(finite_number(meta, 'Energy evaluated at c'), c)
        if not 0 <= state['target_residual'] <= tolerance or state['stages'] < 0: raise ValueError('Unfinished continuation')
        order = sorted((i for i, count in enumerate(pop) if count), key=lambda i: -pop[i])
        if meta['Components'] != str(len(pop)) or meta['Occupied components'] != str(len(order)):
            raise ValueError('Wrong component counts')
        if len(tables['components']) != len(pop): raise ValueError('Missing physical component')
        for i, row in enumerate(tables['components']):
            if (row['component'] != i or row['particles'] != pop[i] or
                    row['nesting_rank'] != (order.index(i) if i in order else None)):
                raise ValueError('Wrong original-component/nesting mapping')
    if model == 'll' and state['gap'] is not None: raise ValueError('Ground energy is not a gap')
    if free:
        rows = tables['free_modes']
        if len(rows) != n: raise ValueError('Missing free occupations')
        expected_modes = []
        for component, count in enumerate(pop):
            for i in range(count):
                mode = 0 if model == 'bf' and component == 0 else i-(count-1)//2
                expected_modes.append((component, i, mode))
        for row, (component, i, mode) in zip(rows, expected_modes):
            if model == 'sun':
                if row['component'] != component: raise ValueError('Wrong free component')
            elif row['species'] != ('boson' if component == 0 else 'fermion') or row['index'] != i:
                raise ValueError('Wrong free species/index')
            if row['mode'] != mode: raise ValueError('Wrong free shell')
            close(row['k'], 2*math.pi*mode/length)
        close(state['energy'], sum(r['k']**2 for r in rows)); close(momentum, sum(r['k'] for r in rows))
        return meta, tables
    if model == 'sun':
        counts = [sum(pop[i] for i in order[a:]) for a in range(len(order))]
        if len(tables['roots']) != sum(counts): raise ValueError('Wrong total nested root count')
        levels, offset = [], 0
        for a, count in enumerate(counts):
            rows = tables['roots'][offset:offset+count]; offset += count
            if any(row['level'] != a for row in rows): raise ValueError('Wrong root level')
            check_labels(rows, count, 'quantum_number')
            levels.append([r['rapidity'] for r in rows])
    elif model == 'bf':
        check_labels(tables['charge_roots'], n, 'I')
        check_labels(tables['auxiliary_roots'], pop[0] if pop[1] else 0, 'J')
        levels = [[r['k'] for r in tables['charge_roots']], [r['lambda'] for r in tables['auxiliary_roots']]]
    else:
        check_labels(tables['roots'], n, 'quantum_number')
        levels = [[r['k'] for r in tables['roots']]]
    for roots in levels:
        if any(a >= b for a, b in zip(roots, roots[1:])): raise ValueError('Unordered/duplicate roots')
    charge = levels[0]
    close(state['energy'], sum(k*k for k in charge)); close(momentum, sum(charge))
    for j, k in enumerate(charge):
        rhs = (math.prod(phase(k-v, 2*c) for i, v in enumerate(charge) if i != j)
               if model == 'll' or (model == 'bf' and pop[1] == 0)
               else math.prod(phase(k-v, c) for v in levels[1]))
        close(cmath.exp(1j*k*length), rhs)
    for a in range(1, len(levels)):
        for j, value in enumerate(levels[a]):
            lhs = math.prod(phase(value-v, c) for v in levels[a-1])
            if a+1 < len(levels): lhs *= math.prod(phase(value-v, c) for v in levels[a+1])
            rhs = math.prod(phase(value-v, 2*c) for i, v in enumerate(levels[a]) if i != j) if model == 'sun' else 1
            close(lhs, rhs)
    tables['levels'] = levels
    return meta, tables


def check_labels(rows, count, key):
    if len(rows) != count: raise ValueError('Wrong sea size')
    for i, row in enumerate(rows):
        if row['index'] != i: raise ValueError('Wrong root index')
        close(row[key], i-(count-1)/2)


def load_cases(solver_dir=None):
    cases, exports = {}, {}
    for name, case in CASES.items():
        model, pop, length, c = case
        _, columns = schemas(case)
        if solver_dir:
            args = [Path(solver_dir)/f'bethe-{PROGRAMS[model]}-pbc']
            if model == 'sun': args += ['--populations', ','.join(map(str, pop))]
            elif model == 'bf': args += ['--bosons', str(pop[0]), '--fermions', str(pop[1])]
            else:
                args += [str(sum(pop))]
                if model == 'gy': args += ['--sz', str((pop[0]-pop[1])/2)]
            args += ['--length', str(length), '--c', str(c)]
            texts = capture_csv_tables(args, columns)
        else: texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in columns}
        cases[name] = read_exports(texts, case); exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items(): (DATA/f'{name}-{table}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-multicomponent-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    energy = lambda name: cases[name][1]['states'][0]['energy']
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    limit = math.pi**2*9*80/(3*9**2)
    axes[0].semilogx(COUPLINGS, [energy(f'sun-c{c}')/limit for c in COUPLINGS], 'o-', color='#0072B2', label='3+3+3 interacting sea')
    axes[0].axhline(energy('sun-free')/limit, color='#009E73', linestyle=':', label='Free 3+3+3')
    axes[0].axhline(energy('sun-polarized')/limit, color='0.5', linestyle='--', label='9+0+0 free / strong limit')
    axes[0].set(xlabel=r'$c/(N/\ell)$', ylabel=r'$E/E_\infty$', title=r'SU(3) gas: N=9, $\ell=9$')
    levels = cases['sun-c1'][1]['levels']
    for a, values in enumerate(levels):
        axes[1].plot([j-(len(values)-1)/2 for j in range(len(values))], values, 'o-', label=f'Level {a}: {len(values)} roots')
    axes[1].set(xlabel='Centered label at each level', ylabel='Physical rapidity', title='c=1: charge and two spin seas')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'sun-seas.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    limit = math.pi**2*5*24/(3*5**2)
    axes[0].semilogx(COUPLINGS, [energy(f'bf-c{c}')/limit for c in COUPLINGS], 'o-', color='#0072B2', label='2 bosons + 3 fermions')
    axes[0].axhline(energy('bf-free')/limit, color='#009E73', linestyle=':', label='Free mixture')
    axes[0].axhline(energy('bf-fermions')/limit, color='0.5', linestyle='--', label='Five free fermions / strong limit')
    axes[0].set(xlabel=r'$c/(N/\ell)$', ylabel=r'$E/E_\infty$', title=r'Equal-coupling mixture: N=5, $\ell=5$')
    names = ['bf-fermions', 'bf-c1', 'bf-one-fermion', 'bf-bosons']
    axes[1].bar(['0 + 5', '2 + 3', '4 + 1', '5 + 0'], [energy(name) for name in names], color=['0.5','#0072B2','#D55E00','#009E73'])
    axes[1].set(xlabel='Bosons + fermions (selected compositions)', ylabel='Total energy', title=r'c=1, same N and $\ell$')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=.15); axes[1].grid(axis='y', alpha=.15)
    save_svg(fig, FIGURES/'bf-composition.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args(); cases = load_cases(args.solver_dir)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} multicomponent and reference calculations.')


if __name__ == '__main__': main()
