#!/usr/bin/env python3
"""Universal Kondo response and finite integrable-ladder field benchmarks."""

import argparse
import itertools
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'docs/tutorials/data'
FIGURES = ROOT/'docs/tutorials/figures'
KONDO_FIELDS = (0, .0001, .01, .1, .5, .99, 1, 1.01, 2, 10, 100, 10000, 1000000)
LADDER_FIELDS = (0, .5, 1, 1.5, 3, 5, 7, 8.5, 9, 10)
# Kondo: b, T_B. Ladder: L, J_r, h, fixed M or None, fixed N_s or None.
KONDO = {f'kondo-b{b}': (b, 1) for b in KONDO_FIELDS}
KONDO.update({'kondo-negative': (-2, 1), 'kondo-scaled': (6, 3),
              'kondo-deriv1-minus': (.999, 1), 'kondo-deriv1-plus': (1.001, 1),
              'kondo-deriv2-minus': (1.999, 1), 'kondo-deriv2-plus': (2.001, 1)})
LADDER = {f'ladder-m{m}': (6, 5, 0, m, None) for m in range(7)}
LADDER.update({f'ladder-h{h}': (6, 5, h, None, None) for h in LADDER_FIELDS})
LADDER.update({
    'ladder-descendant': (6, 0, 0, None, 4),
    'ladder-desc-plus': (6, 0, .125, None, 4), 'ladder-desc-minus': (6, 0, -.125, None, 4),
    'ladder-desc-rung': (6, 1, 0, None, 4),
    'ladder-fixed-field': (6, 5, 3, 2, None), 'ladder-reversed': (6, 5, -3, -2, None),
    'ladder-two': (2, 0, 0, None, None), 'ladder-four': (4, 0, 0, None, None),
    'ladder-triplet-even': (6, 5, 0, None, 5), 'ladder-triplet-odd': (5, 5, 0, None, 4)})
K_COLUMNS = 'energy_change magnetization zero_field_susceptibility magnetization_error scaled_energy_error series_terms lobes evaluations converged status'.split()
L_COLUMNS = 'state_id singlets energy magnetization selected descendant momentum_index residual iterations branches tableaux converged status'.split()
SCHEMAS = {'states': L_COLUMNS, 'representations': 'state_id component population highest_weight'.split(),
           'roots': 'state_id level index quantum_number rapidity'.split()}
A0 = 1/math.sqrt(2*math.pi*math.e)


def close(actual, expected, tol=5e-11):
    if abs(actual-expected) > tol*max(1, abs(expected)):
        raise ValueError(f'Inconsistent observable/identity: {actual} != {expected}')


def integer(row, key):
    value = finite_number(row, key)
    if not value.is_integer(): raise ValueError(f'Noninteger {key}')
    return int(value)


def read_kondo(text, case):
    field, scale = case
    meta, rows = read_csv_export(text, K_COLUMNS, {
        'Program': 'bethe-kondo-response', 'Precision': 'fp64',
        'Calculation': 'T=0; isotropic antiferromagnetic spin-1/2 single-channel scaling limit',
        'Field convention': 'H_field=-b*(S_imp^z+S_host^z); equal g factors',
        'Energy convention': 'Delta E_imp=E_imp(b)-E_imp(0); E_imp=E_with-E_clean_host',
        'Response convention': 'Impurity-induced total magnetization, including host change'})
    close(finite_number(meta, 'Full Zeeman splitting'), field)
    close(finite_number(meta, 'Universal scale T_B'), scale)
    if len(rows) != 1 or rows[0]['converged'] != 'true': raise ValueError('Missing/convergence response')
    row = {k: finite_number(rows[0], k) for k in K_COLUMNS[:-2]}
    tolerance = finite_number(meta, 'Absolute M and Delta E/|b| tolerance')
    if tolerance <= 0: raise ValueError('Invalid tolerance')
    for key in ('magnetization_error', 'scaled_energy_error'):
        if not 0 <= row[key] <= tolerance: raise ValueError('Response error estimate outside tolerance')
    for key in ('series_terms', 'lobes', 'evaluations'):
        if row[key] < 0 or not row[key].is_integer(): raise ValueError('Invalid work counter')
    close(row['zero_field_susceptibility'], A0/scale)
    if field == 0:
        if any(row[k] != 0 for k in K_COLUMNS[:-2] if k != 'zero_field_susceptibility'):
            raise ValueError('Nonzero response/work at zero field')
    elif not (-abs(field)/2 < row['energy_change'] < 0 and 0 < row['magnetization']/field and abs(row['magnetization']) < .5):
        raise ValueError('Wrong energy/magnetization bounds or field sign')
    return meta, row


def read_ladder(texts, case):
    n, rung, field, fixed_m, fixed_s = case
    calculation = ('fixed Sz and singlet count' if fixed_s is not None else 'fixed Sz sector minimum') if fixed_m is not None else (
        'fixed singlet-count sector' if fixed_s is not None else 'global ground state')
    expected = {'Program': 'bethe-ladder-pbc', 'Precision': 'fp64', 'Rungs': str(n), 'Calculation': calculation,
                'Hamiltonian': 'H=sum[S.S_next+T.T_next+4(S.S_next)(T.T_next)]+J_r*sum S.T-h*sum(Sz+Tz)'}
    meta, tables = None, {}
    for table, columns in SCHEMAS.items():
        current, rows = read_csv_export(texts[table], columns, expected, row_status='status' if table == 'states' else None)
        if meta is not None and meta != current: raise ValueError('Mismatched table provenance')
        meta = current; tables[table] = rows
    close(finite_number(meta, 'Rung coupling J_r'), rung)
    close(finite_number(meta, 'Magnetic field h'), field)
    if len(tables['states']) != 1 or len(tables['representations']) != 4: raise ValueError('Missing state/representation')
    raw = tables['states'][0]
    if raw['state_id'] != '0' or raw['selected'] != 'true' or raw['converged'] != 'true':
        raise ValueError('Wrong selected state')
    row = {k: finite_number(raw, k) for k in L_COLUMNS if k not in ('selected', 'descendant', 'converged', 'status')}
    for key in ('state_id', 'singlets', 'magnetization', 'momentum_index', 'iterations', 'branches', 'tableaux'):
        row[key] = integer(row, key)
    if any(row[k] < 0 for k in ('singlets', 'iterations', 'branches', 'tableaux')): raise ValueError('Invalid counts')
    tolerance = finite_number(meta, 'Residual tolerance')
    if tolerance <= 0 or not 0 <= row['residual'] <= tolerance: raise ValueError('Unconverged selected branch')
    pop, shape = [], []
    for i, rep in enumerate(tables['representations']):
        if rep['state_id'] != '0' or rep['component'] != str(i): raise ValueError('Wrong representation join')
        pop.append(integer(rep, 'population')); shape.append(integer(rep, 'highest_weight'))
    if (sum(pop) != n or sum(shape) != n or min(pop+shape) < 0 or shape != sorted(shape, reverse=True) or
            any(sum(shape[:i]) < sum(sorted(pop, reverse=True)[:i]) for i in (1, 2, 3))):
        raise ValueError('Invalid populations or highest-weight dominance')
    if row['singlets'] != pop[0] or row['magnetization'] != pop[1]-pop[3]: raise ValueError('Wrong physical sector')
    if fixed_m is not None and row['magnetization'] != fixed_m: raise ValueError('Wrong fixed Sz')
    if fixed_s is not None and row['singlets'] != fixed_s: raise ValueError('Wrong fixed singlets')
    if raw['descendant'] != ('true' if sorted(pop, reverse=True) != shape else 'false'):
        raise ValueError('Incorrect descendant classification')
    for key, values in [('Populations (s,t+,t0,t-)', pop), ('Highest-weight rows', shape)]:
        if list(map(int, meta[key].split(','))) != values: raise ValueError('Metadata/representation mismatch')
    counts = [sum(shape[a:]) for a in range(1, 4)]
    if len(tables['roots']) != sum(counts): raise ValueError('Wrong nested root count')
    levels, labels, offset = [], [], 0
    for a, count in enumerate(counts, 1):
        roots = tables['roots'][offset:offset+count]; offset += count
        values, numbers = [], []
        for j, root in enumerate(roots):
            if (root['state_id'], root['level'], root['index']) != ('0', str(a), str(j)):
                raise ValueError('Wrong root join/level/index')
            values.append(finite_number(root, 'rapidity')); numbers.append(finite_number(root, 'quantum_number'))
        if any(x >= y for x, y in zip(values, values[1:])): raise ValueError('Unordered roots')
        if count:
            shift = numbers[0]+(count-1)/2
            previous = n if a == 1 else counts[a-2]
            following = counts[a] if a < 3 else 0
            if shift not in ((-.5, .5) if (previous+following) % 2 else (0,)):
                raise ValueError('Wrong packed-sea displacement')
            for j, number in enumerate(numbers): close(number, j-(count-1)/2+shift)
        levels.append(values); labels.append(numbers)
    def theta(value, width): return 2*math.atan(value/width)
    for a, roots in enumerate(levels):
        for j, value in enumerate(roots):
            lhs = n*theta(value, .5) if a == 0 else sum(theta(value-v, .5) for v in levels[a-1])
            lhs -= sum(theta(value-v, 1) for v in roots)
            if a < 2: lhs += sum(theta(value-v, .5) for v in levels[a+1])
            close((lhs-2*math.pi*labels[a][j])/n, 0, max(8*tolerance, 1e-13))
    perm = n-sum(1/(v*v+.25) for v in levels[0])
    close(perm, finite_number(meta, 'Permutation energy'))
    close(row['energy'], perm-n/4+rung*(n/4-pop[0])-field*row['magnetization'])
    close(row['energy'], finite_number(meta, 'Total energy'))
    close(row['energy']/n, finite_number(meta, 'Energy per rung'))
    momentum = (n*counts[0]-sum(2*i for sea in labels for i in sea))/2
    if not momentum.is_integer() or row['momentum_index'] != int(momentum) % n: raise ValueError('Wrong rung momentum')
    if integer(meta, 'Reflected momentum index') != (-row['momentum_index']) % n: raise ValueError('Wrong reflected momentum')
    row.update(populations=pop, shape=shape, levels=levels, labels=labels, permutation_energy=perm)
    return meta, row


def load_cases(solver_dir=None):
    results, exports = {}, {}
    for name, case in {**KONDO, **LADDER}.items():
        kondo = name in KONDO
        tables = ['response'] if kondo else SCHEMAS
        if solver_dir:
            program = 'kondo-response' if kondo else 'ladder-pbc'
            args = [Path(solver_dir)/f'bethe-{program}']
            if kondo: args += ['--field', str(case[0]), '--scale', str(case[1])]
            else:
                n, rung, field, m, s = case
                args += [str(n), '--rung', str(rung), '--field', str(field)]
                if m is not None: args += ['--sz', str(m)]
                if s is not None: args += ['--singlets', str(s)]
            texts = capture_csv_tables(args, tables)
        else: texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in tables}
        results[name] = read_kondo(texts['response'], case) if kondo else read_ladder(texts, case)
        exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items(): (DATA/f'{name}-{table}.csv').write_text(text)
    return results


def ladder_oracle(cases):
    """Minimize independent permutation matrices over every allowed population."""
    import numpy as np
    cache = {}
    def permutation_minimum(pop):
        shape = tuple(sorted(pop, reverse=True))
        if shape not in cache:
            def words(prefix, counts):
                if not any(counts): yield tuple(prefix)
                for color, count in enumerate(counts):
                    if count:
                        left = list(counts); left[color] -= 1
                        yield from words(prefix+[color], left)
            basis = list(words([], shape)); index = {w: i for i, w in enumerate(basis)}
            n = sum(shape); h = np.zeros((len(basis), len(basis)))
            for col, w in enumerate(basis):
                for j in range(n):
                    v = list(w); v[j], v[(j+1) % n] = v[(j+1) % n], v[j]
                    h[index[tuple(v)], col] += 1
            cache[shape] = np.linalg.eigvalsh(h)[0]
        return cache[shape]
    for name, (n, rung, field, fixed_m, fixed_s) in LADDER.items():
        energies = []
        for s, plus, zero in itertools.product(range(n+1), repeat=3):
            minus = n-s-plus-zero
            if minus < 0 or (fixed_s is not None and fixed_s != s): continue
            m = plus-minus
            if fixed_m is not None and fixed_m != m: continue
            energies.append(permutation_minimum((s, plus, zero, minus))-n/4+rung*(n/4-s)-field*m)
        close(cases[name][1]['energy'], min(energies))
    print(f'All ladder minima match independent color matrices ({len(cache)} cached population shapes).')


def ladder_envelope(cases, field):
    energies = [cases[f'ladder-m{m}'][1]['energy']-field*m for m in range(7)]
    m = min(range(7), key=energies.__getitem__)
    return energies[m], m


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-kondo-ladder-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fields = [b for b in KONDO_FIELDS if b > 0]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    axes[0].semilogx(fields, [cases[f'kondo-b{b}'][1]['magnetization'] for b in fields], 'o-', color='#0072B2')
    axes[0].axhline(.5, color='0.5', linestyle='--', label='Saturation limit 1/2')
    axes[0].set(xlabel=r'$|b|/T_B$', ylabel=r'$M_{\rm imp}$ (positive field)', title='Impurity-induced total magnetization')
    axes[1].semilogx(fields, [cases[f'kondo-b{b}'][1]['energy_change']/b for b in fields], 'o-', color='#D55E00')
    axes[1].axhline(-.5, color='0.5', linestyle='--', label='High-field limit −1/2')
    axes[1].set(xlabel=r'$|b|/T_B$', ylabel=r'$\Delta E_{\rm imp}/|b|$', title='Subtracted energy, not absolute binding')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'kondo-response.svg'); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    h = [j/200 for j in range(2001)]
    envelope = [ladder_envelope(cases, v) for v in h]
    for m in range(7):
        axes[0].plot(h, [(cases[f'ladder-m{m}'][1]['energy']-v*m)/6 for v in h], color='0.75', linewidth=.8)
    axes[0].plot(h, [e/6 for e, _ in envelope], color='#0072B2', label='Minimum over fixed-Sz energies')
    axes[0].plot(LADDER_FIELDS, [cases[f'ladder-h{v}'][1]['energy']/6 for v in LADDER_FIELDS], 'o', markersize=4, color='#D55E00', label='Independent field runs')
    axes[0].set(xlabel='Field h', ylabel='Energy per rung', title='L=6, Jr=5: linear sector branches')
    axes[1].step(h, [m/6 for _, m in envelope], where='post', color='#0072B2', label='Sector-envelope magnetization')
    axes[1].plot(LADDER_FIELDS, [cases[f'ladder-h{v}'][1]['magnetization']/6 for v in LADDER_FIELDS], 'o', markersize=4, color='#D55E00', label='One native minimizing representative')
    axes[1].set(xlabel='Field h', ylabel='Total M / number of rungs', title='Finite-size magnetization steps')
    for ax in axes: ax.legend(fontsize=7); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'ladder-field.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--oracle', action='store_true', help='Independent ladder matrices (NumPy)')
    args = parser.parse_args(); cases = load_cases(args.solver_dir)
    if args.oracle: ladder_oracle(cases)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} Kondo/ladder calculations.')


if __name__ == '__main__': main()
