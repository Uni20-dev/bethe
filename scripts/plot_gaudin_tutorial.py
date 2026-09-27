#!/usr/bin/env python3
"""Richardson pairing and central-spin tutorials from regular-variable exports."""

import argparse
import itertools
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'docs/tutorials/data'
FIGURES = ROOT/'docs/tutorials/figures'
COUPLINGS = [('0', 0), ('0.1', .1), ('0.3', .3), ('0.5', .5),
             ('collision', 2/3), ('1', 1), ('2', 2), ('4', 4)]
FIELDS = (0, .25, .5, 1, 2, 4)
LEVELS = (0, 1, 2, 3)
BATH = (1, .7, .3)
# mode, levels/couplings, g/B, pairs/Sz, blocked indices
CASES = {f'rich-{branch}-g{tag}': ('rich', LEVELS, g, m, blocked)
         for tag, g in COUPLINGS for branch, m, blocked in [('paired', 2, ()), ('blocked', 1, (1, 2))]}
CASES.update({
    'rich-two': ('rich', (0, 1), 1, 1, ()), 'rich-full': ('rich', LEVELS, 1, 4, ()),
    'rich-empty': ('rich', LEVELS, 1, 0, (1,)), 'rich-free': ('rich', LEVELS, 0, 1, (1,)),
    'rich-shifted': ('rich', (5, 6, 7, 8), 1, 2, ()),
    'rich-scaled': ('rich', (0, 2, 4, 6), 2, 2, ()),
    'rich-irregular': ('rich', (-2, -.8, .1, .9, 2, 3.5), 1, 2, (1,)),
    'rich-block-order': ('rich', LEVELS, 1, 1, (2, 1)),
    'rich-weak': ('rich', LEVELS, 1e-6, 2, ())})
CASES.update({f'central-s{sz}-b{b}': ('central', BATH, b, sz, ())
              for sz in (-1, 0, 1) for b in FIELDS})
CASES.update({f'central-deriv-{tag}': ('central', BATH, b, 0, ()) for tag, b in
              [('minus', .99), ('plus', 1.01), ('half-minus', .995), ('half-plus', 1.005)]})
CASES.update({
    'central-two': ('central', (1,), 1, 0, ()), 'central-two-zero': ('central', (1,), 0, 0, ()),
    'central-reversed': ('central', BATH, -1, -1, ()),
    'central-permuted': ('central', (.3, 1, .7), 1, 0, ()),
    'central-scaled': ('central', (2, 1.4, .6), 2, 0, ()),
    'central-up': ('central', BATH, 1, 2, ()), 'central-down': ('central', BATH, 1, -2, ()),
    'central-isolated': ('central', (), 1, -.5, ()),
    'central-mixed-zero': ('central', (1, -2, .3, -.1), 0, .5, ()),
    'central-mixed': ('central', (1, -2, .3, -.1), 2, .5, ()),
    'central-high': ('central', BATH, 100, 0, ())})
END = 'iterations stages rejected_stages converged status'.split()


def schema(case):
    if case[0] == 'rich':
        return {'states': 'state_id energy requested_g reached_g residual target_residual pair_number_error'.split()+END,
                'variables': 'state_id index epsilon blocked y'.split()}
    return {'states': 'state_id sz energy requested_field reached_field residual number_error'.split()+END,
            'variables': 'state_id spin coupling v'.split()}


def close(actual, expected, tol=5e-11):
    if abs(actual-expected) > tol*max(1, abs(expected)):
        raise ValueError(f'Inconsistent variable/energy: {actual} != {expected}')


def read_exports(texts, case):
    mode, inputs, coupling, sector, blocked = case
    rich = mode == 'rich'
    expected = {'Program': 'bethe-richardson' if rich else 'bethe-central-spin', 'Precision': 'fp64',
                'Hamiltonian': ('sum epsilon_i*n_i-g*sum_ij b_i^dagger*b_j (including i=j)' if rich else
                                'B*S0^z+sum_j A_j*S0.Sj (all spins 1/2)'),
                'Calculation': ('lowest state in the specified pair/blocked sector' if rich else
                                'lowest state in the specified total-Sz sector')}
    meta, tables = None, {}
    for name, columns in schema(case).items():
        current, rows = read_csv_export(texts[name], columns, expected,
                                       row_status='status' if name == 'states' else None)
        if meta is not None and current != meta: raise ValueError('Mismatched run provenance')
        meta = current; tables[name] = rows
    if len(tables['states']) != 1: raise ValueError('Missing/duplicate state')
    raw = tables['states'][0]
    if raw['state_id'] != '0' or raw['converged'] != 'true': raise ValueError('Wrong state/convergence')
    row = {k: finite_number(raw, k) for k in schema(case)['states'] if k not in ('converged', 'status')}
    tolerance = finite_number(meta, 'Residual tolerance')
    if tolerance <= 0 or not 0 <= row['residual'] <= tolerance: raise ValueError('Invalid backward residual')
    for key in ('iterations', 'stages', 'rejected_stages'):
        if not row[key].is_integer() or row[key] < 0: raise ValueError('Invalid work counter')
    if row['rejected_stages'] > row['stages']: raise ValueError('Invalid rejection count')
    close(row['energy'], finite_number(meta, 'Total energy'))
    close(row['residual'], finite_number(meta, 'Reached backward residual'))
    for key in ('requested_g', 'reached_g') if rich else ('requested_field', 'reached_field'):
        close(row[key], coupling)
    for key in ('Requested coupling g', 'Reached coupling g', 'Energy evaluated at g') if rich else (
            'Requested field B', 'Reached field B', 'Energy evaluated at B'):
        close(finite_number(meta, key), coupling)
    variables = tables['variables']
    if len(variables) != len(inputs)+(not rich): raise ValueError('Missing variables/input levels')
    for i, variable in enumerate(variables):
        if variable['state_id'] != '0' or variable['index' if rich else 'spin'] != str(i):
            raise ValueError('Wrong variable join/order')
    if rich:
        active = [i for i in range(len(inputs)) if i not in blocked]
        l, m, g = len(active), sector, coupling
        for key, value in [('Levels', len(inputs)), ('Unblocked levels', l), ('Pairs', m),
                           ('Blocked levels', len(blocked)), ('Fermions', 2*m+len(blocked))]:
            if meta[key] != str(value): raise ValueError('Wrong pairing sector')
        if not 0 <= row['target_residual'] <= tolerance: raise ValueError('Target g not converged')
        close(row['target_residual'], finite_number(meta, 'Target backward residual'))
        y = []
        for i, v in enumerate(variables):
            close(finite_number(v, 'epsilon'), inputs[i])
            if v['blocked'] != ('true' if i in blocked else 'false'): raise ValueError('Wrong blocked index')
            if i in blocked:
                if v['y'] != '': raise ValueError('Blocked level has an invented variable')
            else: y.append(finite_number(v, 'y'))
        close(sum(y), m)
        if abs(row['pair_number_error']) > 1e-11: raise ValueError('Wrong pair number')
        energies = [2*inputs[i] for i in active]
        close(row['energy'], sum(e*v for e, v in zip(energies, y))-g*m*(l-m+1)+sum(inputs[i] for i in blocked))
        coefficient, linear = g, 1
        coords, values = energies, y
    else:
        n, field = len(inputs)+1, coupling
        close(row['sz'], sector)
        if meta['Bath spins'] != str(len(inputs)) or meta['Up spins'] != str(int(n/2+sector)):
            raise ValueError('Wrong spin sector')
        if meta['Spin reversed'] != ('yes' if field < 0 else 'no'): raise ValueError('Wrong spin frame')
        polarized = abs(sector) == n/2
        for i, v in enumerate(variables):
            if i == 0:
                if v['coupling'] != '': raise ValueError('Invented central self coupling')
            else: close(finite_number(v, 'coupling'), inputs[i-1])
        if polarized:
            if any(v['v'] != '' for v in variables): raise ValueError('Analytic state has invented variables')
            close(row['energy'], sum(inputs)/4+(field/2 if sector > 0 else -field/2))
            row['values'] = []
            return meta, row
        scale, h = max(map(abs, inputs)), abs(field)
        coefficient, linear = scale/(2*h+scale), 2*h/(2*h+scale)
        coords = [0]+[scale/a for a in inputs]
        values = [finite_number(v, 'v') for v in variables]
        m = n/2+(sector if field >= 0 else -sector)
        close(sum(values), m*linear)
        if abs(row['number_error']) > 1e-11: raise ValueError('Wrong up-spin constraint')
        close(row['energy'], (h+scale/2)*values[0]-h/2+sum(inputs)/4)
    # Check the model's quadratic equations independently of exported diagnostics.
    for i, (q, v) in enumerate(zip(coords, values)):
        terms = [coefficient*(v-w)/(q-r) for j, (r, w) in enumerate(zip(coords, values)) if j != i]
        denominator = 1+abs(v*(v-linear))+sum(abs(coefficient/(q-r))*(abs(v)+abs(w))
            for j, (r, w) in enumerate(zip(coords, values)) if j != i)
        if abs(v*(v-linear)-sum(terms))/denominator > max(8*tolerance, 1e-13):
            raise ValueError('Variables violate quadratic equations')
    row['values'] = values
    return meta, row


def load_cases(solver_dir=None):
    cases, exports = {}, {}
    for name, case in CASES.items():
        mode, inputs, coupling, sector, blocked = case
        if solver_dir:
            program = 'richardson' if mode == 'rich' else 'central-spin'
            command = [Path(solver_dir)/f'bethe-{program}']
            if mode == 'rich':
                command += ['--levels', ','.join(map(str, inputs)), '--g', str(coupling), '--pairs', str(sector)]
                if blocked: command += ['--blocked', ','.join(map(str, blocked))]
            else: command += ['--couplings', ','.join(map(str, inputs)), '--field', str(coupling), '--sz', str(sector)]
            texts = capture_csv_tables(command, schema(case))
        else: texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schema(case)}
        cases[name] = read_exports(texts, case); exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items(): (DATA/f'{name}-{table}.csv').write_text(text)
    return cases


def oracle_check(cases):
    """Independent occupation/spin matrices, with no Gaudin equations (NumPy)."""
    import numpy as np
    for name, case in CASES.items():
        mode, inputs, coupling, sector, blocked = case
        if mode == 'rich':
            active = [i for i in range(len(inputs)) if i not in blocked]
            basis = [frozenset(x) for x in itertools.combinations(active, sector)]
            index = {x: i for i, x in enumerate(basis)}
            h = np.zeros((len(basis), len(basis)))
            for j, occupied in enumerate(basis):
                h[j, j] = 2*sum(inputs[i] for i in occupied)+sum(inputs[i] for i in blocked)-coupling*sector
                for old in occupied:
                    for new in set(active)-occupied: h[index[occupied-{old}|{new}], j] -= coupling
        else:
            n = len(inputs)+1
            basis = [w for w in range(1 << n) if w.bit_count() == int(n/2+sector)]
            index = {w: i for i, w in enumerate(basis)}
            h = np.zeros((len(basis), len(basis)))
            for j, word in enumerate(basis):
                central = (word & 1)-.5
                h[j, j] = coupling*central
                for bath, a in enumerate(inputs, 1):
                    spin = ((word >> bath) & 1)-.5
                    h[j, j] += a*central*spin
                    if spin != central: h[index[word ^ 1 ^ (1 << bath)], j] += a/2
        ev, vectors = np.linalg.eigh(h)
        close(cases[name][1]['energy'], ev[0])
        if name == 'central-s0-b1':
            sz0 = sum(abs(vectors[j, 0])**2*((word & 1)-.5) for j, word in enumerate(basis))
    energy = lambda name: cases[name][1]['energy']
    derivative = (energy('central-deriv-half-plus')-energy('central-deriv-half-minus'))/.01
    close(derivative, sz0, 2e-6)
    print('All independent pairing/spin matrix minima and central polarization check passed.')


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-gaudin-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    energy = lambda name: cases[name][1]['energy']
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    gs = [g for _, g in COUPLINGS]
    gaps = [energy(f'rich-blocked-g{tag}')-energy(f'rich-paired-g{tag}') for tag, _ in COUPLINGS]
    axes[0].plot(gs, gaps, 'o-', color='#0072B2')
    axes[0].set(xlabel='Attractive coupling g', ylabel='Blocked-sector energy − paired energy', title='Same four fermions; specified blocked levels')
    for i in range(4):
        axes[1].plot(gs, [cases[f'rich-paired-g{tag}'][1]['values'][i] for tag, _ in COUPLINGS], 'o-', markersize=3, label=f'Level {i}')
    axes[1].axvline(2/3, color='0.5', linestyle='--', label='Pair-root collision')
    axes[1].set(xlabel='Attractive coupling g', ylabel=r'Regular variable $y_i$ (not occupation)', title='Smooth variables through g=2/3')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'richardson-pairing.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    for sz in (-1, 0, 1):
        axes[0].plot(FIELDS, [energy(f'central-s{sz}-b{b}')+b/2 for b in FIELDS], 'o-', label=f'Total Sz={sz}')
    axes[0].set(xlabel='Central field B', ylabel=r'$E_{S^z}(B)+B/2$', title='A=(1, 0.7, 0.3): sector minima')
    slopes = [(energy(f'central-s0-b{b}')-energy(f'central-s0-b{a}'))/(b-a) for a, b in zip(FIELDS, FIELDS[1:])]
    axes[1].errorbar([(a+b)/2 for a, b in zip(FIELDS, FIELDS[1:])], slopes,
                     xerr=[(b-a)/2 for a, b in zip(FIELDS, FIELDS[1:])], fmt='o', capsize=3,
                     color='#0072B2', label='Interval-averaged central Sz')
    derivative = (energy('central-deriv-half-plus')-energy('central-deriv-half-minus'))/.01
    axes[1].plot([1], [derivative], '*', markersize=10, color='#D55E00', label='Local derivative at B=1')
    axes[1].axhline(-.5, color='0.5', linestyle='--', label='Fully down central spin')
    axes[1].set(xlabel='Central field B (bars show averaging interval)', ylabel=r'$\Delta E_{S^z=0}/\Delta B$', title='Energy slopes, not Gaudin variables')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'central-spin-field.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--oracle', action='store_true', help='Check small independent matrices (NumPy)')
    args = parser.parse_args(); cases = load_cases(args.solver_dir)
    if args.oracle: oracle_check(cases)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} pairing/central-spin calculations.')


if __name__ == '__main__': main()
