#!/usr/bin/env python3
"""Finite Potts levels and non-Hermitian XXZ blocks from native exports."""

import argparse
import cmath
import itertools
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'docs/tutorials/data'
FIGURES = ROOT/'docs/tutorials/figures'
POTTS_SIZES = (2, 4, 8, 16, 32, 64)
EVEN = (4, 6, 8, 10, 12)
ODD = (5, 7, 9, 11, 13)
# mode, sites, parameter (Sz at the free endpoint; through-lines in a scan)
CASES = {f'potts-l{n}': ('potts', n, None) for n in POTTS_SIZES}
CASES.update({f'qg-free-n{n}': ('free', n, (n % 2)/2) for n in (2, *EVEN, *ODD)})
CASES.update({'qg-free-plus': ('free', 6, 1), 'qg-free-minus': ('free', 6, -1),
              'qg-free-up': ('free', 6, 3), 'qg-free-down': ('free', 6, -3),
              'qg-scan': ('scan', 8, 4), 'qg-one': ('scan', 8, 6)})
POTTS_COLUMNS = 'branch charge k momentum energy gap x_scaled residual iterations status'.split()
BLOCK_COLUMNS = 'block_id energy block_size zero_occupation modes'.split()
STATE_COLUMNS = 'state_id numbers energy energy_shift gap_from_sea residual iterations converged status'.split()
ROOT_COLUMNS = 'source state_id root_id I z lambda converged'.split()
VELOCITY = 3*math.sqrt(3)/2
BULK = -4/3-2*math.sqrt(3)/math.pi


def close(actual, expected, tolerance=5e-11):
    if abs(actual-expected) > tolerance*max(1, abs(expected)):
        raise ValueError(f'Inconsistent observable/identity: {actual} != {expected}')


def integer(row, key):
    value = finite_number(row, key)
    if not value.is_integer(): raise ValueError(f'Noninteger {key}')
    return int(value)


def schema(case):
    if case[0] == 'potts': return {'levels': POTTS_COLUMNS}
    if case[0] == 'free': return {'blocks': BLOCK_COLUMNS}
    return {'levels': STATE_COLUMNS, 'reference': STATE_COLUMNS, 'roots': ROOT_COLUMNS}


def read_exports(texts, case):
    mode, n, parameter = case
    expected = {'Program': 'bethe-potts-pbc' if mode == 'potts' else 'bethe-xxz-qg-obc',
                'Precision': 'fp64', 'Sites': str(n),
                'Status': {'potts': 'converged', 'free': 'complete', 'scan': 'family converged'}[mode]}
    expected['Hamiltonian'] = {
        'potts': 'H=-sum(X+X^dagger+ZZ^dagger+Z^dagger Z); critical, J=1, PBC',
        'free': 'sum(SxSx+SySy)+i/2*(Sz_1-Sz_N); spin-half J=1',
        'scan': 'sum(SxSx+SySy+Delta SzSz)+i*sqrt(1-Delta^2)/2*(Sz_1-Sz_N); spin-half J=1'}[mode]
    meta, tables = None, {}
    for name, columns in schema(case).items():
        current, rows = read_csv_export(texts[name], columns, expected,
                                       row_status='status' if name in ('levels', 'reference') else None)
        if meta is not None and meta != current: raise ValueError('Mismatched table provenance')
        meta = current
        tables[name] = rows
    if mode == 'potts':
        if (meta['Branches'] != 'all' or meta['Selected charged rows'] != 'both' or
                meta['Vacuum status'] != 'converged'):
            raise ValueError('Wrong Potts selection/reference')
        close(finite_number(meta, 'Velocity'), VELOCITY)
        close(finite_number(meta, 'Bulk energy density'), BULK)
        vacuum = finite_number(meta, 'Vacuum energy')
        tolerance = finite_number(meta, 'Residual tolerance')
        rows = tables['levels']
        labels = [(0, 0)] + [(q, k) for k in range(n) for q in (1, -1)]
        if len(rows) != len(labels): raise ValueError('Wrong selected-family count')
        energies = {}
        for row, (q, k) in zip(rows, labels):
            if (integer(row, 'charge'), integer(row, 'k')) != (q, k):
                raise ValueError('Wrong charge/momentum or duplicate endpoint')
            if row['branch'] != ('ground' if q == 0 else 'charged-one-hole'):
                raise ValueError('Wrong Potts branch')
            for key in POTTS_COLUMNS[3:-1]: row[key] = finite_number(row, key)
            close(row['momentum'], 2*math.pi*k/n)
            close(row['gap'], row['energy']-vacuum)
            close(row['x_scaled'], n*row['gap']/(2*math.pi*VELOCITY))
            if tolerance <= 0 or not 0 <= row['residual'] <= tolerance or row['iterations'] < 0:
                raise ValueError('Invalid solver diagnostics')
            if row['gap'] < 0: raise ValueError('Level below supplied vacuum')
            energies[q, k] = row['energy']
        close(energies[0, 0], vacuum)
        for k in range(n):
            close(energies[1, k], energies[-1, k])
            close(energies[1, k], energies[1, (-k) % n])
        return meta, tables
    close(finite_number(meta, 'Delta'), 0 if mode == 'free' else .6)
    close(finite_number(meta, 'Imaginary end-field coefficient'), .5 if mode == 'free' else .4)
    if mode == 'free':
        down = int(n/2-parameter)
        if integer(meta, 'Down spins') != down: raise ValueError('Wrong signed sector')
        close(finite_number(meta, 'Sz'), parameter)
        dispersive = [k for k in range(1, n) if n % 2 or k != n//2]
        labels = [(z, modes) for z in range(2 if n % 2 else 3)
                  if 0 <= down-z <= len(dispersive)
                  for modes in itertools.combinations(dispersive, down-z)]
        rows = tables['blocks']
        if len(rows) != len(labels): raise ValueError('Missing/extra Jordan blocks')
        dimension, defective = 0, 0
        for i, (row, (z, modes)) in enumerate(zip(rows, labels)):
            actual_modes = tuple(map(int, row['modes'].split(','))) if row['modes'] else ()
            size = 2 if n % 2 == 0 and z == 1 else 1
            if (integer(row, 'block_id'), integer(row, 'zero_occupation'),
                    integer(row, 'block_size'), actual_modes) != (i, z, size, modes):
                raise ValueError('Wrong block identity, occupation or algebraic size')
            row['energy'] = finite_number(row, 'energy')
            close(row['energy'], math.fsum(math.cos(math.pi*k/n) for k in modes))
            row['block_size'] = size
            dimension += size; defective += size == 2
        if (dimension != math.comb(n, down) or integer(meta, 'Sector dimension') != dimension or
                integer(meta, 'Jordan blocks') != len(rows) or integer(meta, 'Size-two blocks') != defective):
            raise ValueError('Incomplete fixed-Sz space')
        return meta, tables
    # The finite positive-real scan is not the full fixed-Sz Hilbert space.
    ell, delta = parameter, .6
    m, gamma = (n-ell)//2, math.acos(delta)
    threshold = n-m+1-(n-2*m+2)*gamma/math.pi
    slots = math.ceil(threshold)-1
    combos = set(itertools.combinations(range(1, slots+1), m))
    close(finite_number(meta, 'Exclusive label threshold'), threshold)
    close(finite_number(meta, 'Sz'), ell/2)
    for key, value in [('ell=N-2M', ell), ('Bethe roots', m), ('Label slots', slots),
                       ('Candidates', len(combos)), ('Converged candidates', len(combos)),
                       ('Returned levels', len(combos)), ('Sea reference converged', 1)]:
        if integer(meta, key) != value: raise ValueError('Incomplete or wrong scan')
    if meta['Requested levels'] != 'all': raise ValueError('Not a complete regular scan')
    if len(tables['reference']) != 1 or len(tables['levels']) != len(combos):
        raise ValueError('Missing scan/reference rows')
    tolerance = finite_number(meta, 'Residual tolerance')
    reference = finite_number(tables['reference'][0], 'energy_shift')
    offset, seen = 0, set()
    for source in ('levels', 'reference'):
        for i, row in enumerate(tables[source]):
            labels = tuple(map(int, row['numbers'].split(',')))
            if source == 'levels':
                if labels not in combos or labels in seen: raise ValueError('Repeated/unsupported labels')
                seen.add(labels)
            elif labels != tuple(range(1, m+1)): raise ValueError('Wrong sea reference')
            if integer(row, 'state_id') != i or row['converged'] != 'true':
                raise ValueError('Wrong state join/convergence')
            for key in ('energy', 'energy_shift', 'gap_from_sea', 'residual', 'iterations'):
                row[key] = finite_number(row, key)
            if tolerance <= 0 or not 0 <= row['residual'] <= tolerance or row['iterations'] < 0:
                raise ValueError('Unconverged diagnostics')
            close(row['energy'], row['energy_shift']+(n-1)*delta/4)
            close(row['gap_from_sea'], row['energy_shift']-reference)
            roots = tables['roots'][offset:offset+m]; offset += m
            if len(roots) != m: raise ValueError('Missing joined roots')
            z, lam = [], []
            for j, (root, label) in enumerate(zip(roots, labels)):
                if (root['source'], integer(root, 'state_id'), integer(root, 'root_id'),
                        integer(root, 'I'), root['converged']) != (source, i, j, label, 'true'):
                    raise ValueError('Wrong root join/label')
                z.append(finite_number(root, 'z')); lam.append(finite_number(root, 'lambda'))
                if not 0 < z[-1] < math.sqrt((1+delta)/(1-delta)):
                    raise ValueError('Not a finite positive root')
                close(z[-1], math.tanh(lam[-1])/math.tan(gamma/2))
            if any(a >= b for a, b in zip(z, z[1:])): raise ValueError('Unordered roots')
            close(row['energy_shift'], -sum(((1+delta)-(1-delta)*v*v)/(1+v*v) for v in z))
            for j, value in enumerate(z):
                phases = sum(math.atan(delta*(value-v)/(1+delta-(1-delta)*value*v)) +
                             math.atan(delta*(value+v)/(1+delta+(1-delta)*value*v))
                             for k, v in enumerate(z) if k != j)
                close((4*n*math.atan(value)-2*math.pi*labels[j]-2*phases)/(2*n), 0,
                      max(8*tolerance, 1e-13))
            for j, value in enumerate(lam):
                def scatter(x, g): return cmath.sinh(x+1j*g)/cmath.sinh(x-1j*g)
                rhs = math.prod(scatter(value-v, gamma)*scatter(value+v, gamma)
                                for k, v in enumerate(lam) if k != j)
                close(scatter(value, gamma/2)**(2*n), rhs)
    if offset != len(tables['roots']): raise ValueError('Extra/unjoined roots')
    shifts = [row['energy_shift'] for row in tables['levels']]
    if shifts != sorted(shifts): raise ValueError('Incorrect energy ranking')
    return meta, tables


def load_cases(solver_dir=None):
    cases, exports = {}, {}
    for name, case in CASES.items():
        mode, n, parameter = case
        if solver_dir:
            program = 'potts-pbc' if mode == 'potts' else 'xxz-qg-obc'
            command = [Path(solver_dir)/f'bethe-{program}', str(n)]
            if mode == 'free': command += ['--delta', '0', '--sz', str(parameter)]
            if mode == 'scan': command += ['--delta', '.6', '--through-lines', str(parameter), '--excitations', 'all']
            texts = capture_csv_tables(command, schema(case))
        else: texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in schema(case)}
        cases[name] = read_exports(texts, case); exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items(): (DATA/f'{name}-{table}.csv').write_text(text)
    return cases


def casimir(case):
    meta, tables = case
    n = int(meta['Sites'])
    if 'blocks' in tables:
        return -24*n*(min(row['energy'] for row in tables['blocks'])+n/math.pi-.5)/math.pi
    return -6*n*(float(meta['Vacuum energy'])-n*BULK)/(math.pi*VELOCITY)


def oracle_check(cases):
    """Optional small spin/clock-matrix checks; requires NumPy, not Bethe roots."""
    import numpy as np
    from reference_potts import clock_block
    from reference_xxz_nonhermitian import hamiltonian
    for n in (2, 4):
        for row in cases[f'potts-l{n}'][1]['levels']:
            q, k = int(row['charge']), int(row['k'])
            size, entries = clock_block(n, q, k, float, lambda r: np.exp(2j*np.pi*k*r/n))
            h = np.zeros((size, size), complex)
            for ij, value in entries.items(): h[ij] = value
            close(row['energy'], np.linalg.eigvalsh(h)[0])
    for name in ('qg-scan', 'qg-one'):
        _, n, ell = CASES[name]
        ev = np.linalg.eigvals(hamiltonian(n, (n-ell)//2, .6))
        for row in cases[name][1]['levels']: close(min(abs(ev-row['energy'])), 0)
    # Rank/nullity, not eigenvalue splitting, tests defective spectra.
    for n in (2, 4, 6, 7):
        rows = cases[f'qg-free-n{n}'][1]['blocks']
        h = hamiltonian(n, n//2, 0)
        groups = []
        for row in rows:
            group = next((g for g in groups if abs(g[0]-row['energy']) < 1e-10), None)
            if group is None: groups.append([row['energy'], [row['block_size']]])
            else: group[1].append(row['block_size'])
        for energy, sizes in groups:
            a = h-energy*np.eye(len(h)); power = np.eye(len(h))
            for k in (1, 2, 3):
                power = power@a
                nullity = len(h)-np.linalg.matrix_rank(power, tol=1e-9)
                if nullity != sum(min(k, s) for s in sizes): raise ValueError('Wrong Jordan nullities')
    print('Independent clock/spin matrix energies and Jordan nullities passed.')


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-potts-qg-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    for n in (8, 16, 64):
        rows = [r for r in cases[f'potts-l{n}'][1]['levels'] if r['charge'] == '1']
        axes[0].plot([r['momentum']/math.pi for r in rows], [r['gap'] for r in rows], 'o-', markersize=3, label=f'L={n}, q=+1')
    p = [2*math.pi*j/400 for j in range(401)]
    axes[0].plot([v/math.pi for v in p], [3*math.sqrt(3)*math.sin(v/2) for v in p], '--', color='black', label='Thermodynamic branch')
    axes[0].set(xlabel=r'$p/\pi$ (one clock site)', ylabel=r'$E-E_0(L)$', title='Selected charged branch')
    ns = POTTS_SIZES[1:]
    axes[1].plot([1/n**2 for n in ns], [casimir(cases[f'potts-l{n}']) for n in ns], 'o-', label=r'$c_L$')
    axes[1].axhline(.8, color='0.5', linestyle='--', label=r'$c=4/5$')
    axes[1].set(xlabel=r'$1/L^2$', ylabel='Vacuum Casimir estimator', title='Finite-size corrections are retained')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'potts-levels.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    rows = cases['qg-free-n6'][1]['blocks']
    for size, marker, color in [(1, 'o', '#0072B2'), (2, 's', '#D55E00')]:
        chosen = [r for r in rows if r['block_size'] == size]
        axes[0].scatter([int(r['block_id']) for r in chosen], [r['energy'] for r in chosen], marker=marker, color=color, label=f'Size-{size} Jordan block')
    axes[0].set(xlabel='Block ID (not energy rank)', ylabel='Hamiltonian energy', title=r'$N=6$, $S^z=0$: 14 blocks, dimension 20')
    for ns, label, color, limit in [(EVEN, 'Even N, Sz=0', '#D55E00', -2), (ODD, 'Odd N, Sz=1/2', '#0072B2', 1)]:
        axes[1].plot([1/n**2 for n in ns], [casimir(cases[f'qg-free-n{n}']) for n in ns], 'o-', color=color, label=label)
        axes[1].axhline(limit, linestyle='--', color=color, alpha=.6)
    axes[1].set(xlabel=r'$1/N^2$', ylabel=r'$C_N=-24N(E_{\min}+N/\pi-1/2)/\pi$', title='Different parity sectors, different limits')
    for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'qg-blocks.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--oracle', action='store_true', help='Also check small independent matrices (NumPy)')
    args = parser.parse_args()
    cases = load_cases(args.solver_dir)
    if args.oracle: oracle_check(cases)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} Potts/QG calculations.')


if __name__ == '__main__': main()
