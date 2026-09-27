#!/usr/bin/env python3
"""Both signs of the free-end biquadratic chain, from validated frontend tables."""

import argparse
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
# kind, N, through-lines (q/real) or mode count (pair/triple), ferro sign
CASES = {f'bq-{sign}-q-l{ell}': ('q', 4, ell, sign == 'ferro')
         for sign in ('af', 'ferro') for ell in (0, 2, 4)}
CASES.update({f'bq-{sign}-real-l0': ('real', 4, 0, sign == 'ferro') for sign in ('af', 'ferro')})
CASES.update({f'bq-af-singlet-n{n}': ('singlet', n, 0, False) for n in (4, 8, 16, 32, 64)})
CASES['bq-ferro-band-n32'] = ('band', 32, 30, True)
CASES.update({f'bq-ferro-{kind}-n{n}': (kind, n, 8 if n == 32 else 1, True)
              for kind in ('pair', 'triple') for n in (16, 32, 64, 128)})
BASE = ['state_id', 'through_lines', 'multiplicity', 'energy', 'gap', 'tl_energy', 'reference_energy', 'residual']
END = ['iterations', 'converged', 'status']
SCHEMAS = {'q': BASE + ['bethe_residual', 'bethe_residual_bound'] + END,
           'real': BASE + END,
           'singlet': BASE + ['phase_residual', 'modulus_residual'] + END,
           'band': ['state_id', 'through_lines', 'defects', 'multiplicity', 'mode', 'wave_number', 'energy', 'gap']}
SCHEMAS['pair'] = SCHEMAS['triple'] = ['state_id', 'mode'] + SCHEMAS['singlet'][1:]


def multiplicity(ell):
    a, b = 1, 3
    for _ in range(ell):
        a, b = b, 3*b-a
    return a


def read_exports(texts, case):
    kind, n, selection, ferro = case
    ell = n-4 if kind == 'pair' else n-6 if kind == 'triple' else selection
    count = ({0: 2, 2: 3, 4: 1}[ell] if kind == 'q' else n-1 if kind == 'band'
             else selection if kind in ('pair', 'triple') else 1)
    expected = {'Program': 'bethe-biquadratic-obc', 'Sites': str(n), 'Precision': 'fp64', 'Spin': '1',
                'TL loop weight': '3', 'XXZ Delta': '1.5',
                'Hamiltonian': f'H={"+" if ferro else "-"}sum_i (S_i.S_(i+1))^2',
                'Gap reference': 'E-(N-1); exact degenerate ferro ground space' if ferro else 'E-E0; global singlet ground state'}
    if kind == 'q':
        expected.update({'Calculation': 'Q-system spectrum search (real and complex roots)',
                         'Expected module dimension': str(count), 'Discovered levels': str(count),
                         'Ordering': 'numerically complete TL module; count matched, not a rigorous certificate'})
    elif kind == 'real':
        expected.update({'Calculation': 'restricted real-root excitations',
                         'Ordering': 'complete within supported family', 'Ground status': 'converged'})
    elif kind == 'singlet':
        expected.update({'Calculation': 'selected complex-root singlet excitation',
                         'Ordering': 'targeted branch, not an exhaustive search or a global first-excitation guarantee'})
    elif kind == 'band':
        expected.update({'Status': 'exact spectral rules', 'Calculation': 'complete one-defect TL module',
                         'Coverage': 'all N-1 eigenvalues of ell=N-2; NOT the full excited spectrum',
                         'Wave number': 'k=pi*j/N is an OBC standing-wave coordinate, not lattice momentum'})
    else:
        word = 'two' if kind == 'pair' else 'three'
        expected.update({'Status': 'converged targeted branches',
                         'Calculation': f'targeted {word}-defect bound-{kind} modes', 'Selected modes': str(count),
                         'Coverage': f'selected {word}-string family only; NOT the full {word}-defect or excited spectrum',
                         'Wave number': 'mode and string center are branch coordinates, not lattice momentum'})
    meta, raw = read_csv_export(texts['states'], SCHEMAS[kind], expected, row_status=None)
    tolerance = finite_number(meta, 'Residual tolerance')
    if len(raw) != count or tolerance <= 0:
        raise ValueError('Wrong number of selected levels or tolerance')
    rows = []
    for i, row in enumerate(raw):
        if row['state_id'] != str(i) or row['through_lines'] != str(ell):
            raise ValueError('Wrong ID or TL sector')
        dimension = multiplicity(ell)
        if row['multiplicity'] != (str(dimension) if dimension < 2**64 else ''):
            raise ValueError('Invalid multiplicity/overflow encoding')
        numbers = {key: finite_number(row, key) for key in ('energy', 'gap')}
        if kind != 'band':
            status = 'converged; numerical admissibility checks passed' if kind == 'q' else 'converged'
            if row['converged'] != 'true' or row['status'] != status:
                raise ValueError('Unverified level')
            for key in SCHEMAS[kind]:
                if key not in ('state_id', 'through_lines', 'multiplicity', 'converged', 'status'):
                    numbers[key] = finite_number(row, key)
            if not 0 <= numbers['residual'] <= tolerance or numbers['iterations'] < 0:
                raise ValueError('Invalid residual/work')
            if any(not 0 <= numbers[key] <= tolerance for key in ('phase_residual', 'modulus_residual') if key in numbers):
                raise ValueError('Unconverged string equation')
            physical = (-2 if ferro else 2)*numbers['reference_energy'] + (1 if ferro else -1)*7*(n-1)/4
            if not math.isclose(numbers['energy'], physical, abs_tol=3e-12):
                raise ValueError('Wrong physical/reference energy map')
            tl = numbers['energy']-(n-1) if ferro else numbers['energy']+(n-1)
            if not math.isclose(numbers['tl_energy'], tl, abs_tol=3e-12):
                raise ValueError('Wrong TL energy shift')
            if kind == 'q' and not 0 <= numbers['bethe_residual'] <= numbers['bethe_residual_bound']:
                raise ValueError('Failed original-equation admissibility')
        else:
            numbers.update(mode=int(row['mode']), wave_number=finite_number(row, 'wave_number'))
            if row['defects'] != '1' or numbers['mode'] != n-1-i:
                raise ValueError('Wrong one-defect mode order')
        if kind in ('pair', 'triple') and row['mode'] != str(i+1):
            raise ValueError('Wrong droplet mode order')
        if ferro and not math.isclose(numbers['gap'], numbers['energy']-(n-1), abs_tol=3e-12):
            raise ValueError('Wrong ferro ground reference')
        rows.append(dict(numbers, through_lines=ell, multiplicity=None if dimension >= 2**64 else dimension))
    if kind == 'singlet':
        columns = ['state_id', 'energy', 'residual', 'converged', 'status']
        ref_meta, ref = read_csv_export(texts['reference'], columns, expected)
        if ref_meta['Command'] != meta['Command'] or len(ref) != 1 or ref[0]['state_id'] != '1' or ref[0]['converged'] != 'true':
            raise ValueError('Invalid singlet reference table')
        ground = finite_number(ref[0], 'energy')
        if not 0 <= finite_number(ref[0], 'residual') <= tolerance:
            raise ValueError('Unconverged ground reference')
        if not math.isclose(rows[0]['gap'], rows[0]['energy']-ground, abs_tol=3e-12):
            raise ValueError('Wrong singlet gap')
    if n == 4 and not ferro:
        ground = -(15+math.sqrt(17))/2
        if any(not math.isclose(r['gap'], r['energy']-ground, abs_tol=3e-12) for r in rows):
            raise ValueError('Wrong four-site AF ground reference')
    content = {}
    if kind == 'q':
        columns = ['state_id', 'through_lines', 'spin', 'multiplets', 'magnetic_states', 'status']
        spin_meta, spins = read_csv_export(texts['spin_content'], columns, expected, row_status=None)
        if spin_meta['Command'] != meta['Command']:
            raise ValueError('Spin table belongs to another calculation')
        content = {i: {} for i in range(count+1)}
        for row in spins:
            i, e, spin, mult, magnetic = (int(row[k]) for k in columns[:-1])
            wanted_ell = ell if i < count else 4 if ferro else 0
            if (i not in content or row['status'] != 'exact' or e != wanted_ell or not 0 <= spin <= e
                    or mult <= 0 or magnetic != (2*spin+1)*mult or spin in content[i]):
                raise ValueError('Invalid physical-spin row')
            content[i][spin] = mult
        for i, counts in content.items():
            e = ell if i < count else 4 if ferro else 0
            if sum((2*s+1)*m for s, m in counts.items()) != multiplicity(e):
                raise ValueError('Incomplete spin decomposition')
    return meta, rows, content


def load_cases(solver=None):
    exports, cases = {}, {}
    for name, case in CASES.items():
        kind, n, selection, ferro = case
        args = [solver, n] + (['--ferromagnetic'] if ferro else [])
        args += (['--q-spectrum', '--through-lines', selection] if kind == 'q' else
                 ['--excitations', 'all', '--through-lines', selection] if kind == 'real' else
                 ['--singlet-excitation'] if kind == 'singlet' else ['--one-defect'] if kind == 'band' else
                 ['--bound-pairs' if kind == 'pair' else '--bound-triples', selection])
        tables = ('states', 'spin_content') if kind == 'q' else ('states', 'reference') if kind == 'singlet' else ('states',)
        texts = (capture_csv_tables(args, tables) if solver else
                 {table: (DATA / f'{name}-{table}.csv').read_text() for table in tables})
        cases[name] = read_exports(texts, case)
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
    plt.rcParams.update({'svg.hashsalt': 'bethe-biquadratic-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    for ax, sign, title in zip(axes, ('af', 'ferro'), ('Negative sign (default)', 'Positive sign (--ferromagnetic)')):
        for ell in (0, 2, 4):
            rows = cases[f'bq-{sign}-q-l{ell}'][1]
            real_energy = cases[f'bq-{sign}-real-l0'][1][0]['energy']
            for row in rows:
                missing = ell == 0 and abs(row['energy']-real_energy) > 1e-10
                color = '#D55E00' if missing else '#0072B2'
                ax.hlines(row['gap'], ell-.35, ell+.35, color=color, linewidth=2)
                ax.text(ell+.42, row['gap'], f"×{row['multiplicity']}", va='center', fontsize=9)
        ax.set(title=title, xlabel=r'TL through-lines $\ell$ (not spin)', ylabel=r'$E-E_0$',
               xticks=[0, 2, 4], xlim=(-.65, 5.1), ylim=(-.4, 7.1))
        ax.grid(axis='y', alpha=.15)
    fig.suptitle('N=4: orange = singlet missed by the real-root family', fontsize=12)
    save_svg(fig, FIGURES / 'bq-four-site.svg'); plt.close(fig)

    sizes = (4, 8, 16, 32, 64)
    fig, ax = plt.subplots(figsize=(7, 3.7), layout='constrained')
    ax.plot(sizes, [cases[f'bq-af-singlet-n{n}'][1][0]['gap'] for n in sizes], 'o-', color='#0072B2')
    ax.set(title='Negative sign: a targeted complex-root singlet', xlabel='Open-chain length N', ylabel=r'$E-E_0$')
    ax.set_xscale('log', base=2); ax.set_xticks(sizes, [str(n) for n in sizes]); ax.grid(alpha=.15)
    save_svg(fig, FIGURES / 'bq-af-singlet.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), layout='constrained')
    rows = cases['bq-ferro-band-n32'][1]
    axes[0].plot([r['wave_number']/math.pi for r in rows], [r['gap'] for r in rows], 'o-', markersize=3, color='#0072B2')
    axes[0].set(title='Complete one-defect module', xlabel=r'Standing-wave $k/\pi$', ylabel=r'$E-(N-1)$', xlim=(0, 1))
    for kind, label, color in [('pair', 'Bound pair', '#D55E00'), ('triple', 'Bound triple', '#009E73')]:
        rows = cases[f'bq-ferro-{kind}-n32'][1]
        axes[1].plot([r['mode'] for r in rows], [r['gap'] for r in rows], 'o-', color=color, markersize=4, label=label)
    axes[1].set(title='Selected droplet modes', xlabel='Mode number (not momentum)', ylabel=r'$E-(N-1)$', xticks=[1, 2, 4, 6, 8])
    axes[1].legend(fontsize=9)
    for ax in axes: ax.grid(alpha=.15)
    fig.suptitle('Positive sign, N=32: different TL modules', fontsize=12)
    save_svg(fig, FIGURES / 'bq-ferro-branches.svg'); plt.close(fig)

    sizes = (16, 32, 64, 128)
    fig, ax = plt.subplots(figsize=(7, 3.8), layout='constrained')
    for kind, limit, label, color in [('single', 1, r'One defect, $g_1=1$', '#0072B2'),
                                    ('pair', 5/3, r'Pair, $g_2=5/3$', '#D55E00'),
                                    ('triple', 2, r'Triple, $g_3=2$', '#009E73')]:
        gaps = [(float(cases[f'bq-ferro-pair-n{n}'][0]['Exact gap above ground space']) if kind == 'single'
                 else cases[f'bq-ferro-{kind}-n{n}'][1][0]['gap']) for n in sizes]
        ax.plot(sizes, [gap-limit for gap in gaps], 'o-', color=color, label=label)
    ax.set_xscale('log', base=2); ax.set_yscale('log')
    ax.set_xticks(sizes, [str(n) for n in sizes])
    ax.set(title='Lowest selected modes approach different bulk edges', xlabel='Open-chain length N', ylabel=r'Gap minus $g_M$')
    ax.legend(fontsize=9); ax.grid(alpha=.15)
    save_svg(fig, FIGURES / 'bq-ferro-size.svg'); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} signed biquadratic calculations.')
