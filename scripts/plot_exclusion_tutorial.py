#!/usr/bin/env python3
"""Periodic TASEP/ASEP relaxation tutorials from native, root-checked exports."""

import argparse
import cmath
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
SIZES = (8, 16, 32, 64, 128)
BIASES = (1, .75, .5, .25, .1, .01, 0)
# model, sites, particles, right rate, left rate
CASES = {}
for size in SIZES:
    for family, particles in [('half', size//2), ('quarter', size//4), ('single', 1)]:
        CASES[f'tasep-{family}-l{size}'] = ('tasep', size, particles, 1, 0)
CASES.update({
    'tasep-scaled': ('tasep', 16, 4, 2, 0),
    'tasep-holes': ('tasep', 16, 12, 1, 0),
    'tasep-empty': ('tasep', 8, 0, 1, 0),
    'tasep-full': ('tasep', 8, 8, 1, 0),
    'tasep-small': ('tasep', 5, 2, 1, 0)})
CASES.update({f'asep-bias{b}': ('asep', 16, 4, (1+b)/2, (1-b)/2) for b in BIASES})
CASES.update({
    'asep-scaled': ('asep', 16, 4, 1.5, .5),
    'asep-reflected': ('asep', 16, 4, .25, .75),
    'asep-holes': ('asep', 16, 12, .75, .25),
    'asep-single': ('asep', 16, 1, .75, .25),
    'asep-symmetric-half': ('asep', 16, 8, .5, .5),
    'asep-empty': ('asep', 8, 0, .75, .25),
    'asep-full': ('asep', 8, 8, .75, .25),
    'asep-small': ('asep', 5, 2, 1, .5)})
SCHEMA = 'has_mode lambda_real lambda_imag gap frequency residual iterations seed_iterations converged status'.split()


def close(actual, expected, tolerance=3e-11):
    if abs(actual-expected) > tolerance*max(1, abs(expected)):
        raise ValueError(f'Inconsistent eigenvalue/root identity: {actual} != {expected}')


def read_exports(texts, case):
    model, size, particles, right, left = case
    reduced = min(particles, size-particles)
    status = 'converged' if reduced else 'stationary_only'
    expected = {'Program': f'bethe-{model}-pbc', 'Precision': 'fp64', 'Sites': str(size),
                'Particles': str(particles), 'Reduced root count': str(reduced), 'Status': status,
                'Units': 'inverse time; not quantum energies', 'Stationary eigenvalue': '0',
                'Calculation': 'leading relaxation pair; representative Im(lambda)>=0',
                'Has relaxation mode': 'true' if reduced else 'false',
                'Generator': 'dP/dt=M P; '+('right hops only' if model == 'tasep' else 'right/left hops')+'; columns sum to zero'}
    meta, raw = read_csv_export(texts['relaxation'], SCHEMA, expected, row_status=None)
    coords = ['index', 'Z_real', 'Z_imag'] if model == 'tasep' else ['index', 'v_real', 'v_imag']
    root_meta, root_rows = read_csv_export(texts['roots'], coords, expected, row_status=None)
    if root_meta != meta or len(raw) != 1:
        raise ValueError('Missing mode row or mismatched run provenance')
    close(finite_number(meta, 'Right-hop rate'), right)
    if model == 'asep':
        close(finite_number(meta, 'Left-hop rate'), left)
    row = raw[0]
    if row['status'] != status or row['converged'] != 'true' or row['has_mode'] != expected['Has relaxation mode']:
        raise ValueError('Unconverged or mislabelled mode')
    numbers = {k: finite_number(row, k) for k in ('residual', 'iterations', 'seed_iterations')}
    tolerance = finite_number(meta, 'Residual tolerance')
    if tolerance <= 0 or any(v < 0 for v in numbers.values()) or numbers['residual'] > tolerance:
        raise ValueError('Invalid work/residual diagnostics')
    if not reduced:
        if root_rows or any(row[k] != '' for k in ('lambda_real', 'lambda_imag', 'gap', 'frequency')):
            raise ValueError('Stationary-only sector has invented relaxation observables')
        return meta, {**numbers, 'has_mode': False, 'roots': []}
    numbers.update({k: finite_number(row, k) for k in ('lambda_real', 'lambda_imag', 'gap', 'frequency')})
    if numbers['gap'] <= 0 or numbers['frequency'] < 0:
        raise ValueError('Wrong decay/frequency sign')
    close(numbers['gap'], -numbers['lambda_real'])
    close(numbers['frequency'], numbers['lambda_imag'])
    for column, label in [('gap', 'Relaxation gap'), ('frequency', 'Angular frequency'),
                          ('lambda_real', 'Eigenvalue real part'), ('lambda_imag', 'Eigenvalue imaginary part')]:
        close(numbers[column], finite_number(meta, label))
    lam = complex(numbers['lambda_real'], numbers['lambda_imag'])
    analytic = model == 'asep' and (reduced == 1 or right == left)
    if model == 'asep' and meta.get('Analytic result') != ('true' if analytic else 'false'):
        raise ValueError('Wrong analytic/root branch')
    if len(root_rows) != (0 if analytic else reduced):
        raise ValueError('Wrong number of reduced roots')
    roots = []
    for i, root in enumerate(root_rows):
        if root['index'] != str(i):
            raise ValueError('Wrong root index')
        roots.append(complex(finite_number(root, coords[1]), finite_number(root, coords[2])))
    if reduced == 1 or right == left:
        exact = complex(-2*(right+left)*math.sin(math.pi/size)**2,
                        abs(right-left)*math.sin(2*math.pi/size) if reduced == 1 else 0)
        close(lam, exact)
    if model == 'tasep':
        close(lam, right*sum((z-1)/2 for z in roots))
        # Original fugacity equations, evaluated through logarithms to avoid large powers.
        log_y = size*math.log(2)+1j*math.pi+sum(cmath.log((z-1)/(z+1)) for z in roots)
        for z in roots:
            d = reduced*cmath.log(1-z)+(size-reduced)*cmath.log(1+z)-log_y
            close(complex(d.real, math.remainder(d.imag, 2*math.pi))/size, 0)
        wave_roots = [2/(z+1) for z in roots]
    elif not analytic:
        ratio = min(right, left)/max(right, left)
        close(finite_number(meta, 'Last reached rate ratio'), ratio)
        delta = 1-ratio
        index = int(meta['Wave root index'])
        if not 0 <= index < reduced:
            raise ValueError('Invalid wave-root index')
        base = complex(finite_number(meta, 'Wave base real'), finite_number(meta, 'Wave base imaginary'))
        close(abs(1+base), 1)
        wave_roots = [1+delta*v+(base if i == index else 0) for i, v in enumerate(roots)]
        close(lam, max(right, left)*sum(ratio*z+1/z-1-ratio for z in wave_roots))
        for i, zi in enumerate(wave_roots):
            rhs = complex((-1)**(reduced-1))
            for j, zj in enumerate(wave_roots):
                if i != j:
                    rhs *= (ratio*zi*zj-(1+ratio)*zi+1)/(ratio*zi*zj-(1+ratio)*zj+1)
            close(zi**size, rhs, 2e-8)  # Reconstructed z loses conditioning near symmetry.
    else:
        wave_roots = []
    if wave_roots:
        product = math.prod(wave_roots)
        if min(abs(product-cmath.exp(sign*2j*math.pi/size)) for sign in (-1, 1)) > 1e-9:
            raise ValueError('Wrong first-harmonic translation factor')
    return meta, {**numbers, 'has_mode': True, 'roots': roots}


def load_cases(solver_dir=None):
    cases, exports = {}, {}
    for name, case in CASES.items():
        model, size, particles, right, left = case
        if solver_dir:
            args = [Path(solver_dir)/f'bethe-{model}-pbc', str(size), '--particles', str(particles)]
            args += ['--rate', str(right)] if model == 'tasep' else ['--right-rate', str(right), '--left-rate', str(left)]
            texts = capture_csv_tables(args, ('relaxation', 'roots'))
        else:
            texts = {t: (DATA/f'{name}-{t}.csv').read_text() for t in ('relaxation', 'roots')}
        cases[name] = read_exports(texts, case)
        exports[name] = texts
    if solver_dir:
        for name, texts in exports.items():
            for table, text in texts.items():
                (DATA/f'{name}-{table}.csv').write_text(text)
    return cases


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.hashsalt': 'bethe-exclusion-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    families = [('half', r'$N=L/2$', '#0072B2'), ('quarter', r'$N=L/4$', '#D55E00'),
                ('single', r'$N=1$', '#009E73')]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), layout='constrained')
    for family, label, color in families:
        gaps = [cases[f'tasep-{family}-l{l}'][1]['gap'] for l in SIZES]
        axes[0].loglog(SIZES, gaps, 'o-', color=color, label=label)
        slopes = [math.log(a/b, 2) for a, b in zip(gaps, gaps[1:])]
        axes[1].semilogx(SIZES[1:], slopes, 'o-', color=color, label=label, base=2)
    axes[0].set(xlabel='Ring size L', ylabel=r'$g/r$', title='TASEP: decay gaps from finite-size roots')
    axes[1].set(xlabel='L (compared with L/2)', ylabel=r'$z_{\mathrm{eff}}=\log_2[g(L/2)/g(L)]$', title='Fixed density is not a single particle')
    for z in (1.5, 2): axes[1].axhline(z, color='0.5', linestyle=':', linewidth=1)
    for ax in axes: ax.legend(fontsize=9); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'tasep-scaling.svg'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    for family, label, color in families:
        values = [l*cases[f'tasep-{family}-l{l}'][1]['frequency']/(2*math.pi) for l in SIZES]
        ax.semilogx(SIZES, values, 'o-', color=color, label=label, base=2)
    for value in (0, .5, 1): ax.axhline(value, color='0.5', linestyle=':', linewidth=.8)
    ax.set(xlabel='Ring size L', ylabel=r'$L\omega/(2\pi r)$', title='Oscillation is separate from relaxation')
    ax.legend(loc='center right', bbox_to_anchor=(1, .28)); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'tasep-frequency.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    rows = [cases[f'asep-bias{b}'][1] for b in BIASES]
    axes[0].plot(BIASES, [r['gap'] for r in rows], 'o-', color='#0072B2', label='Selected relaxation branch')
    axes[0].axhline(2*math.sin(math.pi/16)**2, color='0.5', linestyle=':', label='Exact symmetric gap')
    axes[1].plot(BIASES, [r['frequency'] for r in rows], 'o-', color='#D55E00')
    axes[0].set(ylabel=r'$g/(r+s)$', title='Decay stays finite on this finite ring')
    axes[1].set(ylabel=r'$\omega/(r+s)$', title='Frequency vanishes at symmetry')
    for ax in axes:
        ax.set(xlabel=r'Bias $(r-s)/(r+s)$', xlim=(0, 1)); ax.grid(alpha=.15)
    axes[0].legend(fontsize=8)
    fig.suptitle('ASEP: L=16, N=4, fixed total rate r+s=1', fontsize=12)
    save_svg(fig, FIGURES/'asep-bias.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path, help='regenerate with the native TASEP and ASEP frontends')
    parser.add_argument('--check', action='store_true', help='validate exports without plotting dependencies')
    args = parser.parse_args()
    cases = load_cases(args.solver_dir)
    if not args.check: plot(cases)
    print(f'Validated {len(cases)} periodic exclusion calculations.')


if __name__ == '__main__':
    main()
