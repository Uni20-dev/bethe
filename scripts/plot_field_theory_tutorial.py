#!/usr/bin/env python3
"""Sine-Gordon and Lee-Yang tutorials from native particle and finite-volume tables."""

import argparse
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
SG_LENGTHS = (.1, .2, .5, 1, 2, 5, 10)
LY_LENGTHS = (.001, .003, .01, .03, .1, .3, 1, 2, 5, 10, 20)
LY_EXCITED_LENGTHS = (5, 6, 8, 10, 15, 20)
# kind, mass, circumference (or p for particles), p (unused for Lee-Yang)
CASES = {f'sg-particles-p{p}': ('particles', 1, p, p) for p in (.4, .5, 1)}
CASES['sg-particles-scaled'] = ('particles', 2, .4, .4)
for length in SG_LENGTHS:
    CASES[f'sg-nlie-l{length}'] = ('sg', 1, length, 2)
    CASES[f'sg-by-l{length}'] = ('by', 1, length, 2)
CASES.update({'sg-free-nlie': ('sg', 1, 1, 1), 'sg-free-by': ('by', 1, 1, 1),
              'sg-scaled-nlie': ('sg', 2, .5, 2)})
CASES.update({f'ly-vacuum-l{length}': ('vacuum', 1, length, None) for length in LY_LENGTHS})
CASES['ly-vacuum-scaled'] = ('vacuum', 2, .5, None)
CASES.update({f'ly-excited-l{length}': ('ly', 1, length, None) for length in LY_EXCITED_LENGTHS})
CASES['ly-excited-scaled'] = ('ly', 2, 2.5, None)

LY_DIAGNOSTICS = 'nonlinear_residual nonlinear_error mesh_error cutoff_error direct_tail_bound cutoff intervals iterations cutoffs kernel_products'.split()
SG_DIAGNOSTICS = 'nonlinear_residual kernel_error mesh_error cutoff_error contour_error cutoff intervals iterations kernel_evaluations cutoffs'.split()
END = ['converged', 'status']
SCHEMAS = {
    'particles': {'dispersion': 'branch kind charge rest_energy momentum energy status'.split()},
    'by': {'levels': 'number charge rapidity energy residual quadrature_error tail_bound iterations evaluations phase_status status'.split()},
    'vacuum': {'vacuum': ['casimir_energy', 'scaling_function', 'effective_central_charge'] + LY_DIAGNOSTICS + END},
    'sg': {'levels': ['state', 'casimir_energy', 'scaling_function'] + SG_DIAGNOSTICS + END,
           'source': 'number charge rapidity quantization_residual hole_error source_quadrature_error source_tail_bound root_iterations source_evaluations converged status'.split(),
           'gap': ['gap', 'scaled_gap'] + END},
    'ly': {'levels': ['state', 'casimir_energy', 'scaling_function'] + LY_DIAGNOSTICS + END,
           'source': 'beta pole_displacement quantization_residual source_error root_iterations converged status'.split(),
           'gap': ['gap', 'scaled_gap', 'gap_error'] + END}}
PROGRAMS = {'particles': 'sine-gordon-dispersion', 'by': 'sine-gordon-bethe-yang',
            'sg': 'sine-gordon-excited', 'vacuum': 'lee-yang-vacuum', 'ly': 'lee-yang-excited'}


def close(actual, expected, tolerance=2e-12):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(f'Inconsistent observable: {actual} != {expected}')


def read_exports(texts, case):
    kind, mass, length, p = case
    expected = {'Program': 'bethe-' + PROGRAMS[kind], 'Precision': 'fp64'}
    tables, metas = {}, {}
    for table, columns in SCHEMAS[kind].items():
        meta, raw = read_csv_export(texts[table], columns, expected)
        metas[table] = meta
        parsed = []
        for row in raw:
            if 'converged' in row and row['converged'] != 'true':
                raise ValueError('Unconverged observable')
            numbers = {k: finite_number(row, k) for k in columns
                       if k not in ('state', 'branch', 'kind', 'status', 'converged', 'phase_status')}
            for key, value in numbers.items():
                if (key.endswith(('error', 'residual', 'bound')) or key in
                        ('iterations', 'intervals', 'cutoffs', 'kernel_products', 'kernel_evaluations',
                         'root_iterations', 'source_evaluations', 'evaluations')) and value < 0:
                    raise ValueError('Negative error/work diagnostic')
            parsed.append({**row, **numbers})
        tables[table] = parsed
    meta = next(iter(metas.values()))
    if any(m != meta for m in metas.values()):
        raise ValueError('Tables are not from the same run')
    mass_key = 'Particle mass' if kind in ('ly', 'vacuum') else 'Soliton mass' if kind == 'sg' else 'Soliton mass M'
    close(finite_number(meta, mass_key), mass)
    if p is not None:
        close(finite_number(meta, 'Coupling p'), p)
    if kind == 'particles':
        if meta.get('Energy reference') != 'excitation energy above the infinite-volume vacuum':
            raise ValueError('Wrong dispersion energy reference')
        species = {'soliton': (mass, 1), 'antisoliton': (mass, -1)}
        species.update({f'B{n}': (2*mass*math.sin(n*math.pi*p/2), 0)
                        for n in range(1, 10) if n*p < 1})
        branches = {name: ('particle', m, q) for name, (m, q) in species.items()}
        items = list(species.items())
        for i, (a, (ma, qa)) in enumerate(items):
            for b, (mb, qb) in items[i:]:
                branches[a+'+'+b] = ('threshold', ma+mb, qa+qb)
        if meta.get('Points per branch') != '65' or int(meta['Stable species']) != len(species):
            raise ValueError('Wrong particle enumeration')
        close(finite_number(meta, 'Maximum momentum'), 5*mass)
        grouped = {name: [] for name in branches}
        for row in tables['dispersion']:
            if row['branch'] not in branches:
                raise ValueError('Unexpected particle/threshold')
            k, rest, charge = branches[row['branch']]
            if row['kind'] != k or row['charge'] != charge:
                raise ValueError('Wrong particle kind/charge')
            close(row['rest_energy'], rest)
            close(row['energy'], math.hypot(rest, row['momentum']))
            grouped[row['branch']].append(row)
        for rows in grouped.values():
            if len(rows) != 65:
                raise ValueError('Incomplete dispersion branch')
            for i, row in enumerate(rows):
                close(row['momentum'], 5*mass*i/64)
        tables['branches'] = grouped
        return meta, tables
    close(finite_number(meta, 'Circumference L' if kind == 'by' else 'Circumference'), length)
    if kind == 'by':
        if len(tables['levels']) != 1 or 'wrapping corrections omitted' not in meta.get('Approximation', ''):
            raise ValueError('Wrong asymptotic scope')
        row = tables['levels'][0]
        if row['number'] != .5 or row['charge'] != 2 or row['phase_status'] != 'converged':
            raise ValueError('Wrong Bethe-Yang sector/status')
        close(row['energy'], 2*mass*math.cosh(row['rapidity']))
        if sum(row[k] for k in ('residual', 'quadrature_error', 'tail_bound')) > finite_number(meta, 'Counting tolerance'):
            raise ValueError('Unresolved Bethe-Yang equation')
        return meta, tables
    tolerance = finite_number(meta, 'Absolute Y tolerance' if kind == 'vacuum' else 'Absolute Y tolerance per state')
    if tolerance <= 0:
        raise ValueError('Invalid tolerance')
    levels = tables['vacuum' if kind == 'vacuum' else 'levels']
    if kind == 'vacuum':
        if len(levels) != 1 or meta.get('CFT central charge (theory)') != '-22/5':
            raise ValueError('Wrong vacuum scope')
        close(levels[0]['effective_central_charge'], -6*levels[0]['scaling_function']/math.pi)
    else:
        if kind == 'sg' and meta.get('Calculation') != 'exact continuum two-hole NLIE; not Bethe-Yang':
            raise ValueError('Wrong finite-volume calculation')
        names = ['vacuum', 'two_soliton' if kind == 'sg' else 'one_particle']
        if [r['state'] for r in levels] != names or len(tables['source']) != 1 or len(tables['gap']) != 1:
            raise ValueError('Missing or mislabelled level/source/gap')
        gap, source = tables['gap'][0], tables['source'][0]
        close(gap['gap'], levels[1]['casimir_energy']-levels[0]['casimir_energy'])
        close(gap['scaled_gap'], length*gap['gap']/(2*math.pi if kind == 'sg' else 1))
        if gap['gap'] <= 0:
            raise ValueError('Nonpositive selected excitation gap')
        if kind == 'sg':
            if source['number'] != .5 or source['charge'] != 2 or source['rapidity'] <= 0:
                raise ValueError('Wrong two-soliton sector')
            if any(source[k] > tolerance for k in ('quantization_residual', 'hole_error', 'source_quadrature_error', 'source_tail_bound')):
                raise ValueError('Unresolved NLIE source')
        else:
            close(source['beta']-source['pole_displacement'], math.pi/6)
            if (not 0 < source['pole_displacement'] < math.pi/6
                    or max(source['source_error'], source['quantization_residual']) > tolerance):
                raise ValueError('Unresolved Lee-Yang source')
            if gap['gap_error'] < sum(r['cutoff_error'] for r in levels)/length:
                raise ValueError('Gap error omits level errors')
    for row in levels:
        close(row['scaling_function'], length*row['casimir_energy'])
        for key in ('mesh_error', 'cutoff_error', 'contour_error', 'kernel_error', 'nonlinear_error', 'direct_tail_bound'):
            if key in row and row[key] > tolerance:
                raise ValueError('Unresolved level refinement')
        if row['cutoffs'] < 2 or row['intervals'] <= 0 or row['cutoff'] <= 0:
            raise ValueError('Missing mesh/cutoff verification')
    return meta, tables


def load_cases(solver_dir=None):
    exports, cases = {}, {}
    for name, case in CASES.items():
        kind, mass, length, p = case
        if solver_dir:
            command = [Path(solver_dir)/('bethe-'+PROGRAMS[kind]), '--mass', str(mass)]
            if p is not None:
                command += ['--p', str(p)]
            command += (['--points', '65', '--max-momentum', str(5*mass)] if kind == 'particles'
                        else ['--length', str(length)])
            if kind == 'sg':
                command += ['--tolerance', '1e-7']
            texts = capture_csv_tables(command, SCHEMAS[kind])
        else:
            texts = {table: (DATA/f'{name}-{table}.csv').read_text() for table in SCHEMAS[kind]}
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
    plt.rcParams.update({'svg.hashsalt': 'bethe-field-theory-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(7, 4.3), layout='constrained')
    branches = cases['sg-particles-p0.4'][1]['branches']
    for name, label, color, style in [('soliton', r'$s,\bar s$ (charges $\pm1$)', '#0072B2', '-'),
                                     ('B1', r'$B_1$ (charge 0)', '#009E73', '-'),
                                     ('B2', r'$B_2$ (charge 0)', '#D55E00', '-'),
                                     ('soliton+antisoliton', r'$s+\bar s$ threshold (charge 0)', '0.3', '--')]:
        rows = branches[name]
        ax.plot([r['momentum'] for r in rows], [r['energy'] for r in rows], style, color=color, label=label)
    ax.set(xlabel=r'$k/M$', ylabel=r'$E/M$', title='Sine-Gordon: particles and a neutral threshold, p=0.4', xlim=(0, 5))
    ax.legend(fontsize=9); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'sg-particles.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), layout='constrained')
    exact = [cases[f'sg-nlie-l{l}'][1]['gap'][0] for l in SG_LENGTHS]
    by = [cases[f'sg-by-l{l}'][1]['levels'][0]['energy'] for l in SG_LENGTHS]
    axes[0].plot(SG_LENGTHS, [r['scaled_gap'] for r in exact], 'o-', color='#0072B2', label='NLIE gap')
    axes[0].plot(SG_LENGTHS, [l*e/(2*math.pi) for l, e in zip(SG_LENGTHS, by)], 's--', color='#D55E00', label='Bethe-Yang')
    axes[0].axhline(.75, color='0.5', linestyle=':', label='UV value 3/4')
    axes[0].set(xscale='log', xlabel=r'$ML$', ylabel=r'$L\Delta E/(2\pi)$', title='Same-charge pair: p=2, I=1/2')
    axes[1].plot(SG_LENGTHS, [r['gap']-e for r, e in zip(exact, by)], 'o-', color='#009E73')
    axes[1].axhline(0, color='0.5', linewidth=.8)
    axes[1].set(xscale='log', xlabel=r'$ML$', ylabel=r'$(\Delta E_{\mathrm{NLIE}}-E_{\mathrm{BY}})/M$', title='Finite-volume correction, not root error')
    axes[0].legend(fontsize=9)
    for ax in axes: ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'sg-finite-volume.svg'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
    ceff = [cases[f'ly-vacuum-l{l}'][1]['vacuum'][0]['effective_central_charge'] for l in LY_LENGTHS]
    ax.semilogx(LY_LENGTHS, ceff, 'o-', color='#0072B2', label=r'Computed $-6LE_{0,C}/\pi$')
    ax.axhline(.4, color='0.5', linestyle=':', label=r'UV $c_{\mathrm{eff}}=2/5$, not $c=-22/5$')
    ax.set(xlabel=r'$r=mL$', ylabel=r'$c_{\mathrm{eff}}(r)$', title='Scaling Lee-Yang: periodic vacuum', ylim=(0, .43))
    ax.legend(fontsize=9); ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'ly-vacuum.svg'); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), layout='constrained')
    values = [cases[f'ly-excited-l{l}'][1] for l in LY_EXCITED_LENGTHS]
    axes[0].plot(LY_EXCITED_LENGTHS, [v['gap'][0]['gap'] for v in values], 'o-', label=r'Gap $(E_{1,C}-E_{0,C})/m$', color='#0072B2')
    axes[0].plot(LY_EXCITED_LENGTHS, [v['levels'][1]['casimir_energy'] for v in values], 's--', label=r'Level $E_{1,C}/m$', color='#D55E00')
    axes[0].axhline(1, color='0.5', linestyle=':', label='Infinite-volume mass')
    axes[0].set(xlabel=r'$mL$', ylabel='Energy / mass', title='Level and gap are different')
    axes[0].legend(fontsize=9)
    axes[1].semilogy(LY_EXCITED_LENGTHS, [v['gap'][0]['gap']-1 for v in values], 'o-', color='#009E73', label=r'$\Delta E/m-1$')
    axes[1].semilogy(LY_EXCITED_LENGTHS, [v['gap'][0]['gap_error'] for v in values], 's:', color='0.5', label='Estimated numerical gap error / m')
    axes[1].set(xlabel=r'$mL$', ylabel='Dimensionless difference / estimate', title='Resolved approach to the particle mass')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'ly-gap.svg'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver-dir', type=Path, help='regenerate using the five native executables here')
    parser.add_argument('--check', action='store_true', help='validate exports without plotting dependencies')
    args = parser.parse_args()
    cases = load_cases(args.solver_dir)
    if not args.check:
        plot(cases)
    print(f'Validated {len(cases)} field-theory calculations.')


if __name__ == '__main__':
    main()
