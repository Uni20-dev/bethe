#!/usr/bin/env python3
"""Checked native-solver exports; unbroadened thermodynamic XXX density."""
import argparse
import math
from pathlib import Path
from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg, solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'docs/tutorials/data'
FIGURE = ROOT/'docs/tutorials/figures/xxx-structure-factor-thermo.svg'
NQ, NW = 129, 401
COLUMNS = ['q', 'omega', 'density', 'estimated_error', 'evaluations', 'status']


def read_case(texts):
    expected = {'Program': 'bethe-xxx-structure-factor-thermo', 'Precision': 'fp64',
                'Exchange J': '1', 'Channel': 'zz', 'Spectral points': str(NQ*NW), 'Failed points': '0',
                'Family': 'Thermodynamic two-spinon contribution; partial DSF; no rescaling',
                'Density convention': 'S(q,w) itself; integrals use dw/(2*pi); raising=2*zz'}
    meta, raw = read_csv_export(texts['spectrum'], COLUMNS, expected, row_status=None)
    other, edges = read_csv_export(texts['continuum'], ['q', 'lower', 'upper'], expected, row_status=None)
    if meta != other or len(raw) != NQ*NW or len(edges) != NQ:
        raise ValueError('Inconsistent grid/metadata')
    values, singular = [], 0
    for iq, edge in enumerate(edges):
        q = 2*math.pi*iq/(NQ-1)
        lo, hi = math.pi/2*abs(math.sin(q)), math.pi*math.sin(q/2)
        if any(not math.isclose(finite_number(edge, key), value, abs_tol=2e-14)
               for key, value in [('q', q), ('lower', lo), ('upper', hi)]):
            raise ValueError('Incorrect continuum')
        column = []
        for iw in range(NW):
            row = raw[iq*NW+iw]
            w = math.pi*iw/(NW-1)
            if (not math.isclose(finite_number(row, 'q'), q, abs_tol=2e-14)
                    or not math.isclose(finite_number(row, 'omega'), w, abs_tol=2e-14)):
                raise ValueError('Incorrect grid ordering')
            if row['status'] == 'lower_threshold':
                if (iq in (0, NQ-1) or not math.isclose(w, lo, abs_tol=2e-14)
                        or row['density'] or row['estimated_error']):
                    raise ValueError('Incorrect singular point')
                singular += 1
                column.append(None)
                continue
            if row['status'] not in ('converged', 'outside_continuum'):
                raise ValueError('Numerical failure')
            density = finite_number(row, 'density')
            error = finite_number(row, 'estimated_error')
            if density < 0 or error < 0 or error > 1e-9*max(1, density):
                raise ValueError('Invalid density/error')
            outside = w < lo-1e-14 or w > hi+1e-14 or iq in (0, NQ-1)
            if outside and density != 0:
                raise ValueError('Nonzero density outside continuum')
            if lo+1e-13 < w < hi-1e-13 and density <= 0:
                raise ValueError('Missing interior density')
            if row['status'] == 'outside_continuum' and density != 0:
                raise ValueError('Invalid outside status')
            column.append(density)
        values.append(column)
    # A grid can land exactly on further thresholds depending on libm rounding.
    # Always require the q=pi, omega=0 singularity and truthful metadata.
    if values[(NQ-1)//2][0] is not None or singular != finite_number(meta, 'Divergent threshold points'):
        raise ValueError('Missing threshold marker')
    for iq, column in enumerate(values):
        for a, b in zip(column, values[NQ-1-iq]):
            if (a is None) != (b is None) or (a is not None and not math.isclose(a, b, rel_tol=2e-10, abs_tol=1e-12)):
                raise ValueError('Broken momentum reflection')
    return values


def load_case(solver=None):
    names = ('spectrum', 'continuum')
    texts = (capture_csv_tables([solver, '--momentum-points', NQ, '--points', NW, '--threads', 4], names)
             if solver else {name: (DATA/f'xxx-thermo-{name}.csv').read_text() for name in names})
    values = read_case(texts)
    if solver:
        for name, text in texts.items():
            (DATA/f'xxx-thermo-{name}.csv').write_text(text)
    return values


def plot(values):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    plt.rcParams.update({'svg.hashsalt': 'bethe-xxx-thermo-dsf', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    density = np.array([[np.nan if x is None else x for x in c] for c in values]).T
    q = np.linspace(0, 2, NQ)
    w = np.linspace(0, np.pi, NW)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout='constrained', width_ratios=[1.25, 1])
    mesh = axes[0].pcolormesh(q, w, np.ma.masked_invalid(np.ma.masked_less_equal(density, 0)),
                              shading='nearest', rasterized=True, cmap='magma', norm=LogNorm(.01, 30))
    fig.colorbar(mesh, ax=axes[0], label=r'$S^{zz}_2(q,\omega)$', extend='both')
    fine = np.linspace(0, 2, 601)
    axes[0].plot(fine, np.pi/2*np.abs(np.sin(np.pi*fine)), '--', color='white', lw=.8)
    axes[0].plot(fine, np.pi*np.sin(np.pi*fine/2), '--', color='white', lw=.8)
    axes[0].set(xlabel=r'$q/\pi$', ylabel=r'$\omega/J$', xlim=(0, 2), ylim=(0, np.pi),
                title='Infinite chain: no artificial broadening', xticks=[0,.5,1,1.5,2])
    axes[0].set_facecolor('#f0f0f0')
    for index, label, color in [((NQ-1)//4, r'$q=\pi/2$', '#0072B2'),
                                ((NQ-1)//2, r'$q=\pi$', '#D55E00')]:
        axes[1].plot(w, np.ma.masked_less_equal(density[:, index], 0), label=label, color=color)
    axes[1].set(xlabel=r'$\omega/J$', xlim=(0, np.pi),
                yscale='log', ylim=(.005, 1000), title='Frequency cuts: singular lower thresholds')
    axes[1].legend()
    axes[1].grid(alpha=.15)
    save_svg(fig, FIGURE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable)
    args = parser.parse_args()
    plot(load_case(args.solver))


if __name__ == '__main__':
    main()
