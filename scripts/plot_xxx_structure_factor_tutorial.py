#!/usr/bin/env python3
"""Finite XXX spectral lines, Gaussian plotting only, and checked moment sums."""
import argparse
import math
from pathlib import Path

from tutorial_common import capture_csv_tables, finite_number, read_csv_export, save_svg, solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data'
FIGURES = ROOT / 'docs/tutorials/figures'
N = 64
SPECTRUM = ['state_id', 'momentum_index', 'q', 'gap', 'weight']
MOMENTS = ['momentum_index', 'q', 'weight', 'first_moment', 'full_first_moment', 'first_moment_fraction']


def read_case(texts, n=N):
    count = n*(n+2)//8
    expected = {'Program': 'bethe-xxx-structure-factor', 'Sites': str(n), 'Channel': 'zz',
                'Precision': 'fp64', 'Candidates': str(count), 'Accepted lines': str(count),
                'Converged root states': str(count), 'Full integrated sum rule': '0.25',
                'Family': 'Two-spinon real-root S=1; partial DSF; no rescaling',
                'Fourier convention': 'S_q^+=sum exp(-iqj)S_j^+/sqrt(N); q=P0-Pn',
                'Weight convention': 'S(q,w)=2*pi*sum weight*delta(w-gap); zz=raising/2'}
    meta, raw = read_csv_export(texts['spectrum'], SPECTRUM, expected, row_status=None)
    moment_meta, raw_moments = read_csv_export(texts['moments'], MOMENTS, expected, row_status=None)
    if meta != moment_meta:
        raise ValueError('Mismatched export metadata')
    if len(raw) != count or len(raw_moments) != n:
        raise ValueError('Incomplete two-spinon family/momentum table')
    lines = [{k: finite_number(row, k) for k in SPECTRUM} for row in raw]
    ids = set()
    for line in lines:
        q = line['momentum_index']
        state = line['state_id']
        if (q != int(q) or not 0 < q < n or state != int(state) or not 1 <= state <= count
                or state in ids or line['gap'] <= 0 or line['weight'] <= 0
                or not math.isclose(line['q'], 2*math.pi*q/n, abs_tol=1e-13)):
            raise ValueError('Invalid spectral line')
        ids.add(state)
    e0 = finite_number(meta, 'Ground energy')
    moments = []
    for q, raw_row in enumerate(raw_moments):
        row = {k: finite_number(raw_row, k) for k in MOMENTS[:-1]}
        selected = [line for line in lines if line['momentum_index'] == q]
        full = -4*e0/(3*n)*math.sin(math.pi*q/n)**2
        if (row['momentum_index'] != q or not math.isclose(row['q'], 2*math.pi*q/n, abs_tol=1e-13)
                or not math.isclose(row['weight'], math.fsum(x['weight'] for x in selected), abs_tol=1e-12)
                or not math.isclose(row['first_moment'], math.fsum(x['weight']*x['gap'] for x in selected), abs_tol=1e-12)
                or not math.isclose(row['full_first_moment'], full, abs_tol=1e-12)
                or row['first_moment'] > full+1e-11):
            raise ValueError('Inconsistent moment/sum rule')
        if q == 0:
            if raw_row['first_moment_fraction'] != '':
                raise ValueError('Zero-momentum fraction must be absent')
        elif not math.isclose(finite_number(raw_row, 'first_moment_fraction'), row['first_moment']/full, abs_tol=1e-11):
            raise ValueError('Bad first-moment fraction')
        moments.append(row)
    integrated = math.fsum(x['weight'] for x in lines)/n
    first_fraction = math.fsum(x['first_moment'] for x in moments)/(-2*e0/3)
    if (not 0 < integrated <= .25+1e-12
            or not math.isclose(integrated, finite_number(meta, 'Integrated weight (sum_q/N)'), abs_tol=1e-12)
            or not math.isclose(4*integrated, finite_number(meta, 'Integrated sum-rule fraction'), abs_tol=1e-12)
            or not math.isclose(first_fraction, finite_number(meta, 'Integrated first-moment fraction'), abs_tol=1e-12)):
        raise ValueError('Inconsistent integrated sum rules')
    return lines, moments


def load_case(solver=None):
    names = ('spectrum', 'moments')
    texts = (capture_csv_tables([solver, N], names) if solver else
             {name: (DATA/f'xxx-dsf-n{N}-{name}.csv').read_text() for name in names})
    result = read_case(texts)
    if solver:
        for name, text in texts.items():
            (DATA/f'xxx-dsf-n{N}-{name}.csv').write_text(text)
    return result


def gaussian_spectrum(lines, q_index, omega, eta):
    """S itself, including 2*pi. A normalized kernel, not fitted peak heights."""
    if not math.isfinite(eta) or eta <= 0:
        raise ValueError('Gaussian width must be positive and finite')
    return [2*math.pi*math.fsum(line['weight']*math.exp(-.5*((w-line['gap'])/eta)**2)
                              /(math.sqrt(2*math.pi)*eta)
                              for line in lines if line['momentum_index'] == q_index) for w in omega]


def plot(lines, moments):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    plt.rcParams.update({'svg.hashsalt': 'bethe-xxx-dsf-tutorial', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    points = axes[0].scatter([x['q']/math.pi for x in lines], [x['gap'] for x in lines],
                             c=[x['weight'] for x in lines], s=13, cmap='viridis', norm=LogNorm())
    fig.colorbar(points, ax=axes[0], label=r'Line weight $w_n$ (zz)')
    q = np.linspace(0, 2*math.pi, 401)
    axes[0].plot(q/math.pi, math.pi/2*np.abs(np.sin(q)), '--', color='0.4', lw=1)
    axes[0].plot(q/math.pi, math.pi*np.abs(np.sin(q/2)), '--', color='0.4', lw=1)
    axes[0].set(title=f'N={N}: unbroadened two-spinon lines', xlabel=r'$q/\pi$', ylabel=r'$\Delta E$',
                xlim=(0, 2), ylim=(0, 3.6), xticks=[0,.5,1,1.5,2])
    omega = np.linspace(-.4, 3.6, 1601)
    axes[1].plot(omega, gaussian_spectrum(lines, N//4, omega, .08), color='#0072B2',
                 label=r'Gaussian $\eta=0.08$')
    selected = [x for x in lines if x['momentum_index'] == N//4]
    axes[1].plot([x['gap'] for x in selected], [0]*len(selected), '|', color='#D55E00', ms=10,
                 label='Unbroadened line positions')
    axes[1].set(title=r'Fixed $q=\pi/2$: plotting broadening', xlabel=r'$\omega$',
                ylabel=r'$S^{zz}_{2,\eta}(q,\omega)$', xlim=(-.4, 3.6))
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.15)
    save_svg(fig, FIGURES/'xxx-structure-factor.svg')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=solver_executable, help='Regenerate from this executable on PATH or explicit path')
    args = parser.parse_args()
    plot(*load_case(args.solver))
