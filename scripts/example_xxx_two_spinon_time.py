#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Fresh finite-chain and thermodynamic two-spinon time signals at q=pi/2."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from plot_xxx_structure_factor_tutorial import read_case as read_finite
from spectral_time_transform import coarser_indices, weights
from tutorial_common import capture_csv_tables, finite_number, read_csv_export, solver_executable

Q = math.pi/2
LOWER, UPPER = math.pi/2, math.pi/math.sqrt(2)
COLUMNS = ['q', 'omega', 'sector', 'density', 'estimated_error', 'evaluations', 'status']


def read_segment(texts, count, begin, end):
    """Reject failed/missing rows rather than silently changing the support."""
    expected = {'Program': 'bethe-xxx-structure-factor-thermo', 'Precision': 'fp64',
                'Channel': 'zz', 'Exchange J': '1', 'Spectral points': str(count),
                'Failed points': '0', 'Divergent threshold points': '0',
                'Family': 'Thermodynamic two-spinon contribution; partial DSF; no rescaling',
                'Density convention': 'S(q,w) itself; integrals use dw/(2*pi); raising=2*zz'}
    meta, rows = read_csv_export(texts['spectrum'], COLUMNS, expected, row_status=None)
    other, edges = read_csv_export(texts['continuum'], ['q', 'lower', 'upper'],
                                   expected, row_status=None)
    if meta != other or len(rows) != count or len(edges) != 1:
        raise ValueError('Incomplete or mismatched spectrum/continuum export')
    if any(not math.isclose(finite_number(edges[0], key), value, rel_tol=0, abs_tol=2e-14)
           for key, value in [('q', Q), ('lower', LOWER), ('upper', UPPER)]):
        raise ValueError('Unexpected continuum')
    values = []
    for i, row in enumerate(rows):
        omega, density, error = (finite_number(row, key)
                                 for key in ('omega', 'density', 'estimated_error'))
        target = begin+(end-begin)*i/(count-1)
        if (row['sector'] != '2' or density < 0 or error < 0
                or not math.isclose(finite_number(row, 'q'), Q, rel_tol=0, abs_tol=2e-14)
                or not math.isclose(omega, target, rel_tol=0, abs_tol=2e-14)):
            raise ValueError('Unexpected spectral sector, grid or value')
        if row['status'] == 'outside_continuum':
            if abs(omega-UPPER) > 2e-14 or density != 0:
                raise ValueError('Unexpected outside-continuum point')
        elif row['status'] != 'converged':
            raise ValueError('Failed or singular spectral value')
        if omega < UPPER-2e-14 and density == 0:
            raise ValueError('Missing interior spectral weight')
        values.append((omega, density, error))
    return meta, np.asarray(values)


def time_signals(lines, spectrum, times, phase_order=16):
    """Keep discrete weights and continuum density measures distinct."""
    omega, density, errors = np.asarray(spectrum).T
    times = np.asarray(times, dtype=float)
    gaps = np.array([line['gap'] for line in lines])
    strengths = np.array([line['weight'] for line in lines])
    inside = (gaps >= omega[0]) & (gaps <= omega[-1])
    finite = np.exp(-1j*times[:, None]*gaps) @ strengths
    finite_window = np.exp(-1j*times[:, None]*gaps[inside]) @ strengths[inside]
    options = dict(offsets=omega-LOWER, order=phase_order)
    W = weights(omega, times, **options)
    continuum = W @ density
    higher = weights(omega, times, offsets=omega-LOWER, order=2*phase_order) @ density
    coarse = weights(omega, times, indices=coarser_indices(len(omega), 2), **options) @ density
    moments = {}
    for name, gap, weight in [('finite_all', gaps, strengths),
                               ('finite_window', gaps[inside], strengths[inside])]:
        moments[name] = {'m0': float(math.fsum(weight)),
                         'm1': float(math.fsum(gap*weight))}
    moments['thermo_window'] = {
        f'm{r}': float((weights(omega, [0.], energy_power=r, **options) @ density)[0].real)
        for r in (0, 1)}
    columns = dict(time=times, finite_all_real=finite.real, finite_all_imag=finite.imag,
                   finite_window_real=finite_window.real, finite_window_imag=finite_window.imag,
                   thermo_window_real=continuum.real, thermo_window_imag=continuum.imag,
                   phase_order_change_abs=np.abs(higher-continuum),
                   coarse_grid_change_abs=np.abs(coarse-continuum),
                   propagated_kernel_error_bound=np.abs(W) @ errors)
    return columns, moments


def calculate(args):
    directory = args.output
    directory.mkdir(parents=True, exist_ok=False)

    def capture(command, stem):
        texts = capture_csv_tables(command, ['spectrum', 'moments'] if stem == 'finite'
                                   else ['spectrum', 'continuum'])
        for name, text in texts.items():
            (directory/f'{stem}-{name}.csv').write_text(text)
        return texts

    finite_texts = capture([args.finite_solver, args.sites, '--channel', 'zz', '--threads', 1], 'finite')
    all_lines, finite_moments = read_finite(finite_texts, args.sites)
    lines = [line for line in all_lines if line['momentum_index'] == args.sites//4]
    # Each CLI call has a uniform grid; geometrically spaced segment endpoints
    # resolve the lower threshold without changing either solver's public API.
    segments = max(1, math.ceil(math.log10((UPPER-LOWER)/args.lower_offset)))
    offsets = np.geomspace(args.lower_offset, UPPER-LOWER, segments+1)
    bounds = LOWER+offsets
    bounds[-1] = UPPER
    parts, provenance = [], []
    for i, (begin, end) in enumerate(zip(bounds[:-1], bounds[1:])):
        texts = capture([args.thermo_solver, '--sectors', 2, '--channel', 'zz', '--momentum', Q,
                         '--omega-min', begin, '--omega-max', end, '--points', args.points,
                         '--threads', 1], f'thermo-{i:02}')
        meta, part = read_segment(texts, args.points, begin, end)
        if parts:
            if not np.array_equal(parts[-1][-1], part[0]):
                raise ValueError('Inconsistent repeated segment endpoint')
            part = part[1:]
        parts.append(part)
        provenance.append(meta)
    spectrum = np.concatenate(parts)
    times = np.linspace(0, args.time_max, args.time_points)
    columns, moments = time_signals(lines, spectrum, times, args.phase_order)
    with (directory/'time.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    for r, key in [(0, 'weight'), (1, 'first_moment')]:
        if not math.isclose(moments['finite_all'][f'm{r}'], finite_moments[args.sites//4][key],
                            rel_tol=2e-12, abs_tol=1e-14):
            raise ValueError('Finite-line moments disagree with the solver table')
    report = {
        'model': 'XXX, J=1, zero field and temperature; longitudinal two-spinons only',
        'sites': args.sites, 'momentum': Q, 'finite_lines': len(lines),
        'density_measure': 'C(t)=integral S(q,omega)*exp(-i*omega*t)*domega/(2*pi)',
        'line_measure': 'C_N(t)=sum weight*exp(-i*gap*t)',
        'interpolation': 'piecewise linear 2*x*S in x=sqrt(omega-omega_lower)',
        'thermo_window': [float(spectrum[0, 0]), float(spectrum[-1, 0])],
        'omitted_threshold_interval': [LOWER, float(spectrum[0, 0])],
        'missing_threshold_weight': None,
        'frequency_nodes': len(spectrum), 'phase_order': args.phase_order,
        'moments': moments,
        'finite_full_first_moment': finite_moments[args.sites//4]['full_first_moment'],
        'thermo_full_first_moment': (2/3)*(math.log(2)-1/4)*(1-math.cos(Q)),
        'max_phase_order_change': float(max(columns['phase_order_change_abs'])),
        'max_coarse_grid_change': float(max(columns['coarse_grid_change_abs'])),
        'max_propagated_kernel_error_bound': float(max(columns['propagated_kernel_error_bound'])),
        'accuracy_note': 'Order/grid changes are diagnostics, not bounds on interpolation or omitted weight. '
                         'The kernel bound propagates reported pointwise error estimates only. '
                         'Both spectra omit higher spinon sectors; no finite-size extrapolation is made.',
        'source_exports': sorted(path.name for path in directory.glob('*-*.csv')),
        'thermo_provenance': provenance,
    }
    (directory/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--finite-solver', type=solver_executable, default='bethe-xxx-structure-factor')
    parser.add_argument('--thermo-solver', type=solver_executable, default='bethe-xxx-structure-factor-thermo')
    parser.add_argument('--output', type=Path, default=Path('xxx-time-example'), help='new output directory')
    parser.add_argument('--sites', type=int, default=16, help='multiple of four, at least four')
    parser.add_argument('--points', type=int, default=129, help='frequency points per geometric segment')
    parser.add_argument('--lower-offset', type=float, default=1e-6, help='omitted threshold width in J units')
    parser.add_argument('--phase-order', type=int, default=16)
    parser.add_argument('--time-max', type=float, default=20., help='maximum Jt')
    parser.add_argument('--time-points', type=int, default=81)
    args = parser.parse_args()
    if (args.sites < 4 or args.sites % 4 or args.points < 3 or args.time_points < 2
            or args.phase_order < 2 or not 1e-10 <= args.lower_offset < UPPER-LOWER
            or not math.isfinite(args.time_max) or args.time_max <= 0):
        parser.error('invalid sites/grid; require 1e-10 <= lower-offset < upper-lower')
    report = calculate(args)
    print(json.dumps({key: report[key] for key in ('moments', 'max_phase_order_change',
                     'max_coarse_grid_change', 'omitted_threshold_interval')}, indent=2))
    print(f'Wrote {args.output}/time.csv and report.json; source exports retained.')


if __name__ == '__main__':
    main()
