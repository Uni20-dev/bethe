#!/usr/bin/env python3
"""Join checked two/four/six exports and plot physical sectors and partial sums.

Raw exports are checked before joining, including six-spinon replicate means
and errors. This establishes identities and completeness, not convergence.
The tutorial's separate order, cutoff and sampling study supplies that evidence.
"""
import argparse
from decimal import localcontext
import hashlib
import json
from pathlib import Path
import re
import shlex

from spinon_tutorial_exports import CONVENTION, close, decimal, read_spectrum
from tutorial_common import save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data/xxx-spinon-slice.json'
FIGURE = ROOT / 'docs/tutorials/figures/xxx-spinon-slice.svg'
MOMENTUM = '1.5707963267948966192313216916397514420985846996875529104874722961539082031431045'
SCHEMA = 'bethe.xxx-spinon-slice.v1'


def join_points(four_runs, six_runs):
    """Join read_spectrum/sampling_points results without interpolating a grid.

    Each run has already passed its sector-specific reader. Different precision
    and integration domains are allowed and retained; q, J and source revisions
    must agree. Four-spinon's T and six-spinon's R are deliberately not equated.
    """
    if not four_runs or not six_runs:
        raise ValueError('Both checked two/four and six inputs are required')
    sources = four_runs + six_runs
    meta = sources[0]['metadata']
    identity = ('Bethe revision', 'Uni20 revision', 'Density convention', 'Channel', 'Exchange J')
    if any(run['metadata'][key] != meta[key] for run in sources for key in identity):
        raise ValueError('Cannot combine different source revisions or physical conventions')
    four, six = {}, {}
    for run_index, run in enumerate(four_runs):
        if run['metadata']['Selected spinon sectors'] != '2,4':
            raise ValueError('Expected the checked two/four sector pair')
        grouped = {}
        for row in run['spectrum']:
            grouped.setdefault(decimal(row['omega']), {})[row['sector']] = row
        for omega, rows in grouped.items():
            if omega in four:
                raise ValueError('Duplicate two/four frequency')
            four[omega] = (rows, run_index)
    for run_index, run in enumerate(six_runs, len(four_runs)):
        for point in run['points']:
            omega = decimal(point['omega'])
            if omega in six:
                raise ValueError('Duplicate six-spinon frequency')
            six[omega] = (point, run_index)
    if four.keys() != six.keys():
        raise ValueError('Sector frequency grids differ; interpolation would change the comparison')
    bits = max(int(run['metadata']['Geometry and output bits']) for run in sources)
    points = []
    with localcontext() as context:
        context.prec = max(80, bits // 3 + 50)
        for omega in sorted(four):
            rows, fi = four[omega]
            sampled, si = six[omega]
            q = rows['2']['q']
            for replica in sampled['replicate_records']:
                common_bits = min(int(sources[i]['metadata']['Geometry and output bits']) for i in (fi, si))
                if not close(q, replica['q'], common_bits):
                    raise ValueError('Sector momenta differ')
            s2, s4, s6 = [decimal(value) for value in (rows['2']['density'], rows['4']['density'], sampled['s6'])]
            points.append({'omega': str(omega), 's2': str(s2), 's4': str(s4), 's6': str(s6),
                           's2_plus_s4': str(s2 + s4), 's2_plus_s4_plus_s6': str(s2 + s4 + s6),
                           's6_sampling_se': sampled['sampling_se'],
                           'largest_draw_fraction': sampled['largest_draw_fraction'],
                           'largest_replicate_fraction': sampled['largest_replicate_fraction'],
                           'four_run': fi, 'six_run': si})
    return points


def provenance(run, raw):
    meta = dict(run['metadata'])
    # Keep the exact command arguments, but an installation command should not
    # contain the original developer's absolute build-directory path.
    command = shlex.split(meta.pop('Command', ''))
    if not command:
        raise ValueError('Missing regeneration command')
    command[0] = Path(command[0]).name
    for index, token in enumerate(command[:-1]):
        if token == '--json':
            command[index + 1] = Path(command[index + 1]).name
    return {'metadata': meta, 'summary': run['summary'], 'command': command,
            'raw_sha256': hashlib.sha256(raw).hexdigest()}


def collect(four_paths, six_paths, *, momentum=MOMENTUM):
    # Owned by the six-spinon tutorial: do not duplicate its statistical checks.
    from plot_xxx_six_spinon_tutorial import sampling_points
    four, six, records = [], [], []
    for paths, destination, reader in ((four_paths, four, read_spectrum),
                                        (six_paths, six, sampling_points)):
        for path in paths:
            raw = path.read_bytes()
            result = reader(json.loads(raw), momentum=momentum, exchange='1')
            destination.append(result)
            record = provenance(result, raw)
            if reader is sampling_points:
                # These are small replicate aggregates, not individual draws.
                record['sampling_points'] = result['points']
            else:
                record['spectrum'] = result['spectrum']
                record['diagnostics'] = result.get('diagnostics', [])
            records.append(record)
    data = {'schema': SCHEMA, 'momentum': momentum, 'exchange': '1', 'channel': 'zz',
            'density_convention': CONVENTION, 'source_runs': records,
            'scope': 'Finite-domain partial spectra; sampling SE excludes cutoff and numerical bias',
            'points': join_points(four, six)}
    return validate(data)


def validate(data):
    if (data['schema'] != SCHEMA or data['channel'] != 'zz' or decimal(data['exchange']) != 1
            or data['density_convention'] != CONVENTION):
        raise ValueError('Wrong slice schema or physical convention')
    decimal(data['momentum'])
    four, six = [], []
    reached_six = False
    for run in data['source_runs']:
        if not re.fullmatch('[0-9a-f]{64}', run['raw_sha256']):
            raise ValueError('Missing raw-export identity')
        if run['metadata']['Channel'] != data['channel'] or run['metadata']['Density convention'] != CONVENTION:
            raise ValueError('Inconsistent source convention')
        meta = run['metadata']
        bits = int(meta['Geometry and output bits'])
        if bits < 24 or not close(meta['Exchange J'], data['exchange'], bits):
            raise ValueError('Inconsistent exchange or precision')
        if any(not re.fullmatch('[0-9a-f]{40}', meta[key]) for key in ('Bethe revision', 'Uni20 revision')):
            raise ValueError('Unidentifiable source revision')
        if run['summary']['Outcome'] != 'success' or run['summary']['Failed points'] != '0':
            raise ValueError('Incomplete source calculation')
        # Rebuild from retained checked records so edits to the plotted sums or
        # individual weights cannot silently disconnect them from their source.
        if 'sampling_points' in run:
            reached_six = True
            if meta['Selected spinon sectors'] != '6':
                raise ValueError('Expected a standalone six-spinon calculation')
            for point in run['sampling_points']:
                if not point['replicate_records']:
                    raise ValueError('Missing six-spinon replicate records')
                if any(not close(row['q'], data['momentum'], bits) for row in point['replicate_records']):
                    raise ValueError('Source momentum differs from declared slice')
            six.append({**run, 'points': run['sampling_points']})
        else:
            if reached_six:
                raise ValueError('Source ordering must be two/four runs followed by six runs')
            if any(not close(row['q'], data['momentum'], bits) for row in run['spectrum']):
                raise ValueError('Source momentum differs from declared slice')
            four.append(run)
    expected = join_points(four, six)
    if data['points'] != expected:
        raise ValueError('Plotted values differ from retained sector records')
    if not expected:
        raise ValueError('Empty spectral slice')
    for point in expected:
        for key in ('s2', 's4', 's6', 's6_sampling_se'):
            if decimal(point[key]) < 0:
                raise ValueError('Negative spectral weight or sampling error')
    return data


def plot(data, output=FIGURE):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import math
    validate(data)
    plt.rcParams.update({'svg.hashsalt': 'bethe-spinon-slice', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout='constrained')
    x = [float(p['omega']) for p in data['points']]
    values = lambda key: [float(p[key]) for p in data['points']]
    axes[0].plot(x, values('s2'), 'o-', markersize=3, color='#0072B2', label='2 spinons')
    axes[0].plot(x, values('s4'), 's-', markersize=3, color='#D55E00', label='4 spinons')
    axes[0].errorbar(x, values('s6'), yerr=values('s6_sampling_se'), fmt='^-', markersize=3,
                     color='#009E73', label='6 spinons; sampling SE')
    for key, label, color in [('s2_plus_s4', '2+4', '#D55E00'),
                               ('s2_plus_s4_plus_s6', '2+4+6', '#009E73')]:
        # Physical zeros are absent on a logarithmic axis. Do not replace them
        # with an invented positive floor or connect through them.
        y = [v if v > 0 else float('nan') for v in values(key)]
        axes[1].plot(x, y, 'o-', markersize=3, color=color, label=label)
    axes[1].set_yscale('log')
    q = float(data['momentum'])
    lower, upper = math.pi / 2 * abs(math.sin(q)), math.pi * abs(math.sin(q / 2))
    for axis in axes:
        axis.axvline(lower, color='#777777', linestyle=':', linewidth=1)
        axis.axvline(upper, color='#777777', linestyle='--', linewidth=1)
        axis.set_xlabel(r'$\omega/J$')
        axis.set_ylabel(r'$J S^{zz}(q,\omega)$')
        axis.grid(alpha=.2)
        axis.legend(frameon=False)
    axes[0].set_title('Individual sectors')
    axes[1].set_title('Partial sums; logarithmic scale')
    save_svg(fig, output)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--four-json', nargs='+', type=Path, help='Complete two/four spectral exports')
    parser.add_argument('--six-json', nargs='+', type=Path, help='Complete six spectral exports with replicates')
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--figure', type=Path, default=FIGURE)
    parser.add_argument('--momentum', default=MOMENTUM)
    parser.add_argument('--check', action='store_true', help='Check data without plotting')
    args = parser.parse_args()
    if bool(args.four_json) != bool(args.six_json):
        parser.error('Import both --four-json and --six-json together')
    if args.four_json:
        data = collect(args.four_json, args.six_json, momentum=args.momentum)
        with args.data.open('x') as output:
            json.dump(data, output, indent=2)
            output.write('\n')
    else:
        data = validate(json.loads(args.data.read_text()))
    if not args.check:
        plot(data, args.figure)
    print(f'Checked {len(data["points"])} matched frequencies; no interpolation or rescaling')


if __name__ == '__main__':
    main()
