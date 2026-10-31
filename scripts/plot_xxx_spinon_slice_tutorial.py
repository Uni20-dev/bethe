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

from spinon_tutorial_exports import CONVENTION, close, decimal, plotting_float, read_spectrum
from tutorial_common import save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/tutorials/data/xxx-spinon-slice.json'
FIGURE = ROOT / 'docs/tutorials/figures/xxx-spinon-slice.svg'
MOMENTUM = '1.5707963267948966192313216916397514420985846996875529104874722961539082031431045'
SCHEMA = 'bethe.xxx-spinon-slice.v1'


def join_points(four_runs, six_points):
    """Join checked sectors on the declared input grid, retaining explicit gaps.

    Six records come from read_compact_slice; their metadata, replicate checks
    and declared missing nodes are owned by that reader. T and R are not equated.
    """
    from xxx_six_spinon_tutorial_data import match_requested_frequency
    if not four_runs or not six_points:
        raise ValueError('Both checked two/four and six inputs are required')
    sources = four_runs + [p for p in six_points if p['s6'] is not None]
    meta = sources[0]['metadata']
    identity = ('Bethe revision', 'Uni20 revision', 'Density convention', 'Channel', 'Exchange J')
    if any(run['metadata'][key] != meta[key] for run in sources for key in identity):
        raise ValueError('Cannot combine different source revisions or physical conventions')
    four, six = {}, {}
    for point in six_points:
        omega = decimal(point['omega'])
        if omega in six:
            raise ValueError('Duplicate six-spinon frequency')
        six[omega] = point
    for run_index, run in enumerate(four_runs):
        if run['metadata']['Selected spinon sectors'] != '2,4':
            raise ValueError('Expected the checked two/four sector pair')
        grouped = {}
        for row in run['spectrum']:
            grouped.setdefault(decimal(row['omega']), {})[row['sector']] = row
        for computed, rows in grouped.items():
            omega = match_requested_frequency(computed, six, int(run['metadata']['Geometry and output bits']))
            if omega in four:
                raise ValueError('Duplicate two/four frequency')
            four[omega] = (rows, run_index)
    if four.keys() != six.keys():
        raise ValueError('Sector frequency grids differ; interpolation would change the comparison')
    bits = max(int(run['metadata']['Geometry and output bits']) for run in sources)
    points = []
    with localcontext() as context:
        context.prec = max(80, bits // 3 + 50)
        for omega in sorted(four):
            rows, fi = four[omega]
            sampled = six[omega]
            s2, s4 = [decimal(rows[sector]['density']) for sector in ('2', '4')]
            point = {'omega': sampled['omega'], 'four_computed_omega': str(rows['2']['omega']),
                     'six_computed_omega': sampled['computed_omega'],
                     's2': str(s2), 's4': str(s4), 's2_plus_s4': str(s2 + s4),
                     'four_run': fi, 'six_run': sampled['run_index'], 's6_status': sampled['status']}
            if sampled['s6'] is None:
                point.update(s6=None, s6_sampling_se=None, s2_plus_s4_plus_s6=None,
                             reason=sampled['reason'])
            else:
                common_bits = min(int(sources[fi]['metadata']['Geometry and output bits']),
                                  int(sampled['metadata']['Geometry and output bits']))
                if any(not close(rows['2']['q'], replica['q'], common_bits)
                       for replica in sampled['replicate_records']):
                    raise ValueError('Sector momenta differ')
                s6 = decimal(sampled['s6'])
                point.update(s6=str(s6), s2_plus_s4_plus_s6=str(s2 + s4 + s6),
                             s6_sampling_se=sampled['sampling_se'],
                             largest_draw_fraction=sampled['largest_draw_fraction'],
                             largest_replicate_fraction=sampled['largest_replicate_fraction'])
            points.append(point)
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


def collect(four_paths, six_document, *, momentum=MOMENTUM):
    from xxx_six_spinon_tutorial_data import read_compact_slice
    six = read_compact_slice(six_document, momentum=momentum, exchange='1')
    four, records = [], []
    for path in four_paths:
        raw = path.read_bytes()
        result = read_spectrum(json.loads(raw), momentum=momentum, exchange='1',
                               expected_metadata={'Selected spinon sectors': '2,4'})
        four.append(result)
        records.append({**provenance(result, raw), 'spectrum': result['spectrum'],
                        'diagnostics': result.get('diagnostics', [])})
    data = {'schema': SCHEMA, 'momentum': momentum, 'exchange': '1', 'channel': 'zz',
            'density_convention': CONVENTION, 'four_runs': records, 'six_slice': six_document,
            'scope': 'Finite-domain partial spectra; sampling SE excludes cutoff and numerical bias',
            'points': join_points(four, six['points'])}
    return validate(data)


def checked_four_run(run, *, momentum, exchange):
    """Restore the small table envelope and reuse the common export checks."""
    if not re.fullmatch('[0-9a-f]{64}', run['raw_sha256']):
        raise ValueError('Missing raw-export identity')
    tables = {}
    for name in ('spectrum', 'diagnostics'):
        rows = run[name]
        names = list(rows[0]) if rows else []
        if any(set(row) != set(names) for row in rows):
            raise ValueError('Inconsistent compact table columns')
        tables[name] = {'metadata': run['metadata'], 'summary': run['summary'],
                        'columns': [{'id': key} for key in names],
                        'rows': [[row[key] for key in names] for row in rows]}
    return read_spectrum({'tables': tables}, momentum=momentum, exchange=exchange,
                         expected_metadata={'Selected spinon sectors': '2,4'})


def validate(data):
    from xxx_six_spinon_tutorial_data import read_compact_slice
    if (data['schema'] != SCHEMA or data['channel'] != 'zz' or decimal(data['exchange']) != 1
            or data['density_convention'] != CONVENTION):
        raise ValueError('Wrong slice schema or physical convention')
    four = [checked_four_run(run, momentum=data['momentum'], exchange=data['exchange'])
            for run in data['four_runs']]
    six = read_compact_slice(data['six_slice'], momentum=data['momentum'], exchange=data['exchange'])
    expected = join_points(four, six['points'])
    if data['points'] != expected:
        raise ValueError('Plotted values differ from retained sector records')
    return data


def plot(data, output=FIGURE):
    validate(data)
    x = [plotting_float(p['omega']) for p in data['points']]
    columns = {key: [plotting_float(p[key]) if p[key] is not None else float('nan')
                     for p in data['points']]
               for key in ('s2', 's4', 's6', 's6_sampling_se', 's2_plus_s4', 's2_plus_s4_plus_s6')}
    q = plotting_float(data['momentum'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import math
    plt.rcParams.update({'svg.hashsalt': 'bethe-spinon-slice', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True, layout='constrained')
    values = columns.__getitem__
    # Physical zeros are absent on a logarithmic axis. Do not replace them
    # with an invented positive floor or connect through them. Missing values
    # remain gaps as well; their distinct statuses stay in the compact record.
    positive = lambda key: [v if v > 0 else float('nan') for v in values(key)]
    axes[0].plot(x, positive('s2'), 'o-', markersize=3, color='#0072B2', label='2 spinons')
    axes[0].plot(x, positive('s4'), 's-', markersize=3, color='#D55E00', label='4 spinons')
    missing = sum(p['s6'] is None for p in data['points'])
    sampling_label = '; unavailable' if missing == len(x) else '; sampling SE'
    axes[0].errorbar(x, positive('s6'), yerr=values('s6_sampling_se'), fmt='^-', markersize=3,
                     color='#009E73', label='6 spinons' + sampling_label)
    for key, label, color in [('s2_plus_s4', '2+4', '#D55E00'),
                               ('s2_plus_s4_plus_s6', '2+4+6', '#009E73')]:
        y = positive(key)
        if key == 's2_plus_s4_plus_s6':
            axes[1].errorbar(x, y, yerr=values('s6_sampling_se'), fmt='o-', markersize=3,
                             color=color, capsize=2, elinewidth=.7,
                             label=label + sampling_label)
        else:
            axes[1].plot(x, y, 'o-', markersize=3, color=color, label=label)
    axes[0].set_yscale('log')
    lower, upper = math.pi / 2 * abs(math.sin(q)), math.pi * abs(math.sin(q / 2))
    for axis in axes:
        axis.axvline(lower, color='#777777', linestyle=':', linewidth=1)
        axis.axvline(upper, color='#777777', linestyle='--', linewidth=1)
        axis.set_xlabel(r'$\omega/J$')
        axis.grid(alpha=.2)
        axis.legend(frameon=False)
    axes[0].set_ylabel(r'$J S^{zz}(q,\omega)$')
    axes[0].set_title('Individual sectors')
    axes[1].set_title('Partial sums')
    title = f'q / π = {q / math.pi:.6g}; finite-domain estimates'
    if missing < len(x):
        title += '; sampling errors only'
    if missing:
        title += f'\nSix-spinon results unavailable at {missing}/{len(x)} requested frequencies'
    fig.suptitle(title, fontsize=11)
    save_svg(fig, output)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--four-json', nargs='+', type=Path, help='Complete two/four spectral exports')
    parser.add_argument('--six-data', type=Path, help='Validated compact six-spinon slice')
    parser.add_argument('--data', type=Path, default=DATA)
    parser.add_argument('--figure', type=Path, default=FIGURE)
    parser.add_argument('--momentum', default=MOMENTUM)
    parser.add_argument('--check', action='store_true', help='Check data without plotting')
    args = parser.parse_args()
    if bool(args.four_json) != bool(args.six_data):
        parser.error('Import both --four-json and --six-data together')
    if args.four_json:
        data = collect(args.four_json, json.loads(args.six_data.read_text()), momentum=args.momentum)
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
