#!/usr/bin/env python3
"""Pack and read a finite-domain six-spinon tutorial slice.

The compact file retains per-replicate summaries, not per-observation payloads.
Every read reconstructs the export tables and independently recomputes the
sampling estimates. Raw hashes identify the archived sources; they are not a
claim that an archive is available or has been verified by this reader.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

from plot_xxx_six_spinon_tutorial import sampling_points
from spinon_tutorial_exports import CONVENTION, decimal, read_spectrum

SCHEMA = 'bethe-six-spinon-slice-v1'
CUTOFF = 'Six-spinon absolute rapidity cutoff'
SCOPE = 'Finite-domain estimates; sampling SE excludes rapidity-tail and numerical bias'


def scalar_string(value):
    if not isinstance(value, str):
        raise ValueError('Compact physical values must be decimal strings')
    return decimal(value)


def source_hash(value):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{64}', value):
        raise ValueError('Expected a raw-export SHA256')
    return value


def restored_export(run):
    """Reconstruct only the table envelope, leaving numerical checks shared."""
    tables = {}
    for name in ('spectrum', 'diagnostics', 'replicates'):
        rows = run[name]
        names = list(rows[0]) if rows else []
        if any(set(row) != set(names) for row in rows):
            raise ValueError('Inconsistent compact table columns')
        tables[name] = {'metadata': run['metadata'], 'summary': run['summary'],
                        'columns': [{'id': key} for key in names],
                        'rows': [[row[key] for key in names] for row in rows]}
    return {'tables': tables}


def read_compact_slice(document, *, momentum, exchange='1'):
    """Return checked points in the complete, declared frequency order.

    Different sampling/precision controls at different frequencies are retained
    per run. A slice has one physical convention, cutoff and source revision.
    Missing values must be explicit unavailable records, never survivor means.
    """
    if document.get('schema') != SCHEMA or document.get('density_convention') != CONVENTION:
        raise ValueError('Unknown compact schema or physical convention')
    if (scalar_string(document['momentum']) != decimal(momentum)
            or scalar_string(document['exchange']) != decimal(exchange)
            or decimal(exchange) <= 0):
        raise ValueError('Wrong compact slice momentum or exchange')
    cutoff = scalar_string(document['rapidity_cutoff'])
    if cutoff <= 0:
        raise ValueError('Expected a positive rapidity cutoff')
    frequencies = [scalar_string(value) for value in document['frequencies']]
    if (not frequencies or any(value < 0 for value in frequencies)
            or any(a >= b for a, b in zip(frequencies, frequencies[1:]))):
        raise ValueError('Frequency grid must be nonnegative and strictly increasing')
    by_frequency = {}
    revisions = None
    for index, run in enumerate(document['runs']):
        raw_hash = source_hash(run['raw_sha256'])
        meta = run['metadata']
        if meta['Selected spinon sectors'] != '6' or decimal(meta[CUTOFF]) != cutoff:
            raise ValueError('A compact six-spinon slice must have one sector and cutoff')
        current = (meta['Bethe revision'], meta['Uni20 revision'])
        if revisions is not None and current != revisions:
            raise ValueError('Mixed source revisions in one slice')
        revisions = current
        checked = sampling_points(restored_export(run), momentum=momentum, exchange=exchange)
        for point in checked['points']:
            omega = scalar_string(point['omega'])
            if omega in by_frequency:
                raise ValueError('Duplicate frequency; combining runs requires a separate analysis')
            by_frequency[omega] = dict(point, run_index=index, raw_sha256=raw_hash,
                                       metadata=meta, summary=run['summary'])
    for record in document['unavailable']:
        omega = scalar_string(record['omega'])
        if omega in by_frequency:
            raise ValueError('Duplicate available/unavailable frequency')
        if (record['status'] not in ('failed', 'unavailable')
                or not isinstance(record['reason'], str) or not record['reason'].strip()):
            raise ValueError('Unavailable points need an explicit status and reason')
        if record['raw_sha256'] is not None:
            source_hash(record['raw_sha256'])
        by_frequency[omega] = dict(record, s6=None, sampling_se=None,
                                   replicate_records=[], run_index=None)
    if set(by_frequency) != set(frequencies):
        raise ValueError('Missing or undeclared frequency in compact slice')
    return {'momentum': document['momentum'], 'exchange': document['exchange'],
            'rapidity_cutoff': document['rapidity_cutoff'], 'scope': SCOPE,
            'points': [by_frequency[value] for value in frequencies]}


def compact_slice(raw_exports, *, momentum, frequencies, exchange='1',
                  rapidity_cutoff, unavailable=()):
    """Pack complete six-only exports without copying local command paths.

    The caller declares the full intended grid. Failed/unavailable records must
    be supplied explicitly; omitting a failed export does not remove its node.
    No failed export is converted to an estimated point.
    """
    document = {'schema': SCHEMA, 'density_convention': CONVENTION,
                'momentum': momentum, 'exchange': exchange,
                'rapidity_cutoff': rapidity_cutoff, 'frequencies': list(frequencies),
                'runs': [], 'unavailable': deepcopy(list(unavailable))}
    for raw in raw_exports:
        parsed = json.loads(raw)
        sampling_points(parsed, momentum=momentum, exchange=exchange)
        data = read_spectrum(parsed, momentum=momentum, exchange=exchange)
        # This is the only omitted metadata item: it embeds local output paths.
        # Numerical controls, source identities and timing summaries remain.
        metadata = dict(data['metadata'])
        metadata.pop('Command', None)
        document['runs'].append({'raw_sha256': hashlib.sha256(raw).hexdigest(),
                                 'metadata': metadata, 'summary': data['summary'],
                                 **{key: data[key] for key in ('spectrum', 'diagnostics', 'replicates')}})
    read_compact_slice(document, momentum=momentum, exchange=exchange)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, action='append', default=[],
                        help='Complete six-only export; may be repeated')
    parser.add_argument('--momentum', required=True)
    parser.add_argument('--exchange', default='1')
    parser.add_argument('--cutoff', required=True)
    parser.add_argument('--frequency', action='append', required=True,
                        help='Every requested frequency, in increasing order; repeat per point')
    parser.add_argument('--unavailable', type=Path,
                        help='JSON list of omega/status/reason/raw_sha256 records')
    parser.add_argument('--output', type=Path, required=True, help='New compact file; never overwrite')
    args = parser.parse_args()
    unavailable = json.loads(args.unavailable.read_text()) if args.unavailable else []
    result = compact_slice([path.read_bytes() for path in args.raw], momentum=args.momentum,
                           exchange=args.exchange, rapidity_cutoff=args.cutoff,
                           frequencies=args.frequency, unavailable=unavailable)
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print(f'Packed {len(result["frequencies"])} requested frequencies; '
          f'{len(unavailable)} explicitly unavailable')


if __name__ == '__main__':
    main()
