#!/usr/bin/env python3
"""Pack and read a finite-domain six-spinon tutorial slice.

The compact file retains per-replicate summaries, not per-observation payloads.
Every read reconstructs the export tables and independently recomputes the
sampling estimates. Raw hashes identify the archived sources; they are not a
claim that an archive is available or has been verified by this reader.
"""
import argparse
from copy import deepcopy
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import re

from plot_xxx_six_spinon_tutorial import sampling_points
from spinon_tutorial_exports import CONVENTION, decimal, read_spectrum

SCHEMA = 'bethe-six-spinon-slice-v1'
SET_SCHEMA = 'bethe-six-spinon-slice-set-v1'
CUTOFF = 'Six-spinon absolute rapidity cutoff'
SCOPE = 'Finite-domain estimates; sampling SE excludes rapidity-tail and numerical bias'
NODE_CONTROLS = ('Six-spinon seed', 'Six-spinon replicates', 'Samples per replicate')


def cohort_metadata(metadata):
    """Bind all exported settings except invocation, date and node-specific counts."""
    return {key: value for key, value in metadata.items()
            if key not in ('Command', 'Date', *NODE_CONTROLS)}


def scalar_string(value):
    if not isinstance(value, str):
        raise ValueError('Compact physical values must be decimal strings')
    return decimal(value)


def source_hash(value):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{64}', value):
        raise ValueError('Expected a raw-export SHA256')
    return value


def match_requested_frequency(value, requested, bits):
    """Identify one declared input node within its recorded rounding allowance.

    Allow eight relative machine epsilons for grid arithmetic and serialization.
    Return a Decimal key. No interpolation, nearest-neighbour fallback or unit
    absolute floor: ambiguous grids and genuinely different frequencies fail.
    """
    if isinstance(bits, bool) or not isinstance(bits, int) or bits < 24:
        raise ValueError('Invalid frequency precision')
    candidates = [decimal(node) for node in requested]
    if len(set(candidates)) != len(candidates):
        raise ValueError('Duplicate requested frequency')
    actual = decimal(value)
    with localcontext() as context:
        context.prec = max(50, bits // 3 + 25)
        allowance = Decimal(2) ** (4 - bits)
        matches = [node for node in candidates
                   if abs(actual - node) <= allowance * max(abs(actual), abs(node))]
    if len(matches) != 1:
        raise ValueError('Ambiguous requested frequency' if matches else 'Undeclared frequency')
    return matches[0]


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
    per run. Returned omega is the declared input; computed_omega preserves its
    rounded export coordinate. A slice has one convention, cutoff and source
    revision. Missing values are explicit records, never survivor means.
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
    requested_text = dict(zip(frequencies, document['frequencies']))
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
            computed_omega = point['omega']
            omega = match_requested_frequency(computed_omega, frequencies, int(meta['Geometry and output bits']))
            if omega in by_frequency:
                raise ValueError('Duplicate frequency; combining runs requires a separate analysis')
            by_frequency[omega] = dict(point, omega=requested_text[omega], computed_omega=computed_omega,
                                       run_index=index, raw_sha256=raw_hash,
                                       metadata=meta, summary=run['summary'])
    for record in document['unavailable']:
        omega = scalar_string(record['omega'])
        if omega not in requested_text:
            raise ValueError('Undeclared unavailable frequency')
        if omega in by_frequency:
            raise ValueError('Duplicate available/unavailable frequency')
        if (record['status'] not in ('failed', 'unavailable')
                or not isinstance(record['reason'], str) or not record['reason'].strip()):
            raise ValueError('Unavailable points need an explicit status and reason')
        if record['raw_sha256'] is not None:
            source_hash(record['raw_sha256'])
        by_frequency[omega] = dict(record, omega=requested_text[omega], computed_omega=None,
                                   s6=None, sampling_se=None,
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


def read_slice_set(document, *, momentum, exchange='1'):
    """Read disjoint, individually strict source cohorts on one declared grid.

Compatibility is an explicit reviewed assertion with evidence identities;
it does not waive any per-run numerical, stream, source or convention check.
Unavailable nodes belong to a declared cohort and remain unavailable.
"""
    from spinon_tutorial_provenance import compatibility_sources, source_identity, validate_sources
    if document.get('schema') != SET_SCHEMA or document.get('density_convention') != CONVENTION:
        raise ValueError('Unknown slice-set schema or physical convention')
    if (scalar_string(document['momentum']) != decimal(momentum)
            or scalar_string(document['exchange']) != decimal(exchange) or decimal(exchange) <= 0):
        raise ValueError('Wrong slice-set momentum or exchange')
    cutoff = scalar_string(document['rapidity_cutoff'])
    if cutoff <= 0:
        raise ValueError('Expected a positive rapidity cutoff')
    frequencies = [scalar_string(value) for value in document['frequencies']]
    if (not frequencies or any(value < 0 for value in frequencies)
            or any(a >= b for a, b in zip(frequencies, frequencies[1:]))):
        raise ValueError('Slice-set frequency grid must be nonnegative and strictly increasing')
    requested = dict(zip(frequencies, document['frequencies']))
    compatible = compatibility_sources(document['source_compatibility'])
    contracts = document['source_compatibility'].get('cohort_controls')
    if not isinstance(contracts, dict):
        raise ValueError('Slice-set compatibility needs explicit cohort controls')
    cohorts = document['cohorts']
    if not isinstance(cohorts, list) or len(cohorts) < 2:
        raise ValueError('A slice set needs at least two declared cohorts')
    ids, declared_sources, points = set(), set(), {}
    for cohort in cohorts:
        if not isinstance(cohort, dict):
            raise ValueError('Invalid source cohort')
        name = cohort.get('id')
        if not isinstance(name, str) or not re.fullmatch('[a-z0-9][a-z0-9._-]*', name) or name in ids:
            raise ValueError('Cohort IDs must be unique nonempty identifiers')
        ids.add(name)
        source = source_identity(cohort['source'])
        if source not in compatible:
            raise ValueError('Cohort source is outside the compatibility record')
        declared_sources.add(source)
        contract = contracts.get(name)
        if not isinstance(contract, dict) or source_identity(contract.get('source')) != source:
            raise ValueError('Cohort source differs from its reviewed control record')
        for key in ('binary_sha256', 'manifest_sha256'):
            source_hash(contract.get(key))
        expected_metadata = contract.get('metadata')
        if not isinstance(expected_metadata, dict) or not expected_metadata:
            raise ValueError('Missing bound cohort metadata')
        fixed = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
                 'Density convention': CONVENTION, 'Selected spinon sectors': '6'}
        if (any(expected_metadata.get(key) != value for key, value in fixed.items())
                or scalar_string(expected_metadata['Exchange J']) != decimal(exchange)
                or scalar_string(expected_metadata[CUTOFF]) != cutoff
                or cohort_metadata(expected_metadata) != expected_metadata):
            raise ValueError('Bound cohort physics or common metadata is invalid')
        validate_sources([expected_metadata], document['source_compatibility'])
        if source_identity(expected_metadata) != source:
            raise ValueError('Bound metadata has a different source')
        nodes = contract.get('nodes')
        if not isinstance(nodes, list) or not nodes:
            raise ValueError('Missing declared cohort node controls')
        expected_nodes = {}
        for node in nodes:
            if not isinstance(node, dict):
                raise ValueError('Invalid cohort node controls')
            omega = scalar_string(node['omega'])
            if omega in expected_nodes:
                raise ValueError('Duplicate cohort node controls')
            controls = node.get('metadata')
            if not isinstance(controls, dict) or set(controls) != set(NODE_CONTROLS):
                raise ValueError('Node controls must bind seed, replicates and samples')
            for key, value in controls.items():
                count = scalar_string(value)
                minimum = 0 if key == NODE_CONTROLS[0] else 2 if key == NODE_CONTROLS[1] else 1
                if count != count.to_integral_value() or count < minimum:
                    raise ValueError('Invalid cohort stream/sample controls')
            expected_nodes[omega] = controls
        data = read_compact_slice(cohort['slice'], momentum=momentum, exchange=exchange)
        if scalar_string(data['rapidity_cutoff']) != cutoff:
            raise ValueError('Cohort rapidity cutoffs differ')
        if any(source_identity(run['metadata']) != source for run in cohort['slice']['runs']):
            raise ValueError('Declared cohort source differs from its actual runs')
        if any(cohort_metadata(run['metadata']) != expected_metadata for run in cohort['slice']['runs']):
            raise ValueError('Cohort numerical controls differ from the reviewed record')
        if {decimal(point['omega']) for point in data['points']} != set(expected_nodes):
            raise ValueError('Cohort grid differs from its reviewed node controls')
        for point in data['points']:
            omega = decimal(point['omega'])
            if omega not in requested:
                raise ValueError('Cohort frequency was not requested')
            if omega in points:
                raise ValueError('Duplicate frequency across cohorts; replacement is not allowed')
            if point['s6'] is not None and any(point['metadata'].get(key) != value
                                             for key, value in expected_nodes[omega].items()):
                raise ValueError('Cohort seed or sampling controls differ from the reviewed record')
            points[omega] = dict(point, omega=requested[omega], cohort_id=name,
                                source=dict(zip(('Bethe revision', 'Uni20 revision'), source)))
    if declared_sources != compatible:
        raise ValueError('Compatibility record contains an undeclared cohort source')
    if ids != set(contracts):
        raise ValueError('Compatibility record contains undeclared cohort controls')
    if set(points) != set(frequencies):
        raise ValueError('Missing frequency from source cohorts')
    return {'momentum': document['momentum'], 'exchange': document['exchange'],
            'rapidity_cutoff': document['rapidity_cutoff'], 'scope': SCOPE,
            'source_compatibility': deepcopy(document['source_compatibility']),
            'points': [points[value] for value in frequencies]}


def read_slice(document, *, momentum, exchange='1'):
    """Dispatch one strict source slice or an explicit reviewed slice set."""
    if document.get('schema') == SET_SCHEMA:
        return read_slice_set(document, momentum=momentum, exchange=exchange)
    return read_compact_slice(document, momentum=momentum, exchange=exchange)


def compact_slice_set(cohorts, compatibility, *, momentum, frequencies,
                      exchange='1', rapidity_cutoff):
    """Assemble without replacing points or changing any cohort's records.

Each cohort has an id, exact source metadata and a strict v1 compact slice.
Local run indices remain local to cohort_id in the checked output.
"""
    document = {'schema': SET_SCHEMA, 'density_convention': CONVENTION,
                'momentum': momentum, 'exchange': exchange, 'rapidity_cutoff': rapidity_cutoff,
                'frequencies': list(frequencies), 'cohorts': deepcopy(cohorts),
                'source_compatibility': deepcopy(compatibility)}
    read_slice_set(document, momentum=momentum, exchange=exchange)
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
