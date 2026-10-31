#!/usr/bin/env python3
"""Inspect complete six-spinon slice exports and their sampling diagnostics.

The error recomputed here is the standard error of independent replicate means.
It does not include the omitted rapidity domain or prove that rare configurations
have been adequately sampled. Plotting the resolved tutorial data is separate
from the initial proposal and cutoff calibration.
"""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

from spinon_tutorial_exports import decimal, read_spectrum


def integer(value, *, positive=False):
    number = decimal(value)
    if number != number.to_integral_value() or number < (1 if positive else 0):
        raise ValueError('Expected an exact nonnegative count')
    return int(number)


def same_scalar(actual, expected, bits):
    # A tiny spectral weight still needs a relative check. An absolute unit
    # floor would accept a zeroed small error or an entirely missing faint tail.
    a, b = decimal(actual), decimal(expected)
    return abs(a - b) <= Decimal(2) ** (8 - bits) * max(abs(a), abs(b))


def sampling_points(document, *, momentum, exchange='1', expected_metadata=None):
    """Recompute each mean/SE and retain the evidence for influential samples.

    Requires the installed app's diagnostics and replicate tables. Incomplete
    points raise an error; they must remain failed records in a campaign rather
    than becoming a curve assembled from the successful replicates.
    """
    data = read_spectrum(document, momentum=momentum, exchange=exchange,
                         expected_metadata=expected_metadata)
    meta = data['metadata']
    if '6' not in meta['Selected spinon sectors'].split(','):
        raise ValueError('The export has no six-spinon sector')
    bits = integer(meta['Geometry and output bits'], positive=True)
    count = integer(meta['Six-spinon replicates'], positive=True)
    samples = integer(meta['Samples per replicate'], positive=True)
    if count < 2:
        raise ValueError('At least two independent replicates are needed for an SE')
    stream_rule = 'first_replicate=grid_index*replicates; distinct streams at each point'
    if meta['Six-spinon stream identities'] != stream_rule:
        raise ValueError('Unknown replicate identity rule')
    spectra = [r for r in data['spectrum'] if r['sector'] == '6']
    diagnostics = [r for r in data.get('diagnostics', []) if r['sector'] == '6']
    replicates = data.get('replicates', [])
    if len(diagnostics) != len(spectra) or len(replicates) != count * len(spectra):
        raise ValueError('Missing or duplicate six-spinon detail rows')
    points = []
    with localcontext() as context:
        context.prec = max(50, bits // 3 + 25)
        for index, spectrum in enumerate(spectra):
            omega = decimal(spectrum['omega'])
            checks = [r for r in diagnostics if decimal(r['omega']) == omega]
            groups = [r for r in replicates if decimal(r['omega']) == omega]
            ids = [integer(r['replicate']) for r in groups]
            if (len(checks) != 1 or len(ids) != count
                    or set(ids) != set(range(index * count, (index + 1) * count))):
                raise ValueError('Incomplete or duplicate replicate identities')
            groups.sort(key=lambda r: integer(r['replicate']))
            if spectrum['status'] != 'estimated':
                raise ValueError('Expected a sampled interior point')
            means, maxima = [], []
            for row in groups:
                if integer(row['completed_samples']) != samples or row['failed_sample'] is not None:
                    raise ValueError('Incomplete replicate; cannot estimate over survivors')
                mean, maximum = decimal(row['density']), decimal(row['maximum_absolute_density'])
                if mean < 0 or maximum < mean or maximum > samples * mean:
                    raise ValueError('Invalid nonnegative observation summary')
                for key in ('supported_roots', 'density_attempts'):
                    integer(row[key])
                for key in ('contour_diagnostic', 'arithmetic_diagnostic'):
                    if decimal(row[key]) < 0:
                        raise ValueError('Negative numerical diagnostic')
                means.append(mean)
                maxima.append(maximum)
            mean = sum(means) / count
            error = (sum((value - mean) ** 2 for value in means) / (count * (count - 1))).sqrt()
            check = checks[0]
            if not same_scalar(mean, spectrum['density'], bits) or not same_scalar(error, check['sampling_error'], bits):
                raise ValueError('Exported mean/SE differs from independent replicate calculation')
            attempts = sum(integer(row['density_attempts']) for row in groups)
            if attempts != integer(spectrum['evaluations']) or attempts != integer(check['evaluations']):
                raise ValueError('Inconsistent density-attempt count')
            for key in ('contour_diagnostic', 'arithmetic_diagnostic'):
                if decimal(check[key]) < 0:
                    raise ValueError('Negative numerical diagnostic')
            points.append({
                'omega': spectrum['omega'], 's6': str(mean), 'sampling_se': str(error),
                'relative_sampling_se': str(error / mean) if mean else None,
                'largest_draw_fraction': str(max(maxima) / (samples * sum(means))) if mean else None,
                'largest_replicate_fraction': str(max(means) / sum(means)) if mean else None,
                'contour_diagnostic': check['contour_diagnostic'],
                'arithmetic_diagnostic': check['arithmetic_diagnostic'],
                'supported_roots': sum(integer(row['supported_roots']) for row in groups),
                'density_attempts': attempts, 'replicate_records': groups,
                'status': spectrum['status'],
            })
    return {'metadata': data['metadata'], 'summary': data['summary'], 'points': points,
            'scope': 'Finite-domain estimates; sampling SE excludes rapidity-tail and numerical bias'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True, help='Complete selected-sector JSON export')
    parser.add_argument('--momentum', required=True, help='Exact input momentum as a decimal string')
    parser.add_argument('--exchange', default='1')
    parser.add_argument('--output', type=Path, required=True, help='New analysis JSON; never overwrite')
    args = parser.parse_args()
    raw = args.raw.read_bytes()
    result = sampling_points(json.loads(raw), momentum=args.momentum, exchange=args.exchange)
    result['raw_sha256'] = hashlib.sha256(raw).hexdigest()
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print(f'Checked {len(result["points"])} points and all replicate means/errors')


if __name__ == '__main__':
    main()
