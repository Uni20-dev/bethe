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
                if (mean < 0 or maximum < 0
                        or (maximum < mean and not same_scalar(maximum, mean, bits))
                        or (maximum > samples * mean and not same_scalar(maximum, samples * mean, bits))):
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
            reported_error = decimal(check['sampling_error'])
            # Serialization rounds the means at their own scale. If they are
            # nearly equal, that roundoff can be much larger than eps * SE.
            # Keep the actual data scale, with no absolute unit floor: faint
            # spectra still cannot lose a nonzero SE under this allowance.
            error_roundoff = Decimal(2) ** (8 - bits) * max(*means, error, abs(reported_error))
            if (reported_error < 0 or not same_scalar(mean, spectrum['density'], bits)
                    or abs(error - reported_error) > error_roundoff):
                raise ValueError('Exported mean/SE differs from independent replicate calculation')
            attempts = sum(integer(row['density_attempts']) for row in groups)
            if attempts != integer(spectrum['evaluations']) or attempts != integer(check['evaluations']):
                raise ValueError('Inconsistent density-attempt count')
            for key in ('contour_diagnostic', 'arithmetic_diagnostic'):
                if decimal(check[key]) < 0:
                    raise ValueError('Negative numerical diagnostic')
            points.append({
                'omega': spectrum['omega'], 's6': str(mean), 'sampling_se': str(error),
                'sampling_se_reconstruction_allowance': str(error_roundoff),
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


def paired_cutoffs(lower, upper):
    """Compare checked analyses with identical streams and numerical controls.

    Each argument comes from sampling_points. The proposal's normalized PDF
    changes with its cutoff, so individual differences can have either sign.
    Keep their covariance; do not clamp negative draws or infer an omitted-tail
    bound from a small paired error.
    """
    lm, um = lower['metadata'], upper['metadata']
    cutoff = 'Six-spinon absolute rapidity cutoff'
    required = ('Bethe revision', 'Uni20 revision', 'Program', 'Channel',
                'Density convention', 'Exchange J', 'Selected spinon sectors',
                'Geometry and output bits', 'Samples per replicate',
                'Six-spinon seed', 'Six-spinon uniform grid bits',
                'Six-spinon replicates', 'Six-spinon stream identities')
    prefixes = ('Six-spinon ', 'Density ', 'First density ', 'Second density ')
    keys = set(required) | {key for key in lm.keys() | um.keys() if key.startswith(prefixes)}
    keys.discard(cutoff)
    if any(key not in lm or key not in um or lm[key] != um[key] for key in keys):
        raise ValueError('Cutoff comparison requires common streams, precision, proposal and controls')
    r0, r1 = decimal(lm[cutoff]), decimal(um[cutoff])
    if not 0 < r0 < r1:
        raise ValueError('Require increasing positive rapidity cutoffs')
    if len(lower['points']) != len(upper['points']):
        raise ValueError('Cutoff comparison has different frequency grids')
    pairs = []
    with localcontext() as context:
        context.prec = max(80, integer(lm['Geometry and output bits']) // 3 + 50)
        for a, b in zip(lower['points'], upper['points']):
            if decimal(a['omega']) != decimal(b['omega']):
                raise ValueError('Cutoff comparison has different frequency grids')
            ar, br = a['replicate_records'], b['replicate_records']
            if [x['replicate'] for x in ar] != [x['replicate'] for x in br]:
                raise ValueError('Cutoff comparison has different replicate identities')
            if any(x['q'] != y['q'] for x, y in zip(ar, br)):
                raise ValueError('Cutoff comparison has different momenta')
            x, y = [decimal(row['density']) for row in ar], [decimal(row['density']) for row in br]
            n = len(x)
            mx, my = sum(x) / n, sum(y) / n
            differences = [v - u for u, v in zip(x, y)]
            change = sum(differences) / n
            variance = sum((d - change) ** 2 for d in differences) / (n * (n - 1))
            covariance = sum((u - mx) * (v - my) for u, v in zip(x, y)) / (n * (n - 1))
            vx = sum((u - mx) ** 2 for u in x) / (n * (n - 1))
            vy = sum((v - my) ** 2 for v in y) / (n * (n - 1))
            pairs.append({'omega': a['omega'], 'change': str(change), 'paired_sampling_se': str(variance.sqrt()),
                          'covariance_of_means': str(covariance),
                          'replicate_correlation': str(covariance / (vx * vy).sqrt()) if vx and vy else None,
                          'fraction_of_upper_estimate': str(change / my) if my else None,
                          'replicate_differences': [str(d) for d in differences]})
    return {'lower_cutoff': str(r0), 'upper_cutoff': str(r1), 'points': pairs,
            'scope': 'Paired finite-domain difference; neither an extrapolation nor a bound on the omitted tail'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True, help='Complete selected-sector JSON export')
    parser.add_argument('--momentum', required=True, help='Exact input momentum as a decimal string')
    parser.add_argument('--exchange', default='1')
    parser.add_argument('--compare-raw', type=Path, help='Larger-cutoff export with the same streams and controls')
    parser.add_argument('--output', type=Path, required=True, help='New analysis JSON; never overwrite')
    args = parser.parse_args()
    raw = args.raw.read_bytes()
    result = sampling_points(json.loads(raw), momentum=args.momentum, exchange=args.exchange)
    result['raw_sha256'] = hashlib.sha256(raw).hexdigest()
    if args.compare_raw:
        other = args.compare_raw.read_bytes()
        checked = sampling_points(json.loads(other), momentum=args.momentum, exchange=args.exchange)
        result['cutoff_comparison'] = paired_cutoffs(result, checked)
        result['comparison_raw_sha256'] = hashlib.sha256(other).hexdigest()
    with args.output.open('x') as output:
        json.dump(result, output, indent=2)
        output.write('\n')
    print(f'Checked {len(result["points"])} points and all replicate means/errors')


if __name__ == '__main__':
    main()
