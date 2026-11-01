"""Validate model precision-recovery accounting in successful spectral exports.

This checks arithmetic policy and retained-observation accounting, not physical
convergence. Failed/cancelled runs remain in their raw archives; this reader does
not turn a subset of their successful observations into an estimated spectrum.
"""
from spinon_tutorial_exports import close, decimal, output_precision_bits

STAGES = ('fp80', 'fp128', 'reference')
COUNTERS = ('attempts', 'accepted', 'retried', 'supported_roots', 'density_attempts')
PRECISION_FIELDS = ('arithmetic format', 'geometry bits', 'density comparison bits', 'density guard bits')
COLUMNS = frozenset(('q', 'omega', 'replicate', 'stage', 'kind', 'attempts', 'accepted', 'retried',
                     'supported_roots', 'density_attempts', 'failure', 'density_status', 'density_failure', 'reason'))
DENSITY_FAILURES = frozenset(('none', 'contour_tail', 'arithmetic_precision', 'tail_domain',
                             'exponent_range', 'special_function_budget', 'rounded_collision'))


def count(value, *, minimum=0):
    number = decimal(value)
    if number != number.to_integral_value() or number < minimum:
        raise ValueError('Invalid recovery count or precision')
    return int(number)


def recovery_policy(metadata):
    """Return the explicit stage policy, or None for a fixed-precision export."""
    output_bits = output_precision_bits(metadata)
    stages = metadata.get('Six-spinon precision stages')
    stage_keys = {key for key in metadata if any(key.startswith(f'Six-spinon {name} ') for name in STAGES)}
    reference_keys = ('Reference density comparison bits', 'Reference density guard bits')
    two_four = 'Two/four-spinon geometry bits'
    if stages is None:
        if stage_keys or any(key in metadata for key in (*reference_keys, two_four)):
            raise ValueError('Stage precision metadata without an explicit recovery policy')
        return None
    names = stages.split(',')  # output_precision_bits validates names and order.
    expected = {f'Six-spinon {name} {field}' for name in names for field in PRECISION_FIELDS}
    if stage_keys != expected or any(key in metadata for key in ('Density comparison bits', 'Density guard bits')):
        raise ValueError('Missing, extra or ambiguously labelled stage precision')
    policy = []
    uniform_bits = count(metadata['Six-spinon uniform grid bits'], minimum=1)
    for name in names:
        prefix = f'Six-spinon {name} '
        geometry = count(metadata[prefix + 'geometry bits'], minimum=24)
        lower = count(metadata[prefix + 'density comparison bits'], minimum=24)
        guard = count(metadata[prefix + 'density guard bits'], minimum=25)
        fmt = metadata[prefix + 'arithmetic format']
        if not lower < guard <= geometry <= output_bits or uniform_bits >= geometry:
            raise ValueError('Stage precision cannot retain its comparisons or uniform grid')
        if name == 'reference':
            if fmt != 'arbitrary binary' or geometry != output_bits:
                raise ValueError('Reference stage must use MPFR at the output precision')
        else:
            expected_format, expected_bits = {'fp80': ('extended80', 64), 'fp128': ('binary128', 113)}[name]
            if (fmt != expected_format or geometry != expected_bits or guard != geometry
                    or lower < 53 or (name == 'fp80' and lower != 53)):
                raise ValueError('Native stage format and precision disagree')
        policy.append({'stage': name, 'format': fmt, 'geometry_bits': geometry,
                       'comparison_bits': lower, 'guard_bits': guard})
    reference = policy[-1]
    if (count(metadata[reference_keys[0]]) != reference['comparison_bits']
            or count(metadata[reference_keys[1]]) != reference['guard_bits']
            or any(stage['guard_bits'] > reference['guard_bits'] for stage in policy[:-1])):
        raise ValueError('Reference density policy cannot retain preceding stages')
    other_sectors = set(metadata['Selected spinon sectors'].split(',')) & {'2', '4'}
    if bool(other_sectors) != (two_four in metadata) or (other_sectors and count(metadata[two_four]) != output_bits):
        raise ValueError('Two/four-sector precision must be the reference precision')
    relative, absolute = (decimal(metadata[key]) for key in ('Density relative tolerance', 'Density absolute tolerance'))
    if relative < 0 or absolute < 0 or relative + absolute <= 0:
        raise ValueError('Invalid recovery tolerances')
    first = (count(metadata['First density residue budget'], minimum=2),
             count(metadata['First density tail order'], minimum=1))
    second = (count(metadata['Second density residue budget'], minimum=2),
              count(metadata['Second density tail order'], minimum=1))
    if first == second or any(a > b for a, b in zip(first, second)):
        raise ValueError('Recovery requires distinct increasing contour controls')
    if decimal(metadata['Recovery seconds per configuration']) <= 0:
        raise ValueError('Recovery needs a positive per-configuration time budget')
    return policy


def recovery_records(data):
    """Reconcile stage/failure rows against every successful physical replicate.

    Input is already decoded by read_spectrum; all table provenance must match.
    A successful run has no terminal failure: each unsuccessful stage must lead
    to the next declared stage on the same sample, and every sample is accepted.
    """
    metadata = data['metadata']
    policy = recovery_policy(metadata)
    if policy is None:
        if 'recovery' in data:
            raise ValueError('Recovery rows without an explicit recovery policy')
        return None
    if 'recovery' not in data:
        raise ValueError('Missing recovery table')
    names = [item['stage'] for item in policy]
    bits = output_precision_bits(metadata)
    samples = count(metadata['Samples per replicate'], minimum=1)
    replicas = {}
    for row in data.get('replicates', []):
        key = decimal(row['omega']), count(row['replicate'])
        if key in replicas:
            raise ValueError('Duplicate physical replicate for recovery accounting')
        if count(row['completed_samples']) != samples or row['failed_sample'] is not None:
            raise ValueError('Recovery cannot hide an incomplete physical replicate')
        replicas[key] = row
    grouped = {key: {name: {'stage': None, 'failures': {}} for name in names} for key in replicas}
    for row in data['recovery']:
        if set(row) != COLUMNS:
            raise ValueError('Unexpected or missing recovery columns')
        key = decimal(row['omega']), count(row['replicate'])
        if key not in grouped or row['stage'] not in names:
            raise ValueError('Unmatched recovery replicate, frequency or stage')
        if not close(row['q'], replicas[key]['q'], bits):
            raise ValueError('Unmatched recovery momentum')
        group = grouped[key][row['stage']]
        attempts, roots, density = (count(row[field]) for field in ('attempts', 'supported_roots', 'density_attempts'))
        if roots > 6 * attempts or density > 2 * roots:
            raise ValueError('Recovery work exceeds the declared observation/contour budget')
        if row['kind'] == 'stage':
            if group['stage'] is not None:
                raise ValueError('Duplicate recovery stage row')
            accepted, retried = count(row['accepted']), count(row['retried'])
            if any(row[field] != '' for field in ('failure', 'density_status', 'density_failure')) or row['reason'] != '':
                raise ValueError('Stage totals must not masquerade as a failure bucket')
            if accepted + retried != attempts:
                raise ValueError('Successful-run stage attempts must be accepted or retried')
            group['stage'] = row
        elif row['kind'] == 'failure':
            if row['accepted'] is not None or row['retried'] is not None or not attempts:
                raise ValueError('Failure buckets need positive attempts and no allocated accept/retry counts')
            if not isinstance(row['reason'], str) or not row['reason'].strip():
                raise ValueError('Failure bucket needs a reason')
            failure, status, density_failure = row['failure'], row['density_status'], row['density_failure']
            if failure == 'density':
                if status != 'exhausted' or (density_failure != '' and density_failure not in DENSITY_FAILURES):
                    raise ValueError('Successful recovery cannot retry a terminal or unknown density failure')
            elif failure in ('geometry', 'exponent_range'):
                if status != '' or density_failure != '':
                    raise ValueError('Geometry/range failure cannot invent density failure tags')
            else:
                raise ValueError('Unknown or terminal failure in a successful recovery')
            failure_key = failure, status, density_failure, row['reason']
            if failure_key in group['failures']:
                raise ValueError('Duplicate recovery failure bucket')
            group['failures'][failure_key] = row
        else:
            raise ValueError('Unknown recovery row kind')
    totals = {name: dict.fromkeys(COUNTERS, 0) for name in names}
    for key, groups in grouped.items():
        expected_attempts = samples
        accepted_roots = accepted_density = accepted_observations = 0
        for name in names:
            group = groups[name]
            stage = group['stage']
            if stage is None:
                raise ValueError('Missing recovery stage row, including an unused stage')
            attempts, accepted, retried, roots, density = (count(stage[field]) for field in COUNTERS)
            if attempts != expected_attempts:
                raise ValueError('Recovery retry chain loses or duplicates observations')
            failed = list(group['failures'].values())
            failed_attempts, failed_roots, failed_density = (
                sum(count(row[field]) for row in failed) for field in ('attempts', 'supported_roots', 'density_attempts'))
            if failed_attempts != retried or failed_roots > roots or failed_density > density:
                raise ValueError('Recovery failure buckets disagree with stage totals')
            successful_roots, successful_density = roots - failed_roots, density - failed_density
            if successful_roots > 6 * accepted or successful_roots % 2 or successful_density != 2 * successful_roots:
                raise ValueError('Accepted stage work is inconsistent with labelled roots and two contour controls')
            expected_attempts = retried
            accepted_observations += accepted
            accepted_roots += successful_roots
            accepted_density += successful_density
            for field in COUNTERS:
                totals[name][field] += count(stage[field])
        if expected_attempts or accepted_observations != samples:
            raise ValueError('Successful recovery has an unresolved final stage')
        physical = replicas[key]
        if accepted_roots != count(physical['supported_roots']) or accepted_density != count(physical['density_attempts']):
            raise ValueError('Accepted recovery work differs from the physical replicate')
    return {'policy': policy, 'rows': data['recovery'], 'totals': totals}
