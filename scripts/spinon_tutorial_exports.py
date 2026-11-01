"""Read physical selected-sector exports without losing precision or failures.

These checks establish the export's convention, identity and completion, not
convergence of its quadrature, sampling distribution or finite domain.
"""
from decimal import Decimal, InvalidOperation, localcontext
import math
import re

CONVENTION = 'S(q,w); moments use dw/(2*pi); raising=2*zz'


def decimal(value):
    if isinstance(value, bool) or value is None:
        raise ValueError('Missing numerical value')
    try:
        result = Decimal(str(value))
    except InvalidOperation as error:
        raise ValueError('Invalid numerical value') from error
    if not result.is_finite():
        raise ValueError('Nonfinite numerical value')
    return result


def plotting_float(value):
    """Convert only after unit changes, without turning a nonzero value into zero."""
    number = decimal(value)
    result = float(number)
    if not math.isfinite(result) or (number != 0 and result == 0):
        raise ValueError('Value cannot be represented at plotting precision')
    return result


def table_rows(table, bits):
    names = [column['id'] for column in table['columns']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate table column')
    result = []
    for raw in table['rows']:
        if len(raw) != len(names):
            raise ValueError('Ragged table row')
        values = []
        for value in raw:
            if isinstance(value, dict):
                if set(value) != {'value', 'precision'} or value['precision'] != bits:
                    raise ValueError('Inconsistent scalar precision')
                decimal(value['value'])
                value = value['value']
            values.append(value)
        result.append(dict(zip(names, values)))
    return result


def close(actual, expected, bits, *, scale_floor='1'):
    # Decimal arithmetic must not silently fall back to its 28-digit default.
    with localcontext() as context:
        context.prec = max(50, bits // 3 + 25)
        a, b = decimal(actual), decimal(expected)
        return abs(a - b) <= Decimal(2) ** (8 - bits) * max(Decimal(scale_floor), abs(a), abs(b))


def read_spectrum(document, *, momentum, exchange='1', expected_metadata=None):
    """Return normalized rows plus complete metadata/summary for a successful run.

    The caller selects the momentum and numerical controls explicitly. All
    exported points/replicates are retained; a failed point is an error here,
    not a reason to drop that point from a plot. Numbers remain decimal strings.
    """
    tables = document['tables']
    primary = tables['spectrum']
    metadata, summary = primary['metadata'], primary['summary']
    expected = {'Program': 'bethe-xxx-structure-factor-thermo', 'Channel': 'zz',
                'Density convention': CONVENTION}
    for key, value in (expected_metadata or {}).items():
        if key in expected and expected[key] != value:
            raise ValueError('Requested controls conflict with physical convention')
        expected[key] = value
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError('Wrong physical convention or numerical controls')
    bits = int(metadata['Geometry and output bits'])
    if bits < 24 or not close(metadata['Exchange J'], exchange, bits):
        raise ValueError('Wrong exchange or arithmetic precision')
    for key in ('Bethe revision', 'Uni20 revision'):
        if not re.fullmatch('[0-9a-f]{40}', metadata.get(key, '')):
            raise ValueError('A clean identifiable source revision is required')
    if summary.get('Outcome') != 'success' or summary.get('Failed points') != '0':
        raise ValueError('Incomplete spectral calculation')
    sectors = metadata['Selected spinon sectors'].split(',')
    if not sectors or len(set(sectors)) != len(sectors) or set(sectors) - {'2', '4', '6'}:
        raise ValueError('Invalid selected sectors')
    expected_rows = sectors + (['sum:' + '+'.join(sectors)] if len(sectors) > 1 else [])
    result = {'metadata': dict(metadata), 'summary': dict(summary)}
    for name, table in tables.items():
        if table['metadata'] != metadata or table['summary'] != summary:
            raise ValueError('Mismatched table provenance')
        result[name] = table_rows(table, bits)
    points = {}
    for row in result['spectrum']:
        if not close(row['q'], momentum, bits):
            raise ValueError('Wrong slice momentum')
        omega = decimal(row['omega'])
        point = points.setdefault(omega, {})
        if row['sector'] in point or row['sector'] not in expected_rows:
            raise ValueError('Duplicate or unexpected contribution')
        point[row['sector']] = row
        status = row['status']
        if status not in ('converged', 'estimated', 'physical_zero', 'partial_sum'):
            raise ValueError('Unresolved or singular spectral point')
        value = decimal(row['density'])
        if value < 0 or (status == 'physical_zero' and value != 0):
            raise ValueError('Invalid spectral density or physical zero')
        expected_status = 'partial_sum' if row['sector'].startswith('sum:') else None
        if (status == 'partial_sum') != (expected_status is not None):
            raise ValueError('Incorrect partial-sum label')
    if len(points) != int(metadata['Spectral points']):
        raise ValueError('Missing spectral points')
    with localcontext() as context:
        context.prec = max(50, bits // 3 + 25)
        for point in points.values():
            if set(point) != set(expected_rows):
                raise ValueError('Missing contribution')
            if len(sectors) > 1:
                total = sum(decimal(point[s]['density']) for s in sectors)
                if not close(point[expected_rows[-1]]['density'], total, bits, scale_floor='0'):
                    raise ValueError('Partial sum differs from physical sector sum')
    for name in ('diagnostics', 'replicates'):
        for row in result.get(name, []):
            if not close(row['q'], momentum, bits) or decimal(row['omega']) not in points:
                raise ValueError('Unmatched diagnostic/replicate point')
            if name == 'replicates' and row.get('failed_sample') is not None:
                raise ValueError('Failed replicate in completed spectrum')
    return result
