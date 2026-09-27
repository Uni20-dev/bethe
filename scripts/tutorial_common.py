"""Small shared helpers for tutorials; importing this needs no plotting packages."""

import csv
import io
import math


def read_csv_export(text, columns, expected, *, row_status='status'):
    metadata = dict(line[2:].split(': ', 1) for line in text.splitlines()
                    if line.startswith('# ') and ': ' in line)
    for key, value in {'Outcome': 'success', 'Status': 'converged', **expected}.items():
        if metadata.get(key) != value:
            raise ValueError(f'{key}: expected {value!r}, got {metadata.get(key)!r}')
    for key in ('Bethe revision', 'Uni20 revision', 'Command'):
        if not metadata.get(key):
            raise ValueError(f'Missing provenance: {key}')
    reader = csv.DictReader(io.StringIO('\n'.join(
        line for line in text.splitlines() if not line.startswith('#'))))
    if reader.fieldnames != columns:
        raise ValueError(f'Unexpected schema: {reader.fieldnames}')
    rows = list(reader)
    for row in rows:
        if (None in row or None in row.values()
                or (row_status is not None and row[row_status] != 'converged')):
            raise ValueError(f'Invalid/incomplete row: {row}')
    return metadata, rows


def finite_number(row, key):
    try:
        value = float(row[key])
    except (KeyError, ValueError, TypeError) as error:
        raise ValueError(f'Missing or invalid {key}') from error
    if not math.isfinite(value):
        raise ValueError(f'Non-finite {key}')
    return value


def save_svg(fig, path):
    output = io.StringIO()
    fig.savefig(output, format='svg', metadata={'Date': None})
    # Normalize Matplotlib's trailing spaces in path coordinates for Git.
    path.write_text('\n'.join(line.rstrip() for line in output.getvalue().splitlines()) + '\n')
