"""Small shared helpers for tutorials; importing this needs no plotting packages."""

import argparse
import csv
import io
import math
from pathlib import Path
import shutil
import subprocess
import tempfile


def solver_executable(value):
    """Argparse type: accept a command on PATH or an explicit executable path."""
    found = shutil.which(value)
    if found is None:
        raise argparse.ArgumentTypeError(
            f'Executable not found or not executable: {value!r}; '
            'put it on PATH or supply an explicit path')
    return Path(found).resolve()


def capture_csv_tables(command, names):
    """Capture named frontend tables, including export-only auxiliary tables."""
    with tempfile.TemporaryDirectory(prefix='bethe-tutorial-') as directory:
        paths = {name: Path(directory) / f'{name}.csv' for name in names}
        args = [str(arg) for arg in command] + ['--precision', 'fp64', '--format', 'csv']
        for name, path in paths.items():
            args += ['--csv-table', f'{name}={path}']
        subprocess.run(args, check=True, capture_output=True, text=True)
        return {name: path.read_text() for name, path in paths.items()}


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
