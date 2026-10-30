"""Install into disposable prefixes, relocate, and exercise installed programs."""

import argparse
from contextlib import contextmanager
import csv
from decimal import Decimal
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, timeout=60, **kwargs)
    if result.returncode:
        raise RuntimeError(f'{command!r} exited {result.returncode}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def check_files(prefix, args):
    expected = {Path(args.bindir) / (app + args.suffix) for app in args.apps}
    expected.add(Path(args.docdir) / 'COPYING')
    expected.update(Path(path) for path in args.runtime_files)
    actual = {p.relative_to(prefix) for p in prefix.rglob('*') if p.is_file()}
    if actual != expected:
        raise AssertionError(f'Missing: {expected-actual}; unexpected installed files: {actual-expected}')


def table_rows(document, name):
    table = document['tables'][name]
    columns = [column['id'] for column in table['columns']]
    return [dict(zip(columns, [cell['value'] if isinstance(cell, dict) else cell for cell in row]))
            for row in table['rows']]


def positive_finite(value):
    return value is not None and Decimal(str(value)).is_finite() and Decimal(str(value)) > 0


def check_spinon_calculations(args, prefix, work, env):
    # These small rules exercise the relocated runtime and structured exports.
    # They do not establish sampling, quadrature or rapidity-tail convergence.
    def calculate(app, options, filename):
        executable = str(prefix / args.bindir / (app + args.suffix))
        run([executable, *options, '--threads', '2', '--quiet', '--json', filename], cwd=work, env=env)
        return json.loads((work / filename).read_text())

    moments = 'bethe-xxx-moments-thermo'
    spectrum = 'bethe-xxx-structure-factor-thermo'
    doc = calculate(moments, ['--global', '--diagnostics'], 'two-moments.json')
    rows = table_rows(doc, 'moments')
    assert {(row['sector'], int(row['moment'])) for row in rows} == {('2', 0), ('2', 1)}, rows
    assert all(row['status'] == 'estimated' and positive_finite(row['value']) for row in rows), rows
    assert all(0 < Decimal(str(row['fraction'])) < 1 for row in rows), rows

    doc = calculate(moments, ['--momentum', '1.5', '--sectors', '4', '--precision', 'long-double',
                             '--four-order', '2', '--four-cutoff', '2', '--diagnostics',
                             '--density-relative-tolerance', '1e-4',
                             '--density-absolute-tolerance', '1e-24'], 'four-moments.json')
    rows = table_rows(doc, 'moments')
    assert {(row['sector'], int(row['moment'])) for row in rows} == {('4', 0), ('4', 1)}, rows
    assert all(row['status'] == 'estimated' and positive_finite(row['value']) for row in rows), rows
    assert all(row['sampling_error'] is None for row in table_rows(doc, 'diagnostics'))

    if args.mpreal.upper() not in ('ON', 'TRUE', '1'):
        return
    six = ['--sectors', '6', '--precision', 'bits-384', '--six-cutoff', '8', '--seed', '42',
           '--density-relative-tolerance', '1e-5', '--diagnostics', '--replicate-output']
    momentum = ['--momentum', '1.57079632679489661923132169163975']
    doc = calculate(spectrum, [*six, *momentum, '--omega', '2', '--replicates', '3', '--samples', '8'],
                    'six-spectrum.json')
    rows = table_rows(doc, 'spectrum')
    assert len(rows) == 1 and rows[0]['sector'] == '6', rows
    assert rows[0]['status'] == 'estimated' and positive_finite(rows[0]['density']), rows
    assert positive_finite(table_rows(doc, 'diagnostics')[0]['sampling_error'])
    replicas = table_rows(doc, 'replicates')
    assert len(replicas) == 3 and all(row['failed_sample'] is None for row in replicas), replicas

    for scope, options in [('fixed', momentum), ('global', ['--global'])]:
        doc = calculate(moments, [*six, *options, '--replicates', '2', '--samples', '4', '--covariance'],
                        f'six-{scope}-moments.json')
        rows = table_rows(doc, 'moments')
        assert {(row['sector'], int(row['moment'])) for row in rows} == {('6', 0), ('6', 1)}, rows
        assert all(row['status'] == 'estimated' and positive_finite(row['value']) for row in rows), rows
        assert all(positive_finite(row['sampling_error']) for row in table_rows(doc, 'diagnostics'))
        replicas = table_rows(doc, 'replicates')
        assert len(replicas) == 2 and all(row['failed_sample'] is None for row in replicas), replicas
        covariance = table_rows(doc, 'covariance')
        assert len(covariance) == 4, covariance
        assert all(positive_finite(row['covariance']) for row in covariance
                   if row['moment_i'] == row['moment_j']), covariance


@contextmanager
def preserve_manifests(build_dir):
    # cmake --install records its temporary destinations in the build tree.
    # Do not replace a user's real install history with the test's throwaway one.
    paths = [build_dir / name for name in ('install_manifest.txt', 'install_manifest_BetheRuntime.txt')]
    saved = {path: path.read_bytes() if path.exists() else None for path in paths}
    try:
        yield
    finally:
        for path, content in saved.items():
            if content is None:
                path.unlink(missing_ok=True)
            else:
                path.write_bytes(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cmake', required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--config', default='')
    parser.add_argument('--bindir', required=True)
    parser.add_argument('--docdir', required=True)
    parser.add_argument('--suffix', default='')
    parser.add_argument('--fp128', default='OFF')
    parser.add_argument('--mpreal', default='OFF')
    parser.add_argument('--runtime-files', nargs='*', default=[])
    parser.add_argument('--apps', nargs='+', required=True)
    args = parser.parse_args()
    for destination in (args.bindir, args.docdir):
        if Path(destination).is_absolute() or '..' in Path(destination).parts:
            parser.error('test destinations must stay inside the temporary install prefix')
    # Keep the test hermetic even when launched by a package builder using
    # DESTDIR or symlink-install mode. Never modify a user's install prefix.
    env = dict(os.environ)
    for key in ('DESTDIR', 'CMAKE_INSTALL_MODE', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH'):
        env.pop(key, None)
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    command = [args.cmake, '--install', str(args.build_dir.resolve())]
    if args.config:
        command += ['--config', args.config]
    with preserve_manifests(args.build_dir), tempfile.TemporaryDirectory(prefix='bethe-install-') as directory:
        root = Path(directory)
        prefix = root / 'original prefix'
        run(command + ['--prefix', str(prefix)], env=env)
        check_files(prefix, args)
        relocated = root / 'relocated prefix'
        prefix.rename(relocated)
        work = root / 'calculation'
        work.mkdir()
        # Search only the installed binary directory, not the build tree.
        env['PATH'] = str(relocated / args.bindir)
        for app in args.apps:
            executable = shutil.which(app + args.suffix, path=env['PATH'])
            assert executable, app
            assert app in run([executable, '--help'], cwd=work, env=env)
            assert app in run([executable, '--version'], cwd=work, env=env)
        xxx = str(relocated / args.bindir / ('bethe-xxx-pbc' + args.suffix))
        output = run([xxx, '8', '--excitations', 'all', '--threads', '2'], cwd=work, env=env)
        assert 'converged' in output
        # Each new threaded frontend must find the relocated oneTBB runtime
        # and execute a scan without relying on build-tree library paths.
        threaded_scans = {
            'bethe-lieb-liniger-pbc': ['3', '--length', '4', '--c', '2', '--padding', '1'],
            'bethe-lieb-liniger-obc': ['3', '--length', '4', '--c', '2', '--padding', '1'],
            'bethe-q-boson-pbc': ['4', '--particles', '3', '--eta', '1'],
            'bethe-xxz-qg-obc': ['8', '--delta', '0.6', '--through-lines', '4'],
            'bethe-biquadratic-obc': ['8', '--through-lines', '2'],
        }
        for app, scan in threaded_scans.items():
            executable = str(relocated / args.bindir / (app + args.suffix))
            output = run([executable, *scan, '--excitations', 'all', '--threads', '2'], cwd=work, env=env)
            assert 'converged' in output, (app, output)
        for app, grid in (
            ('bethe-lieb-liniger-dispersion', ['--c', '4']),
            ('bethe-hubbard-dispersion', ['--u', '4']),
            ('bethe-hubbard-dispersion', ['--u', '4', '--density', '0.5']),
        ):
            executable = str(relocated / args.bindir / (app + args.suffix))
            output = run([executable, *grid, '--points', '3', '--threads', '2'], cwd=work, env=env)
            assert 'converged' in output, (app, output)
        precisions = ['fp64', 'long-double'] + (['fp128'] if args.fp128.upper() in ('ON', 'TRUE', '1') else [])
        for precision in precisions:
            output = run([xxx, '4', '--precision', precision, '--format', 'plain'], cwd=work, env=env)
            assert 'Total energy: -2' in output, output
        # Also exercise a BLAS/LAPACK-using model and a relative output path.
        hubbard = str(relocated / args.bindir / ('bethe-hubbard-pbc' + args.suffix))
        for precision in precisions:
            filename = f'hubbard-{precision}.csv'
            output = run([hubbard, '4', '--u', '4', '--precision', precision,
                          '--format', 'plain', '--csv', filename], cwd=work, env=env)
            text = (work / filename).read_text()
            assert '# Outcome: success' in text, output
            rows = list(csv.DictReader(io.StringIO('\n'.join(line for line in text.splitlines() if not line.startswith('#')))))
            assert rows and all(row['status'] == 'converged' for row in rows), rows
        check_spinon_calculations(args, relocated, work, env)
        # The named component must have the same payload as the full install.
        component = root / 'component prefix'
        run(command + ['--prefix', str(component), '--component', 'BetheRuntime'], env=env)
        check_files(component, args)
        if os.name == 'posix':
            staged_env = dict(env, DESTDIR=str(root / 'stage'))
            run(command + ['--prefix', '/bethe-stage', '--component', 'BetheRuntime'], env=staged_env)
            check_files(root / 'stage/bethe-stage', args)
    print(f'Installed and relocated {len(args.apps)} programs; calculations, component and staging checks passed.')


if __name__ == '__main__':
    main()
