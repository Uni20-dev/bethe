"""Install into disposable prefixes, relocate, and exercise installed programs."""

import argparse
from contextlib import contextmanager
import csv
import io
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
    actual = {p.relative_to(prefix) for p in prefix.rglob('*') if p.is_file()}
    if actual != expected:
        raise AssertionError(f'Missing: {expected-actual}; unexpected installed files: {actual-expected}')


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
