#!/usr/bin/env python3
"""Check or reproduce the native four-spinon neighbourhood precision audit."""
import argparse
from decimal import Decimal
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import subprocess
import time

from check_xxx_four_spinon_precision import capture, close, number
from tutorial_common import solver_executable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/validation/four-spinon-neighbourhood-precision-20261031.json'
MOMENTA = ('0.6', '1.0', '1.5', '2.0', '2.5', '3.0')
FREQUENCIES = ('1.0', '1.4', '1.8', '2.2', '2.6', '3.0', '3.5', '4.0', '4.5')
TOLERANCES = ('1e-8', '1e-6', '1e-4')
NATIVE = list(itertools.product(('fp64', 'long-double'), TOLERANCES, MOMENTA, FREQUENCIES))
REFERENCES = [('fp128', '1e-8', '2.5', '1.8'), ('bits-128', '1e-8', '2.5', '1.8'),
              ('bits-128', '1e-8', '1.5', '3.0'), ('bits-128', '1e-8', '1.5', '4.5'),
              ('bits-192', '1e-8', '2.5', '1.8')]
CONTROLS = {'Selected spinon sectors': '4', 'Exchange J': '1', 'Channel': 'zz',
            'Four-spinon order': '4', 'Four-spinon log-radius cutoff': '8',
            'Density absolute tolerance': '1e-30', 'First density residue budget': '64',
            'Second density residue budget': '80', 'First density tail order': '16',
            'Second density tail order': '20', 'Scheduler concurrency limit': '1',
            'Spectral points': '1', 'Recovery seconds per configuration': '60',
            'Density convention': 'S(q,w); moments use dw/(2*pi); raising=2*zz'}
SOURCE_KEYS = ('Bethe revision', 'Uni20 revision', 'Compiler', 'Build type', 'Platform')


def arguments(spec):
    fmt, tol, q, w = spec
    return ['--sectors', '4', '--momentum', q, '--omega', w, '--exchange', '1', '--channel', 'zz',
            '--precision', fmt, '--four-order', '4', '--four-cutoff', '8',
            '--density-relative-tolerance', tol, '--density-absolute-tolerance', '1e-30',
            '--density-terms', '64', '--density-tail-order', '16', '--tolerance', '1e-12',
            '--threads', '1', '--diagnostics', '--format', 'json']


def stem(spec):
    fmt, tol, q, w = spec
    return f'{fmt}-r{tol}-q{q}-w{w}'


def import_case(directory, spec):
    path = directory / (stem(spec)+'.json')
    receipt = json.loads(path.with_suffix('.receipt.json').read_text())
    full = capture(path.read_bytes(), spec[0], receipt['exit_code'])
    m, summary = full['metadata'], full['summary']
    if (full['command'] != ['bethe-xxx-structure-factor-thermo', *arguments(spec)]
            or any(m.get(k) != v for k, v in CONTROLS.items())
            or m['Density relative tolerance'] != spec[1]
            or m['Precision'] != spec[0] or m['Bethe revision'] != receipt['source']):
        raise ValueError('Raw case has different physical or numerical controls')
    if len(full['spectrum']) != 1 or len(full['diagnostics']) != 1:
        raise ValueError('Require one complete selected-sector record')
    row, diagnostic = full['spectrum'][0], full['diagnostics'][0]
    if row['sector'] != '4' or diagnostic['sector'] != '4':
        raise ValueError('Wrong sector')
    if any(row[k] != diagnostic[k] for k in ('q', 'omega', 'evaluations')):
        raise ValueError('Diagnostic belongs to another point')
    case = dict(precision=spec[0], relative_tolerance=spec[1], q=spec[2], omega=spec[3],
                bits=int(m['Geometry and output bits']), comparison_bits=int(m['Density comparison bits']),
                guard_bits=int(m['Density guard bits']), computed_q=row['q'], computed_omega=row['omega'],
                status=row['status'], density=row['density'], evaluations=row['evaluations'],
                contour_diagnostic=diagnostic['contour_diagnostic'],
                arithmetic_diagnostic=diagnostic['arithmetic_diagnostic'], reason=diagnostic['reason'],
                exit_code=full['exit_code'], outcome=summary['Outcome'], failed_points=summary['Failed points'],
                cpu_seconds=summary['Run CPU seconds'], wall_seconds=summary['Elapsed seconds'],
                raw_sha256=full['raw_sha256'])
    source = {k: m[k] for k in SOURCE_KEYS}
    source.update(cpu=receipt['cpu'], binary_sha256=receipt['binary_sha256'])
    return case, source


def validate(data):
    if data['schema'] != 'bethe.four-spinon-neighbourhood.v1' or data['controls'] != CONTROLS:
        raise ValueError('Wrong schema or controls')
    source = data['source']
    for key, length in [('Bethe revision', 40), ('Uni20 revision', 40), ('binary_sha256', 64)]:
        if len(source[key]) != length or any(c not in '0123456789abcdef' for c in source[key]):
            raise ValueError('Missing source identity')
    if not source['cpu']:
        raise ValueError('Missing timing hardware')
    cases = data['cases']
    identities = [(c['precision'], c['relative_tolerance'], c['q'], c['omega']) for c in cases]
    if identities != NATIVE + REFERENCES:
        raise ValueError('Incomplete, duplicate or reordered audit cases')
    policy = data.get('native_comparison_policy', 'binary64-fallback')
    if policy not in ('binary64-fallback', 'intermediate-native'):
        raise ValueError('Unknown native comparison policy')
    native_widths = {c['bits'] for c in cases if c['precision'] == 'long-double'}
    if len(native_widths) != 1:
        raise ValueError('Mixed native arithmetic formats')
    native_bits, = native_widths
    for c in cases:
        fmt, bits = c['precision'], c['bits']
        for field in ('computed_q', 'computed_omega', 'density', 'contour_diagnostic', 'arithmetic_diagnostic'):
            if c[field] is not None and not isinstance(c[field], str):
                raise ValueError('Physical values must retain decimal-string encoding')
        expected = {'fp64': 53, 'fp128': 113, 'bits-128': 128, 'bits-192': 192}
        if bits != expected.get(fmt, bits) or bits < 53:
            raise ValueError('Wrong effective precision')
        low = max(24, bits-64) if fmt.startswith('bits-') else (24 if bits == 53 else 53)
        if policy == 'intermediate-native' and fmt == 'fp128' and 53 < native_bits < 113:
            low = native_bits
        if c['comparison_bits'] != low or c['guard_bits'] != bits:
            raise ValueError('Changed density comparison')
        epsilon = Decimal(2) ** (8-bits)
        if not close(c['computed_q'], c['q'], epsilon) or not close(c['computed_omega'], c['omega'], epsilon):
            raise ValueError('Wrong computed coordinates')
        if len(c['raw_sha256']) != 64 or any(x not in '0123456789abcdef' for x in c['raw_sha256']):
            raise ValueError('Missing raw provenance')
        if number(c['cpu_seconds']) <= 0 or number(c['wall_seconds']) <= 0 or int(c['evaluations']) < 0:
            raise ValueError('Missing measured work')
        if c['status'] == 'numerical_failure':
            if (c['density'] is not None or c['exit_code'] != 2 or c['outcome'] != 'partial'
                    or c['failed_points'] != '1' or not c['reason']):
                raise ValueError('Failed point must remain unavailable')
        elif c['status'] in ('estimated', 'physical_zero'):
            if c['exit_code'] or c['outcome'] != 'success' or c['failed_points'] != '0':
                raise ValueError('Successful point has inconsistent status')
            value = number(c['density'])
            if c['status'] == 'physical_zero':
                if value != 0 or int(c['evaluations']) != 0:
                    raise ValueError('Physical zero has a contribution')
            elif value <= 0 or any(number(c[k]) < 0 for k in ('contour_diagnostic', 'arithmetic_diagnostic')):
                raise ValueError('Missing positive estimate or diagnostics')
        else:
            raise ValueError('Unexpected terminal status')
    # These references check three contrasting points, not every grid value.
    refs = {(c['q'], c['omega']): c for c in cases if c['precision'] == 'bits-128'}
    if any(c['status'] != 'estimated' for c in refs.values()):
        raise ValueError('Missing complete wide reference')
    for c in cases:
        reference = refs.get((c['q'], c['omega']))
        if reference and c['status'] == 'estimated' and not close(c['density'], reference['density'], Decimal('1e-6')):
            raise ValueError('Contrasting-point disagreement exceeds the demonstration budget')
    for q, w in itertools.product(MOMENTA, FREQUENCIES):
        zeros = {c['status'] == 'physical_zero' for c in cases if c['q'] == q and c['omega'] == w}
        if len(zeros) != 1:
            raise ValueError('Physical support changed with precision or tolerance')
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, default=DATA)
    p.add_argument('--raw-dir', type=Path, help='Import retained native cases, or create exports with --solver')
    p.add_argument('--reference-dir', type=Path, help='Separate retained reference directory; defaults to --raw-dir')
    p.add_argument('--solver', type=solver_executable, help='Regenerate all 329 cases, retaining numerical failures')
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    if (a.check and (a.raw_dir or a.reference_dir or a.solver)) or (a.solver and not a.raw_dir) or (a.reference_dir and not a.raw_dir):
        p.error('Use --check alone, or --raw-dir with optional --solver/--reference-dir')
    if a.raw_dir:
        if a.solver:
            if a.reference_dir:
                p.error('Regeneration writes every case to the new --raw-dir')
            a.raw_dir.mkdir(parents=True, exist_ok=False)
            cpu = platform.processor() or platform.machine()
            if Path('/proc/cpuinfo').exists():
                cpu = next((s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                            if s.startswith('model name')), cpu)
            binary_hash = hashlib.sha256(Path(a.solver).read_bytes()).hexdigest()
            for spec in NATIVE + REFERENCES:
                command = [str(a.solver), *arguments(spec)]
                start = time.monotonic()
                path = a.raw_dir / (stem(spec)+'.json')
                timed_out = False
                with path.open('xb') as out, path.with_suffix('.err').open('xb') as err:
                    try:
                        code = subprocess.run(command, stdout=out, stderr=err, timeout=600,
                            env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', UNI20_COLOR='never')).returncode
                    except subprocess.TimeoutExpired:
                        code, timed_out = 124, True
                receipt = dict(exit_code=code, timed_out=timed_out, cpu=cpu,
                    binary_sha256=binary_hash, driver_wall_seconds=time.monotonic()-start)
                try:
                    raw = json.loads(path.read_bytes())
                    receipt['source'] = raw['tables']['spectrum']['metadata']['Bethe revision']
                finally:
                    path.with_suffix('.receipt.json').write_text(json.dumps(receipt)+'\n')
                if code not in (0, 2):
                    raise RuntimeError(f'Incomplete run retained: {path}')
        cases, sources = [], []
        for i, spec in enumerate(NATIVE + REFERENCES):
            directory = a.raw_dir if i < len(NATIVE) else (a.reference_dir or a.raw_dir)
            case, source = import_case(directory, spec)
            cases.append(case)
            sources.append(source)
        if any(s != sources[0] for s in sources):
            raise ValueError('Mixed binaries, sources or timing hardware')
        # Preserve the policy of imported historical exports. Fresh native
        # binary128 uses the available intermediate format; validation still
        # checks that its recorded width matches the native cases in this run.
        fp128 = next(c for c in cases if c['precision'] == 'fp128')
        policy = 'binary64-fallback' if fp128['comparison_bits'] == 53 else 'intermediate-native'
        data = validate(dict(schema='bethe.four-spinon-neighbourhood.v1', controls=CONTROLS,
            source=sources[0], native_comparison_policy=policy,
            scope='Native screen coverage; wide value checks at three contrasting points only. Finite order4/T8 domain, not a physical accuracy bound.', cases=cases))
        a.data.write_text(json.dumps(data, indent=2)+'\n')
    else:
        data = validate(json.loads(a.data.read_text()))
    for fmt, tol in itertools.product(('fp64', 'long-double'), TOLERANCES):
        rows = [c for c in data['cases'][:len(NATIVE)] if c['precision'] == fmt and c['relative_tolerance'] == tol]
        print(fmt, tol, sum(c['status'] == 'estimated' for c in rows), 'estimated;',
              sum(c['status'] == 'numerical_failure' for c in rows), 'failed;',
              sum(c['status'] == 'physical_zero' for c in rows), 'physical zeros')


if __name__ == '__main__':
    main()
