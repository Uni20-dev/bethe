# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Locate the public repository's optional Python reference modules."""
import hashlib
from pathlib import Path
import sys

REFERENCE_DIR = Path(__file__).resolve().parent
REPO_ROOT = REFERENCE_DIR.parents[2]
PUBLIC_SCRIPTS = REPO_ROOT / 'scripts'
if not (PUBLIC_SCRIPTS / 'reference_xxx_vertex_density.py').is_file():
    raise ImportError('QTM controls require the scripts/ references from the same Bethe checkout')

def reference_sources():
    """Reject cached numerical modules from another checkout before reuse."""
    paths = set()
    for name, module in sys.modules.copy().items():
        if not name.startswith(('reference_', 'audit_xxx_qtm')):
            continue
        expected = next((directory / (name + '.py')
                         for directory in (REFERENCE_DIR, PUBLIC_SCRIPTS)
                         if (directory / (name + '.py')).is_file()), None)
        actual = getattr(module, '__file__', None)
        if expected is None or actual is None or Path(actual).resolve() != expected.resolve():
            raise ImportError(f'QTM reference module {name} must come from the same Bethe checkout')
        paths.add(expected.resolve())
    return paths


# A search-path change cannot replace already-imported dependencies: fail
# explicitly instead of silently mixing their functions with this checkout.
reference_sources()
# Prefer these QTM modules, then the same checkout's public controls, ahead
# of PYTHONPATH. These are optional source-checkout references, not a package.
for directory in (PUBLIC_SCRIPTS, REFERENCE_DIR):
    while str(directory) in sys.path:
        sys.path.remove(str(directory))
    sys.path.insert(0, str(directory))


def source_hashes(*extra):
    """Bind executed reference code without local papers or machine paths."""
    paths = {Path(__file__).resolve(), *(Path(p).resolve() for p in extra),
             *reference_sources()}
    return {str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}
