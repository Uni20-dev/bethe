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
# Keep this package first: it contains the QTM implementations. The added
# directory contains already-public vertex controls, not an installed runtime.
if str(PUBLIC_SCRIPTS) not in sys.path:
    sys.path.append(str(PUBLIC_SCRIPTS))


def source_hashes(*extra):
    """Bind executed reference code without local papers or machine paths."""
    paths = {Path(__file__).resolve(), *(Path(p).resolve() for p in extra)}
    for name, module in sys.modules.copy().items():
        if (name.startswith(('reference_', 'audit_xxx_qtm'))
                and getattr(module, '__file__', None)):
            path = Path(module.__file__).resolve()
            if REPO_ROOT in path.parents:
                paths.add(path)
    return {str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}
