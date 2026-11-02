#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ian McCulloch
"""Small, reproducible precision ladder for the two isotropic QTM references."""
import argparse
import json
import platform
from pathlib import Path
import sys
import time

import mpmath as mp

from _reference_support import REFERENCE_DIR, source_hashes as public_source_hashes

import reference_xxx_qtm as qtm


# Decimal strings are parsed afresh INSIDE each requested precision context.
FIXTURES = {
    "ell1": (["-.7"], [".3"], "-.4", ".6"),
    "ell2": (["-.7", ".8"], ["-.2", "1.4"], "-.4", ".6"),
    "ell3": (["-.7", ".8", "1.8"], ["-.2", "1.4", "2.4"], "-.4", ".6"),
    "near-cross-pole": (["-.7", ".8"], ["-.699", "1.4"], "-.4", ".6"),
    "near-sinh-zero": (["-.7", ".8"], ["-.2", ".301"], "-.4", ".6"),
    "clustered": (["-.7", "-.67", ".8"], ["-.2", ".3", "1.4"], "-.4", ".6"),
    "other-contours": (["-.7", ".8"], ["-.2", "1.4"], "-.3", ".45"),
    "symmetric-offset": (["-.7", ".8", "1.8"], ["-.2", "1.4", "2.4"], "-.55", ".55"),
}


def encode(x):
    """Lossless-to-working-precision decimal strings, never JSON float casts."""
    if isinstance(x, mp.mpc):
        return dict(real=mp.nstr(x.real, mp.mp.dps), imag=mp.nstr(x.imag, mp.mp.dps))
    if isinstance(x, mp.mpf):
        return mp.nstr(x, mp.mp.dps)
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    return x


def source_hashes():
    return public_source_hashes(__file__, REFERENCE_DIR/"test_reference_xxx_qtm.py")

def comparison(name, digits, maxdegree):
    x, y, iu, iv = FIXTURES[name]
    with mp.workdps(digits):
        u, v = [mp.mpc(z, iu) for z in x], [mp.mpc(z, iv) for z in y]
        start = time.perf_counter()
        jmp = qtm.jmp(u, v, maxdegree=maxdegree)
        middle = time.perf_counter()
        prl = qtm.prl(u, v)
        finish = time.perf_counter()
        error = abs(jmp["amplitude"]-prl["amplitude"])/max(
            mp.mpf(1), abs(jmp["amplitude"]), abs(prl["amplitude"]))
        tolerance = mp.power(10, -(digits-10))
        row = dict(
            fixture=name, ell=len(u), digits=digits, holes=u, particles=v,
            passed=bool(error < tolerance), scaled_difference=error, tolerance=tolerance,
            scaling="abs(A_JMP-A_PRL)/max(1,abs(A_JMP),abs(A_PRL))",
            measure=qtm.measure(len(u)),
            jmp=dict(amplitude=jmp["amplitude"], density=jmp["density"],
                     literal_first_column_amplitude=jmp["literal_first_column_amplitude"],
                     null_determinants=jmp["null_determinants"],
                     max_quad_error=jmp["max_quad_error"],
                     phi_cache_size=jmp["phi_cache_size"], seconds=middle-start,
                     maxdegree=maxdegree),
            prl=dict(amplitude=prl["amplitude"], density=prl["density"],
                     null_determinants=prl["null_determinants"],
                     determinant_slopes=prl["determinant_slopes"], seconds=finish-middle),
        )
        print(f"{name:17s} dps={digits}: scaled difference={mp.nstr(error, 5)} "
              f"JMP={middle-start:.2f}s PRL={finish-middle:.2f}s "
              f"{'PASS' if row['passed'] else 'FAIL'}", flush=True)
        return encode(row)


def finite_q_check(digits):
    with mp.workdps(digits):
        x, y, iu, iv = FIXTURES["ell2"]
        u, v = [mp.mpc(z, iu) for z in x], [mp.mpc(z, iv) for z in y]
        slope = qtm.prl(u, v)["determinant_slopes"][0]
        rows = []
        for text in (".08", ".04", ".02", ".01"):
            gamma = mp.mpf(text)
            ratio = mp.det(qtm.finite_q_matrix(u, v, gamma))/gamma/slope
            rows.append(dict(gamma=gamma, determinant_over_gamma_over_limit=ratio,
                             relative_error=abs(ratio-1)))
        passed = (rows[0]["relative_error"] < mp.mpf(".01")
                  and all(b["relative_error"] < a["relative_error"]*mp.mpf(".27")
                          for a, b in zip(rows, rows[1:])))
        return encode(dict(fixture="ell2", digits=digits, passed=bool(passed),
                           determinant_slope=slope, rows=rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, nargs="+", default=[25, 40])
    parser.add_argument("--fixture", choices=FIXTURES, nargs="+", default=list(FIXTURES))
    parser.add_argument("--maxdegree", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True,
                        help="new JSON result; existing files are never overwritten")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    if any(d < 20 for d in args.digits) or args.maxdegree < 4:
        parser.error("use at least 20 decimal digits and quadrature maxdegree >= 4")
    hashes = source_hashes()
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    rows = []
    for name in args.fixture:
        for digits in args.digits:
            rows.append(comparison(name, digits, args.maxdegree))
    finite = finite_q_check(max(args.digits))
    # Rebuild complex numbers from their strings at the HIGHER precision.
    # This checks stability against precision, not only two routes agreeing.
    ladder = []
    with mp.workdps(max(args.digits)):
        for name in args.fixture:
            selected = sorted((r for r in rows if r["fixture"] == name), key=lambda r: r["digits"])
            for a, b in zip(selected, selected[1:]):
                for route in ("jmp", "prl"):
                    av, bv = a[route]["amplitude"], b[route]["amplitude"]
                    av, bv = mp.mpc(av["real"], av["imag"]), mp.mpc(bv["real"], bv["imag"])
                    difference = abs(av-bv)/max(1, abs(av), abs(bv))
                    tolerance = mp.power(10, -(a["digits"]-10))
                    ladder.append(encode(dict(
                        fixture=name, route=route, from_digits=a["digits"], to_digits=b["digits"],
                        scaled_difference=difference, tolerance=tolerance,
                        passed=bool(difference < tolerance))))
    if hashes != source_hashes():
        raise RuntimeError("source changed during audit; results not written")
    result = dict(
        schema="xxx-isotropic-qtm-pointwise-v1",
        conventions="Pauli zz symmetric amplitude; density includes ell!^-2 (2*pi)^(-2*ell)",
        scope="pointwise only; no outer integration, physical spinon identification or tail estimate",
        python=platform.python_version(), mpmath=mp.__version__,
        platform=platform.platform(), source_sha256=hashes,
        invocation=sys.argv[1:], wall_seconds=time.perf_counter()-wall_start,
        cpu_seconds=time.process_time()-cpu_start, processes=1,
        passed=all(r["passed"] for r in rows+ladder+[finite]),
        pointwise=rows, precision_ladder=ladder, finite_q=finite,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(f"wrote {args.output}; overall {'PASS' if result['passed'] else 'FAIL'}", flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
