"""Compare two installed Bethe-Yang executables at identical numerical targets.

Run separately from CTest: python benchmark_sine_gordon_bethe_yang.py BEFORE AFTER
The fixed p=2, ML=10 reference is independent of the production Fourier phase.
"""
import argparse
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--precision", nargs="+", default=["fp64", "bits-128", "bits-256"])
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    programs = {name: str(getattr(args, name).resolve()) for name in ("before", "after")}
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")
    reference = Decimal(
        "2.0857356187851404880309779322375079570900382654030915686633279793820126035814552563203616916129861034106042386768200734652"
    )
    observations = []
    for precision in args.precision:
        for repeat in range(args.repeats):
            for name in (("before", "after") if repeat % 2 == 0 else ("after", "before")):
                command = [programs[name], "--p", "2", "--length", "10", "--mass", "1", "--number", "0.5",
                           "--levels", "1", "--threads", "1", "--precision", precision,
                           "--max-iterations", "512", "--format", "json"]
                started = time.monotonic()
                result = subprocess.run(command, env=environment, capture_output=True, text=True, timeout=1800)
                elapsed = time.monotonic() - started
                if result.returncode:
                    raise RuntimeError(f"{name} {precision}: {result.stderr}\n{result.stdout}")
                document = json.loads(result.stdout)
                table = document["tables"]["levels"]
                if document["status"] != "complete" or len(table["rows"]) != 1:
                    raise RuntimeError("incomplete or unexpected benchmark output")
                row = dict(zip([c["id"] for c in table["columns"]],
                               [v["value"] if isinstance(v, dict) else v for v in table["rows"][0]]))
                if row["status"] != "converged" or row["phase_status"] != "converged":
                    raise RuntimeError(f"unconverged benchmark observation: {row}")
                with localcontext() as context:
                    context.prec = 160
                    tolerance = Decimal(table["metadata"]["Counting tolerance"])
                    residual = sum(Decimal(row[k]) for k in ("residual", "quadrature_error", "tail_bound"))
                    error = abs(Decimal(row["energy"]) - reference)
                    if residual > tolerance or error > tolerance:
                        raise RuntimeError(f"numerical target missed: residual={residual}, energy error={error}")
                observations.append({
                    "executable": name, "precision": precision, "repeat": repeat + 1,
                    "iterations": int(row["iterations"]), "phase_evaluations": int(row["evaluations"]),
                    "grid_cpu_seconds": float(table["summary"]["CPU time"].removesuffix(" s")),
                    "process_cpu_seconds": float(table["summary"]["Run CPU seconds"]),
                    "process_wall_seconds": elapsed, "counting_tolerance": str(tolerance),
                    "residual_plus_phase_error": str(residual), "energy_error": str(error),
                    "energy": row["energy"], "rapidity": row["rapidity"],
                    "bethe_revision": table["metadata"]["Bethe revision"],
                    "uni20_revision": table["metadata"]["Uni20 revision"],
                    "compiler": table["metadata"]["Compiler"], "build_type": table["metadata"]["Build type"],
                })
                print(f"{name} {precision} repeat {repeat + 1}: {row['iterations']} integrals, "
                      f"{row['evaluations']} evaluations", file=sys.stderr, flush=True)
    print(json.dumps({"schema": "bethe-bethe-yang-benchmark-v1", "programs": programs,
                      "parameters": {"p": "2", "mass": "1", "length": "10", "number": "0.5",
                                     "levels": 1, "threads": 1, "max_iterations": 512},
                      "observations": observations}, indent=2))


if __name__ == "__main__":
    main()
