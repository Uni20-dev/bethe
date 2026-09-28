#!/usr/bin/env python3
"""Validate and plot the tutorial's actual frontend exports.

With --solver PATH, first regenerate CSVs (only after all three runs succeed).
Without it, use the checked-in data. No Bethe formulas are used to draw curves.
"""

import argparse
import math
from pathlib import Path
import subprocess

from tutorial_common import solver_executable, finite_number, read_csv_export, save_svg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/tutorials/data"
FIGURES = ROOT / "docs/tutorials/figures"
POINTS = 401
CASES = {"xxz-delta1": (1, False), "xxz-delta2": (2, False),
         "xxz-delta2-folded": (2, True)}
COLUMNS = ["branch", "p", "p_over_pi", "cell_momentum", "energy", "lower", "upper", "status"]


def read_export(text, delta, folded):
    expected = {"Program": "bethe-xxz-dispersion", "Delta": str(delta),
                "Exchange J": "1", "Precision": "fp64", "Branches": "all",
                "Points per branch": str(POINTS), "Bulk status": "converged",
                "Spinon status": "converged", "Outcome": "success", "Status": "converged",
                "Continuum momentum": ("Q modulo pi; union of both translation branches"
                                       if folded else "p1+p2=Q modulo 2*pi; unfolded")}
    _, rows = read_csv_export(text, COLUMNS, expected)
    branches = {"spinon": [], "two-spinon": []}
    for row in rows:
        branch = row["branch"]
        if branch not in branches or row["status"] != "converged" or None in row:
            raise ValueError(f"Invalid/incomplete row: {row}")
        required = ("energy",) if branch == "spinon" else ("lower", "upper")
        absent = ("lower", "upper") if branch == "spinon" else ("energy",)
        if any(row[key] != "" for key in absent):
            raise ValueError("Non-applicable energy columns must be empty")
        numbers = {key: finite_number(row, key) for key in ("p", "p_over_pi", "cell_momentum", *required)}
        if branch == "two-spinon" and not 0 <= numbers["lower"] <= numbers["upper"]:
            raise ValueError("Invalid continuum bounds")
        if branch == "spinon" and numbers["energy"] < 0:
            raise ValueError("Negative spinon energy")
        branches[branch].append(numbers)
    for branch, rows in branches.items():
        if len(rows) != POINTS:
            raise ValueError(f"{branch}: expected {POINTS} rows, got {len(rows)}")
        extent = 1 if branch == "spinon" else 2
        for i, row in enumerate(rows):
            if not math.isclose(row["p_over_pi"], extent * i / (POINTS - 1), abs_tol=1e-14):
                raise ValueError("Incomplete/out-of-order momentum grid")
            if not math.isclose(row["p"], math.pi * row["p_over_pi"], abs_tol=1e-14):
                raise ValueError("Momentum units disagree")
    return branches


def load_cases(solver=None):
    exports = {}
    for name, (delta, folded) in CASES.items():
        if solver:
            command = [str(solver), "--delta", str(delta), "--exchange", "1",
                       "--branch", "all", "--precision", "fp64", "--points", str(POINTS),
                       "--format", "csv"]
            if folded:
                command.append("--folded")
            text = subprocess.run(command, check=True, text=True, capture_output=True).stdout
        else:
            text = (DATA / f"{name}.csv").read_text()
        exports[name] = (text, read_export(text, delta, folded))
    if solver:
        DATA.mkdir(parents=True, exist_ok=True)
        for name, (text, _) in exports.items():
            (DATA / f"{name}.csv").write_text(text)
    return {name: rows for name, (_, rows) in exports.items()}


def plot(cases):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"svg.hashsalt": "bethe-xxz-tutorial", "font.size": 11,
                         "axes.spines.top": False, "axes.spines.right": False})
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = ["#0072B2", "#D55E00"]

    def continuum(ax, rows, color, label):
        x = [row["p_over_pi"] for row in rows]
        low = [row["lower"] for row in rows]
        high = [row["upper"] for row in rows]
        ax.fill_between(x, low, high, color=color, alpha=.2, label=label)
        ax.plot(x, low, color=color, linewidth=1.5)
        ax.plot(x, high, color=color, linewidth=1.5)
        ax.set(xlabel=r"Total constituent momentum $Q/\pi$", ylabel=r"$E/J$",
               xlim=(0, 2), ylim=(0, 4.5), xticks=[0, .5, 1, 1.5, 2])
        ax.grid(alpha=.15)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for col, delta in enumerate((1, 2)):
        rows = cases[f"xxz-delta{delta}"]
        ax = axes[0, col]
        ax.plot([r["p_over_pi"] for r in rows["spinon"]],
                [r["energy"] for r in rows["spinon"]], color=colors[col], linewidth=2)
        ax.set(title=rf"$\Delta={delta}$ · single spinon", xlabel=r"Spinon momentum $p/\pi$",
               ylabel=r"$\epsilon/J$", xlim=(0, 1), ylim=(0, 2.3), xticks=[0, .25, .5, .75, 1])
        ax.grid(alpha=.15)
        continuum(axes[1, col], rows["two-spinon"], colors[col], "Kinematically allowed")
        axes[1, col].set_title(rf"$\Delta={delta}$ · two spinons, unfolded")
    save_svg(fig, FIGURES / "xxz-spinons.svg")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for ax, name, label in zip(axes, ("xxz-delta2", "xxz-delta2-folded"),
                              ("Unfolded constituent branch", "Folded union (two-site cell)")):
        continuum(ax, cases[name]["two-spinon"], colors[1], label)
        ax.set_title(label)
        ax.axvline(1, color="0.4", linestyle=":", linewidth=1)
    axes[1].set_xlabel(r"Representative $Q/\pi$ (repeats every 1)")
    save_svg(fig, FIGURES / "xxz-folding.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=solver_executable, help="Regenerate data with this executable")
    parser.add_argument("--check", action="store_true", help="Validate exports without plotting")
    args = parser.parse_args()
    cases = load_cases(args.solver.resolve() if args.solver else None)
    if not args.check:
        plot(cases)
    print(f"Validated {len(cases)} exports ({2 * POINTS} rows each).")


if __name__ == "__main__":
    main()
