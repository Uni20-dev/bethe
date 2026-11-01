# Reproducing the six-spinon reference dataset

<!-- BEGIN GENERATED PAGE LINK -->
[Read this page on the documentation site](https://uni20-dev.github.io/bethe/docs/tutorials/xxx-six-spinon-reproduction/)
<!-- END GENERATED PAGE LINK -->

[Six-spinon tutorial](xxx-six-spinon.md) ·
[Selected-sector options](../xxx-spinon-sectors.md)

This page reproduces the sixteen-frequency reference slice at $`q=\pi/2`$,
$`J=1`$, with longitudinal normalization and rapidity cutoff $`R=1040`$.
The points are unbroadened. They span $`1.6\le\omega/J\le4`$, including both
sides of the two-spinon upper edge. This display grid does not reach the outer
four- or six-spinon edges and is not a grid for calculating global moments.

For a hands-on introduction, use the [small examples](xxx-six-spinon.md) first.
You can redraw the reference figure immediately from the
[compact dataset](data/xxx-six-spinon-slice.json); rerunning the integrations is
optional. The compact file retains the measured values, replicate summaries,
precision stages, numerical controls, costs and original source identities.

## Dataset settings and compute budget

The sixteen queries used 30.34 CPU hours and 2.49 wall
hours in total on a Xeon Gold 6254 with 24 workers. Individual queries took
3.4--24.3 minutes. These are measured costs of the recorded
executable, not a runtime guarantee for another machine. Start with one query
when budgeting a new calculation. The [validation record](../validation/six-spinon-reference-20261102.json)
retains per-point costs, source revisions, checks and sampling diagnostics.

Each frequency has 32 independent replicates. The sample counts below were
chosen before this campaign, using earlier sampling errors where available.
They total 1,613,824 draws. The observed relative sampling errors are
2.49--6.52%; 15 of the sixteen points meet the initial 5% sampling target. These errors describe observed replicate
scatter, not a bound on rare unsampled configurations or the omitted rapidity
domain. The [tutorial](xxx-six-spinon.md#read-the-errors-before-increasing-the-grid)
explains these separate uncertainties.

Most observations use native arithmetic. On the recorded x86 build, the
ordered stage list `fp80,fp128,reference` first attempts extended80 geometry
with a 53/64-bit density comparison, then binary128 geometry with a 64/113-bit
comparison. A final MPFR retry uses 1536-bit geometry and a 256/320-bit density
comparison. Accumulation and output use 1536 bits throughout. Thus the large
reference precision does not imply that every form factor is evaluated in MPFR.

Every stage uses the same contour controls and target tolerances: 3120/3136
residue terms, tail orders 52/56, relative tolerance `1e-6`, and zero absolute
tolerance. A retry repeats the same random draw, including the geometry and
all supported roots. It does not discard a difficult observation or draw a
replacement. See [precision retries](xxx-six-spinon.md#retry-difficult-configurations-at-wider-precision)
for the meaning of the individual checks.

## Run the installed program

Use an [installation](../building.md) with extended80, binary128 and MPFR enabled
for this stage list. The executable is assumed to be on `PATH`; from a build
directory, use `./bethe-xxx-structure-factor-thermo`.

The recipe defaults to four workers. Set `threads` to the CPUs available to your
job, and place `out` on a filesystem with room for the exports and logs. The
source tree need not be your working directory. `--samples-per-task 128` is the
measured choice for the reference machine; scheduling grain does not change
the sample identities or accumulation order.

```sh
(
  set -eu
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  threads=4
  out=six-slice-raw
  mkdir "$out"
  q=1.5707963267948966192313216916397514420985846996875529104874722961539082031431045
  cat > "$out/grid.txt" <<'GRID'
1.6 7400 7936
1.65 7401 2304
1.7 7402 3072
1.8 7403 2048
1.95 7404 1280
2 7405 1792
2.1 7406 1792
2.18 7407 3584
2.215 7408 1024
2.23 7409 1280
2.3 7410 3840
2.5 7411 4096
2.75 7412 4096
3 7413 4096
3.5 7414 4096
4 7415 4096
GRID
  while read -r omega seed samples; do
    bethe-xxx-structure-factor-thermo --sectors 6 \
      --momentum "$q" --omega "$omega" --exchange 1 --channel zz \
      --precision bits-1536 --density-guard-bits 320 \
      --six-precision-stages fp80,fp128,reference --samples-per-task 128 \
      --six-cutoff 1040 --replicates 32 --samples "$samples" --seed "$seed" \
      --uniform-bits 52 --six-soft-fraction 0.5 \
      --six-background-weight 0.5 --six-one-pair-weight 0.25 \
      --six-two-pairs-weight 0.25 --six-one-pair-scale 4 \
      --six-one-pair-gap 0.5 --six-two-pairs-scale 8 --six-two-pairs-gap 1 \
      --density-terms 3120 --density-tail-order 52 --density-seconds 60 \
      --density-relative-tolerance 1e-6 --density-absolute-tolerance 0 \
      --detail-row-limit 1000000 --threads "$threads" \
      --diagnostics --replicate-output --quiet \
      --json "$out/w${omega}.json" \
      > "$out/w${omega}.stdout" 2> "$out/w${omega}.stderr"
  done < "$out/grid.txt"
)
```

The new directory protects previous output. A nonzero exit stops the loop and
retains that query's logs. Preserve failed results and investigate the reason
before rerunning the same input in a different directory. Neither a zero
placeholder nor an average over successful survivors represents a failed query.
The detail-row limit permits the worst-case recovery accounting at the largest
sample count; it does not change numerical precision or tolerances.

The proposal weights and scales are explicit, including their defaults. This
keeps the importance distribution visible when comparing implementations.
Raising the cutoff, changing that distribution or increasing the sample count
is a new numerical experiment; it is not required to redraw this dataset.

## Pack and redraw completed results

Once all sixteen queries complete, run the following from the source directory,
with `raw_dir` pointing to their output. It reads the declared grid and validates
every replicate and recovery record. Missing, failed or inconsistent exports
are errors; an incomplete run cannot silently become a shorter successful curve.

```python
from pathlib import Path
import subprocess
import sys

raw_dir = Path("six-slice-raw")
q = "1.5707963267948966192313216916397514420985846996875529104874722961539082031431045"
nodes = [line.split()[0] for line in (raw_dir / "grid.txt").read_text().splitlines()]
command = [sys.executable, "scripts/xxx_six_spinon_tutorial_data.py",
           "--momentum", q, "--cutoff", "1040", "--output", "six-slice.json"]
for omega in nodes:
    command += ["--frequency", omega, "--raw", str(raw_dir / f"w{omega}.json")]
subprocess.run(command, check=True)
```

With the optional plotting dependencies installed:

```sh
python3 -m pip install -r docs/tutorials/requirements.txt
python3 scripts/plot_xxx_six_spinon_tutorial.py \
  --compact six-slice.json \
  --momentum 1.5707963267948966192313216916397514420985846996875529104874722961539082031431045 \
  --figure six-slice.svg
```

This produces light and dark SVG figures. To redraw the supplied results, use
`--compact docs/tutorials/data/xxx-six-spinon-slice.json` instead. The compact
file contains decimal values, sampling and recovery diagnostics, controls and
raw-export hashes. Per-observation payloads and raw archives are not needed for
the plot. Different program versions or machines need not reproduce every last
digit; compare the physical estimates, numerical checks and sampling errors.
