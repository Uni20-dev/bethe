# Reproducing a large-domain six-spinon slice

<!-- BEGIN GENERATED PAGE LINK -->
[Read this page on the documentation site](https://uni20-dev.github.io/bethe/docs/tutorials/xxx-six-spinon-reproduction/)
<!-- END GENERATED PAGE LINK -->

[Six-spinon tutorial](xxx-six-spinon.md) ·
[Selected-sector options](../xxx-spinon-sectors.md)

This optional recipe specifies sixteen unbroadened frequencies at nominal
$`q=\pi/2`$, with longitudinal normalization, $`J=1`$ and rapidity cutoff
$`R=1040`$. Each frequency uses 32 independent replicates of 768 draws.
The fixed budget was chosen from the [calibration calculations](../validation/six-spinon-slice-calibration-20261031.json)
before production; it does not guarantee a particular sampling error.

**Production is still in progress.** This recipe records the requested grid
and controls; it is not a completed-data or convergence claim. The grid spans
$`1.6\le\omega/J\le4`$, including both sides of the two-spinon upper edge.
It does not reach the four- or six-spinon outer upper edges, and it is not a
grid for evaluating global moments.

## Run the installed program

Use an MPFR-enabled [installation](../building.md). The command below expects
`bethe-xxx-structure-factor-thermo` on `PATH`; from a build directory, replace
that name with `./bethe-xxx-structure-factor-thermo`.

This is a costly calculation, unlike the small example in the tutorial.
Pilot timings suggest about 1.5–2.2 hours **per frequency with sixteen workers**
on the measured Xeon Gold 6254 system; rare difficult observations and machine
load can change that cost. The recipe defaults to four workers. Set `threads`
to match the CPUs available to your job, and choose `out` on a filesystem with
room for the exports and logs. The source tree need not be your working directory.

```sh
(
  set -eu
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  threads=4
  out=six-slice-raw
  mkdir "$out"
  q=1.5707963267948966192313216916397514420985846996875529104874722961539082031431045
  cat > "$out/grid.txt" <<'GRID'
1.6 7400
1.65 7401
1.7 7402
1.8 7403
1.95 7404
2 7405
2.1 7406
2.18 7407
2.215 7408
2.23 7409
2.3 7410
2.5 7411
2.75 7412
3 7413
3.5 7414
4 7415
GRID
  while read -r omega seed; do
    bethe-xxx-structure-factor-thermo --sectors 6 \
      --momentum "$q" --omega "$omega" --exchange 1 --channel zz \
      --precision bits-1536 --density-guard-bits 320 \
      --six-cutoff 1040 --replicates 32 --samples 768 --seed "$seed" \
      --uniform-bits 52 --six-soft-fraction 0.5 \
      --six-background-weight 0.5 --six-one-pair-weight 0.25 \
      --six-two-pairs-weight 0.25 --six-one-pair-scale 4 \
      --six-one-pair-gap 0.5 --six-two-pairs-scale 8 --six-two-pairs-gap 1 \
      --density-terms 1560 --density-tail-order 48 --density-seconds 60 \
      --density-relative-tolerance 1e-6 --density-absolute-tolerance 0 \
      --threads "$threads" --diagnostics --replicate-output --quiet \
      --json "$out/w${omega}.json" \
      > "$out/w${omega}.stdout" 2> "$out/w${omega}.stderr"
  done < "$out/grid.txt"
)
```

The new directory prevents accidental overwriting of a previous run. A nonzero
program exit stops the loop and retains that query's output and logs. Preserve
failed results and investigate their diagnostics before rerunning the same
frequency, seed and sample count in a different directory. Do not replace a
failure by zero or average only successful observations.

The recipe states the proposal parameters explicitly, including defaults, so
they remain visible when comparing implementations. Geometry uses 1536 bits;
the density comparisons use 256/320 bits. The contour comparisons use
1560/1576 residue terms and tail orders 48/52. These controls address distinct
numerical limitations, as explained in the [precision study](xxx-six-spinon.md#distinguish-contour-failures-from-arithmetic-failures).
Sampling error still excludes the omitted rapidity domain.

## Pack and redraw completed results

Once every frequency has completed, run this from the source directory, with
`raw_dir` pointing to the calculation directory above. It reads the recorded
grid, calls the shared validator and packs all requested frequencies. Failed,
missing or inconsistent exports are errors, so a partial run cannot silently
become a shorter successful curve.

```python
from pathlib import Path
import subprocess
import sys

raw_dir = Path("six-slice-raw")  # Change this if the calculation ran elsewhere.
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

This produces light and dark SVG figures. The compact data retain decimal
values, source revisions, controls, costs, sampling diagnostics, replicate
summaries and raw-export hashes. They suffice to redraw and independently
check the plotted estimates without distributing per-observation payloads.
Different program versions or machines need not reproduce every last digit;
compare the physical estimates, numerical checks and observed sampling errors.
