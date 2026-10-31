# Reproducing the two-, four- and six-spinon comparison

<!-- BEGIN GENERATED PAGE LINK -->
[Read this page on the documentation site](https://uni20-dev.github.io/bethe/docs/tutorials/xxx-spinon-slice-reproduction/)
<!-- END GENERATED PAGE LINK -->

[Combined comparison](xxx-spinon-slice.md) ·
[Six-spinon reproduction](xxx-six-spinon-reproduction.md) ·
[Selected-sector options](../xxx-spinon-sectors.md)

This optional calculation produces two- and four-spinon spectra on the same
sixteen-frequency grid as the six-spinon recipe. We use $`q=\pi/2`$, $`J=1`$
and longitudinal normalization throughout. No normalization multiplier is
applied when the sectors are combined.

**All sixteen two/four queries are complete; six-spinon production is still
in progress.** The recipe records the requested grid and controls. Completing
these queries does not bound the omitted integration domains. The frequency range
$`1.6\le\omega/J\le4`$ samples both sides of the two-spinon upper edge. It
omits the outer four- and six-spinon upper edges and is not a moment grid.

## Calculate the two- and four-spinon sectors

Use an MPFR-enabled [installation](../building.md), with
`bethe-xxx-structure-factor-thermo` on `PATH`. From a build directory, use
`./bethe-xxx-structure-factor-thermo` instead. The source directory need not be
your working directory. Set `threads` to match the CPUs available to the job
and choose a fresh output directory on a filesystem with room for the results.

These are substantial calculations. For example, the $`\omega/J=1.8`$ query
took about 24 minutes with four workers on the measured Xeon Gold 6254 system.
The [calibration record](../validation/four-spinon-slice-calibration-20261031.json)
retains timings for the representative controls; cost varies with the domain,
order, hardware and load.

The columns below specify frequency, quadrature order, log-radius cutoff,
geometry bits, density guard bits, residue terms and tail order. More residue
terms and wider geometry accompany the deep soft domain below the two-spinon
upper edge. Above it, the selected domain is smaller and order refinement
determines the cost. The two sets of controls follow the
[separate refinement checks](xxx-spinon-slice.md#measured-four-spinon-controls).

The separate `--tolerance 1e-40` below applies only to the two-spinon kernel.
It preserves the setting used for these production exports; it does not request
or demonstrate forty-digit accuracy for the four-spinon integral. For an ordinary
two-spinon point, the [introductory command](xxx-structure-factor-thermo.md#generate-the-data)
uses validated binary64 arithmetic and a `1e-12` kernel tolerance. The wider
geometry here is needed for the four-spinon soft domain, as explained in the
calibration above.

```sh
(
  set -eu
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  threads=4
  out=four-slice-raw
  mkdir "$out"
  q=1.5707963267948966192313216916397514420985846996875529104874722961539082031431045
  cat > "$out/grid.txt" <<'GRID'
1.6 8 768 1536 256 1536 48
1.65 8 768 1536 256 1536 48
1.7 8 768 1536 256 1536 48
1.8 8 768 1536 256 1536 48
1.95 8 768 1536 256 1536 48
2 8 768 1536 256 1536 48
2.1 8 768 1536 256 1536 48
2.18 8 768 1536 256 1536 48
2.215 8 768 1536 256 1536 48
2.23 32 24 192 192 64 16
2.3 32 24 192 192 64 16
2.5 64 24 192 192 64 16
2.75 64 24 192 192 64 16
3 64 24 192 192 64 16
3.5 64 24 192 192 64 16
4 64 24 192 192 64 16
GRID
  while read -r omega order cutoff bits guard terms tail; do
    bethe-xxx-structure-factor-thermo --sectors 2,4 \
      --momentum "$q" --omega "$omega" --exchange 1 --channel zz \
      --precision "bits-$bits" --density-guard-bits "$guard" \
      --four-order "$order" --four-cutoff "$cutoff" \
      --density-terms "$terms" --density-tail-order "$tail" \
      --density-relative-tolerance 1e-8 --density-absolute-tolerance 1e-30 \
      --tolerance 1e-40 --threads "$threads" --diagnostics --quiet \
      --json "$out/w${omega}.json" \
      > "$out/w${omega}.stdout" 2> "$out/w${omega}.stderr"
  done < "$out/grid.txt"
)
```

The new directory prevents overwriting a previous run. A failed query stops
the loop and leaves its output and logs for inspection. Keep it, investigate
the failure, and rerun the complete query in another directory. Do not drop
the frequency or replace a failed sector by zero.

## Recorded production cost

The [production record](../validation/four-spinon-slice-production-20261031.json)
retains each frequency's command, raw-export hash, executable hash, CPU model,
worker count and timings. The sixteen selected two/four queries used
52,478.8 CPU seconds, or about 14.6 CPU hours. This excludes the calibration
studies and all six-spinon work. Thirteen queries ran on an Intel Core
i7-12700KF and three on a Xeon Gold 6254; their differing queries and controls
do not constitute a hardware benchmark.

The measurements used four workers except at $`\omega/J=2.3,3,4`$, where
two were used. The recipe above uses four consistently. CPU time sums work
across threads; wall time measures how long a query takes. The record separates
application CPU and elapsed time from the external driver's wall time. Its
wall-time totals sum individual queries, some of which overlapped, so they
are not the elapsed duration of the campaign. The source revisions and
executable hashes identify the original calculations; later code improvements
do not alter these historical measurements.

## Calculate six spinons and combine the results

Run the [six-spinon recipe](xxx-six-spinon-reproduction.md), which supplies
the same declared frequencies and produces the checked compact file
`six-slice.json`. Use the same Bethe and Uni20 revisions for both sectors:
the joining script rejects mixed source revisions or physical conventions.

Once every requested frequency has completed, run the following from the
source directory, adjusting `raw_dir` and `six_data` to your output locations.
It passes the entire declared grid to the validator. Missing, failed or
inconsistent two/four exports stop the import; failed or unavailable six-spinons
in an explicitly incomplete compact record remain visible gaps.

```python
from pathlib import Path
import subprocess
import sys

raw_dir = Path("four-slice-raw")
six_data = Path("six-slice.json")
nodes = [line.split()[0] for line in (raw_dir / "grid.txt").read_text().splitlines()]
command = [sys.executable, "scripts/plot_xxx_spinon_slice_tutorial.py",
           "--four-json", *(str(raw_dir / f"w{omega}.json") for omega in nodes),
           "--six-data", str(six_data), "--data", "spinon-slice.json",
           "--figure", "spinon-slice.svg"]
subprocess.run(command, check=True)
```

The plotting dependencies are listed in `docs/tutorials/requirements.txt`.
The command writes a compact record and light/dark SVG figures. Later redraws
use only that record:

```sh
python3 scripts/plot_xxx_spinon_slice_tutorial.py \
  --data spinon-slice.json --figure spinon-slice.svg
```

Both the individual six-spinon estimates and the $`2+4+6`$ partial sum carry
the six-spinon sampling standard error. Numerical and omitted-domain effects
are recorded separately, so these bars do not represent total uncertainty.
No missing value is filled with zero, and no sector is rescaled.
