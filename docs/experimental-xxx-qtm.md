# Isotropic QTM reference amplitudes

[Experimental components](experimental.md) ·
[Derivation](experimental-xxx-qtm-derivation.md) ·
[XXX spinon densities](xxx-spinon-density.md)

The optional Python package in `experimental/reference/qtm/` compares two
representations of the XXX quantum-transfer-matrix (QTM) particle–hole
amplitude at the **same complex rapidities**. It then compares their full
partition sum at the real boundary with the public two-, four- and six-spinon
vertex references. This separates special-function and determinant errors
from errors in an outer integration.

| Route | Calculation | Source |
| --- | --- | --- |
| `jmp` | Ordinary Gamma functions inside one-dimensional real-line quadratures | [Babenko et al., JMP 2021](../CITATIONS.md#babenko-gohmann-kozlowski-suzuki-2021) |
| `prl` | Analytic isotropic limit of the basic-hypergeometric representation; ordinary hypergeometric functions at unit argument | [Babenko et al., PRL 2021](../CITATIONS.md#babenko-gohmann-kozlowski-sirker-suzuki-2021) |

The two routes share mpmath Gamma/Barnes functions and physical conventions.
Their quadrature and hypergeometric evaluations provide different numerical
controls, rather than fully independent derivations of physical normalization.
The partition comparison follows the correspondence proposed by
[Dugave et al. (2016)](../CITATIONS.md#dugave-gohmann-kozlowski-suzuki-2016).
General higher-sector equivalence and a global contour deformation are not
proved by these pointwise checks.

## Run a small comparison

Use Python 3.10 or later and mpmath. No Bethe C++ build or CMake option is
required. These scripts are available in the source checkout and are not
installed as a Python package. From the repository root:

```sh
python3 -m venv /path/to/qtm-venv
/path/to/qtm-venv/bin/python -m pip install -r experimental/reference/qtm/requirements.txt
/path/to/qtm-venv/bin/python experimental/reference/qtm/audit_xxx_qtm_partitions.py \
  --fixture two --output /path/to/scratch/qtm-partitions.json
```

The output path must be new. This command generates its own boundary reference,
checks both labelled partitions at two interior contour offsets and two
precisions, and compares the full complex results with the vertex amplitude.
It uses no saved research data. A nonzero exit status reports a failed check;
the output retains the numerical diagnostics when the audit completes.

For direct use, construct decimal inputs **inside** the requested precision:

```python
import sys
sys.path.insert(0, "experimental/reference/qtm")
import mpmath as mp
import reference_xxx_qtm as qtm

with mp.workdps(35):
    holes = [mp.mpc("-.7", "-.4")]
    particles = [mp.mpc(".3", ".6")]
    a = qtm.jmp(holes, particles)
    b = qtm.prl(holes, particles)
    error = abs(a["amplitude"] - b["amplitude"])
    print(mp.nstr(a["amplitude"], 30), mp.nstr(error, 5))
```

The caller owns the mpmath context. Passing a Python binary64 value cannot
recover digits already lost before conversion. The package does not alter
the process-wide default precision on import. Results retain complex values;
taking an absolute value or real part before the partition sum changes the
calculation.

## Operator and integration measures

For $`\ell`$ holes $`U`$ and $`\ell`$ particles $`V`$, both routes return the
symmetric **Pauli-zz amplitude** $`A_\ell(U,V)`$, without a space/time phase.
The separate `density` result includes its integration measure:

```math
\rho_\ell(U,V)=\frac{A_\ell(U,V)}{(\ell!)^2(2\pi)^{2\ell}}.
```

For spin operators $`S^z=\sigma^z/2`$, divide this correlator contribution by
four. A generally complex contour amplitude is not a positive spectral weight.
The package exposes these conventions directly; it never fits a multiplier to
a measured spinon fraction or sum rule.

For real vertex rapidities $`\beta_j`$, let $`n=2\ell`$ and
$`t_j=\beta_j/\pi+i\eta`$, with $`0\lt\eta\lt1`$. For each size-$`\ell`$
subset $`H`$ and its complement $`P`$, use

```math
U_H=\{t_j-i:j\in H\},\qquad V_P=\{t_j:j\in P\}.
```

There are **2, 6 and 20 ordered assignments** for two, four and six spinons.
Complementary assignments are distinct terms, even when their boundary
coefficients become conjugates. Discarding one member loses its contribution.

The public vertex reference returns $`W_n`$ in the convention

```math
S_n^{zz}(q,\omega)=\int d^np\,W_n(\boldsymbol p)\,
\delta_{2\pi}\!\left(q-\sum_jp_j\right)
\delta\!\left(\omega-\sum_j\epsilon_j\right).
```

Both conservation deltas are still unintegrated and the periodic delta
integrates to one over a momentum period. Inverse Fourier measures contribute
$`1/(2\pi)^2`$, spin-to-Pauli conversion contributes four, and
$`dp/d\beta=1/\cosh\beta`$. Summing all QTM assignments changes its common
measure to $`1/[n!(2\pi)^n]`$; then $`dt=d\beta/\pi`$. Thus the absolute
boundary comparison is

```math
\mathcal A_n(\boldsymbol\beta)=\sum_{|H|=\ell}A_H(0)
=\frac{n!(2\pi^2)^n}{\pi^2}\,
W_n(\boldsymbol p)\prod_j\frac1{\cosh\beta_j}.
```

The factorial, operator factors and Jacobians are tested separately.
The phase mapping gives $`p(\beta)=2\arctan(e^\beta)`$ and the physical
spin-Hamiltonian energy $`\epsilon(\beta)=\pi/[2\cosh\beta]`$ at $`J=1`$.

## Contours and numerical limits

The published isotropic particle/hole contours use distinct infinitesimal
prescriptions:

```math
\widehat C_p=\mathbb R+i(1+0_+)/2,\qquad
\widehat C_h=\mathbb R-i(1+0_-)/2.
```

The reference permits nonsingular points with
$`-1\lt\Im u\lt0\lt\Im v\lt1`$. This wider test domain does not certify
deforming an entire outer contour. The fixture `symmetric-offset` uses
imaginary parts −0.55 and +0.55; it is a finite-offset test, not the published
infinitesimal prescription.

Exact same-species collisions, a hole/particle separation of $`-i`$, and the
PRL total-rapidity denominator zero need separate limiting formulas and are
rejected. Increasing arithmetic precision does not supply a missing limit.
Near a supported but ill-conditioned point, compare precision and quadrature
settings instead of trusting one successful evaluation.

The JMP result reports mpmath error estimates for matrix-entry quadratures.
They are not certified bounds on the determinant product. The audits report
cancellation indicators, determinant zeros and refinement differences; none
is an interval enclosure or an error estimate for an unevaluated outer integral.

## Choose an audit

All commands below run from a complete public checkout. Each output path must
be new; change the example paths to your writable scratch directory.

```sh
# Compare identical complex rapidities and the unexpanded finite-q determinant.
python3 experimental/reference/qtm/audit_xxx_qtm_pointwise.py \
  --fixture ell1 ell2 --digits 25 40 --output /path/to/scratch/qtm-pointwise.json

# Refine contour offset, arithmetic precision and the vertex residue tail.
python3 experimental/reference/qtm/audit_xxx_qtm_vertex.py \
  --fixture two four --output /path/to/scratch/qtm-vertex.json

# Compare every complex partition and remove the common offset analytically.
python3 experimental/reference/qtm/audit_xxx_qtm_partitions.py \
  --fixture two four --output /path/to/scratch/qtm-partitions-new.json

# Optionally reuse a passed vertex audit produced by these exact sources.
python3 experimental/reference/qtm/audit_xxx_qtm_partitions.py \
  --fixture two four --checkpoint /path/to/scratch/qtm-vertex.json \
  --output /path/to/scratch/qtm-partitions-reused.json

python3 -m unittest discover -s experimental/reference/qtm -p 'test_reference_*.py' -v
```

Activate the virtual environment, or replace `python3` with its interpreter.
`--help` lists all fixture and precision controls. `six-bulk` and
`six-asymmetric` exercise the twenty assignments; `six-clustered` probes a
more difficult neighborhood in the vertex/partition audits.

The partition audit compares two offsets, precision settings and quadrature
splittings. It forms $`A_H(0)=K_0B_H`$ by removing the explicit sine factor
algebraically, **without dividing by a small sine**. The
[derivation](experimental-xxx-qtm-derivation.md) explains why $`B_H`$ is
unchanged by the common imaginary translation. Offset extrapolation remains
a separate check in the vertex audit, not a required normalization step.

JSON results use decimal strings for arbitrary-precision reals and explicit
real/imaginary objects for complex values. They retain input rapidities,
controls, source hashes, timings and all check outcomes. A reused checkpoint
must be passed, contain the requested fixture names and exact decimal rapidity
lists, and match the numerical sources. Its recorded precision and vertex-tail
controls may differ; the new audit compares the saved complex terms with its
fresh refined calculation. Source or input-file changes during execution
reject the result.

The optional full test suite takes about two minutes on the recorded reference
CPU. It is not added to routine C++ CI. Full six-spinons-by-partition quadrature
audits can take substantially longer; no large calculation is needed to run
the small examples. A compact [reference fixture](../experimental/reference/qtm/pointwise-fixture.json)
and [validation record](validation/experimental-xxx-qtm.json) identify the
shipped controls. To regenerate the fixture in a new scratch file:

```sh
python3 experimental/reference/qtm/generate_pointwise_fixture.py \
  --digits 32 --refined-digits 42 --output /path/to/scratch/qtm-fixture.json
```

This checks all 2, 6 and 20 PRL partitions at both precisions, the full boundary
sums against the vertex target, and one independent JMP partition per sector.
It takes about three minutes on the recorded CPU. These checks support the
stated pointwise comparisons, not
general sector equivalence, a complete dynamical spectrum or a production
six-dimensional quadrature.
