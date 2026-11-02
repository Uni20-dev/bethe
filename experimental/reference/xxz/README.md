# Independent finite-q vertex reference

Requires Python 3 and mpmath; the fixed-grid controls also require NumPy and
SciPy. These are optional reference dependencies, not C++ build dependencies. This directory runs on its own; it imports neither
the C++ implementation nor research data. From the repository root:

```sh
python3 -m pip install -r experimental/reference/xxz/requirements.txt
python3 -m unittest discover -s experimental/reference/xxz -p 'test_*.py' -v
python3 experimental/reference/xxz/vertex_reference.py --epsilon .4 --dps 65 --terms 96
python3 experimental/reference/xxz/vertex_quadrature.py moments --spinons 4 --order 4
python3 experimental/reference/xxz/vertex_quadrature.py cut --spinons 4 --order 4
```

`Vertex(epsilon)` uses epsilon=acosh(Δ). Its rapidity argument β is related to
the C++ spectral angle by φ=epsilon β/π. `amplitudes(beta, vacuum, terms)` returns
the raw complex two/four-spinon amplitudes before completeness, Fourier and
energy Jacobians. It evaluates literal q-products and independently shifted
residues; it does not call the C++ cancelled recurrence.

`vertex-point.json` records all eight four-spinon amplitudes at Δ=1.5, with
65/100-digit and 40/64-level refinement. This is a pointwise control, not a
converged structure factor or an outer-integral fixture. Complex phase agreement
is stronger than checking only a squared density.

The command-line scan prints raw rapidity densities for examining the approach
to the isotropic limit. Neither those values nor a finite-order contour sum
prove normalization or convergence at arbitrary rapidities. Coincident nodes
and confluent limits require separate analysis.

See [the API guide](../../../docs/experimental-xxz-vertex.md) and
[Abada–Bougourzi–Si-Lakhal](../../../CITATIONS.md#abada-bougourzi-si-lakhal-1997),
Eqs. (2.6)–(2.7).

`vertex_quadrature.py` uses SciPy elliptic functions and vectorized NumPy
amplitudes to reproduce the fixed-order C++ checks. `vertex_numpy.py` provides
literal q-products and the cancelled residue recurrence; this recurrence is
shared mathematically with C++, so it is not an independent derivation. The
mpmath literal-residue tests provide that additional control. The quadrature
checks cover winding/sector accounting, full-zone moments, fixed-k cuts and
the full-axis Gaussian normalization. The independent periodic sech sum checks
the spinon energy.

These binary64 reference modules deliberately reject unresolved elliptic
parameters and momentum shifts. The C++ evaluator retains the complementary
modulus directly and can use wider precision. Neither implementation discards
failed amplitude nodes or renormalizes a finite frequency window. All tests
use small explicit orders; agreement at one order does not certify a
thermodynamic integral or its product/residue tails.
