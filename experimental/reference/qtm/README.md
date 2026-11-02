# Isotropic quantum-transfer-matrix references

Optional pointwise Python controls for two representations of XXX
particle–hole amplitudes and their spinon boundary correspondence.
Start with [the usage and normalization guide](../../../docs/experimental-xxx-qtm.md).
The [derivation](../../../docs/experimental-xxx-qtm-derivation.md) gives the
distinguished-column sum and nonuniform hypergeometric tail limit.

Run from a complete public Bethe checkout; the vertex comparison imports the
already-public `scripts/reference_*.py` controls from this checkout. No local
papers, private checkpoint, C++ build or cluster account is needed.
These controls take precedence over `PYTHONPATH`. If numerical reference
modules from another checkout are already imported, the audit fails explicitly;
start a fresh Python process instead of mixing implementations. Source records
include the loaded transitive controls and reject foreign reference modules.

```sh
python3 -m venv /path/to/qtm-venv
/path/to/qtm-venv/bin/python -m pip install -r experimental/reference/qtm/requirements.txt
/path/to/qtm-venv/bin/python -m unittest discover -s experimental/reference/qtm -p 'test_reference_*.py' -v
```

The full optional suite takes about two minutes on the recorded reference CPU;
it is not part of routine C++ CI. Each audit uses one process and writes a new
JSON file, preserving decimal precision, source hashes, controls and failures.
Use a writable scratch directory for outputs. Detailed multi-partition audits
can take longer than the test suite.

| File | Purpose |
| --- | --- |
| `reference_xxx_qtm.py` | JMP Gamma quadrature and analytic isotropic PRL hypergeometric limit |
| `reference_xxx_qtm_vertex.py` | All labelled partitions, physical measure conversion and vertex target |
| `audit_xxx_qtm_pointwise.py` | Precision ladder at identical complex rapidities and unexpanded finite-q control |
| `audit_xxx_qtm_vertex.py` | Boundary-offset, arithmetic and vertex-tail refinement |
| `audit_xxx_qtm_partitions.py` | Full complex partition comparison and analytic removal of the common offset |
| `generate_pointwise_fixture.py` | Compact two/four/six-spinon refinement and selected independent JMP controls |

Regenerate the committed `pointwise-fixture.json` into a new scratch file:

```sh
python3 experimental/reference/qtm/generate_pointwise_fixture.py \
  --digits 32 --refined-digits 42 --output /path/to/scratch/qtm-fixture.json
```

This checks every PRL partition at both precisions, each full boundary sum
against the vertex target, and one independent JMP partition per sector.
It takes about three minutes on the recorded CPU; it does not run every
six-spinon JMP partition.

These are reference evaluators, not converged outer integrals or positive
Lehmann spectra. Coincident rapidities and singular denominator limits are
explicitly unsupported. Neither pointwise agreement nor a passing fixture
proves a global contour deformation or sector-completeness theorem.
