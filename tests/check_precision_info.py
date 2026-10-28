"""Build introspection and strict fixed-format/arbitrary precision selection."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

info, native, structure_factor, fp128, arbitrary, mpfr, mplapack_fp128 = sys.argv[1:]
enabled = lambda value: value.upper() in ("ON", "TRUE", "1")
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def run(program, args, code=0):
    result = subprocess.run([program, *args], capture_output=True, text=True, env=env, timeout=60)
    assert result.returncode == code, (args, result.stdout, result.stderr)
    return result


def document(program, args):
    return json.loads(run(program, [*args, "--format", "json"]).stdout)["tables"]


screen = run(info, []).stdout
assert "Bethe build and arithmetic" in screen and "long-double" in screen
assert "bethe-info" in run(info, ["--help"]).stdout
assert "bethe-info" in run(info, ["--version"]).stdout
table = document(info, [])["precisions"]
columns = [c["id"] for c in table["columns"]]
rows = {row[0]: dict(zip(columns, row)) for row in table["rows"]}
for row in rows.values():
    for key in ("significand_bits", "decimal_digits", "min_exponent", "max_exponent", "storage_bytes"):
        if row[key] is not None:
            row[key] = int(row[key])
assert set(rows) == {"fp64", "fp80", "long-double", "fp128", "bits-N", "digits-N"}
assert rows["fp64"]["available"] and rows["fp64"]["significand_bits"] == 53
assert rows["fp128"]["available"] == enabled(fp128)
assert rows["bits-N"]["available"] == rows["digits-N"]["available"] == enabled(mpfr)
metadata = table["metadata"]
assert (metadata["Uni20 float128 configured"] == "true") == enabled(fp128)
assert (metadata["MPLAPACK binary128 configured"] == "true") == enabled(mplapack_fp128)
assert (metadata["Uni20 arbitrary-precision real/complex configured"] == "true") == enabled(mpfr)
assert (metadata["MPLAPACK MPFR configured"] == "true") == enabled(arbitrary)
assert ("bethe-xxx-pbc" in metadata["Arbitrary-precision frontends"]) == enabled(mpfr)
assert ("bethe-xxx-structure-factor" in metadata["Arbitrary-precision frontends"]) == enabled(arbitrary)
for name in ("bethe-richardson", "bethe-central-spin", "bethe-xxx-obc"):
    assert (name in metadata["Arbitrary-precision frontends"]) == enabled(mpfr)
for name in ("bethe-lieb-liniger-pbc", "bethe-lieb-liniger-obc", "bethe-lieb-liniger-dispersion",
             "bethe-su3-pbc", "bethe-xxz-pbc", "bethe-xxz-obc"):
    assert (name in metadata["Arbitrary-precision frontends"]) == enabled(arbitrary)
if enabled(fp128):
    assert rows["fp128"]["significand_bits"] == 113
long_double = rows["long-double"]
is80 = (long_double["significand_bits"], long_double["min_exponent"], long_double["max_exponent"]) == (64, -16381, 16384)
assert rows["fp80"]["available"] == is80
if is80:
    for program in (native, structure_factor):
        a = document(program, ["4", "--precision", "fp80"])
        b = document(program, ["4", "--precision", "long-double"])
        assert set(a) == set(b)
        for name in a:
            assert a[name]["rows"] == b[name]["rows"]
            metadata = a[name]["metadata"]
            assert metadata["Precision"] == "fp80"
            assert metadata["Arithmetic format"] == "extended80"
            assert int(metadata["Working precision (bits)"]) == 64
else:
    assert "fp80 is unavailable" in run(native, ["4", "--precision", "fp80"], 1).stderr

with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/"preserve.json"
    path.write_text("sentinel")
    invalid = ["mpfr", "bits-0", "digits--100"]
    if not enabled(mpfr):
        invalid += ["bits-512", "digits-100"]
    for value in invalid:
        run(native, ["4", "--precision", value, "--json", str(path), "--force"], 1)
        assert path.read_text() == "sentinel"
    if not enabled(fp128):
        assert "fp128 is unavailable" in run(native, ["4", "--precision", "fp128"], 1).stderr
    if not enabled(arbitrary):
        for value in ("bits-512", "digits-100"):
            result = run(structure_factor, ["4", "--precision", value, "--json", str(path), "--force"], 1)
            assert "arbitrary-precision form factors require" in result.stderr
            assert "UNI20_ENABLE_MPREAL=ON" in result.stderr
            assert "UNI20_BACKEND_MPLAPACK=ON" in result.stderr
            assert path.read_text() == "sentinel"

print("Precision selection and build information contracts passed")
