"""Shared momentum inputs through both production frontends and export validation."""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

spectrum, moments, mpfr = sys.argv[1:]
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", UNI20_COLOR="never")


def invoke(program, args, success=True):
    result = subprocess.run([program, *map(str, args)], capture_output=True,
                            text=True, env=env, timeout=60)
    assert (result.returncode == 0) == success, (args, result.returncode, result.stderr)
    return result


def document(program, args):
    return json.loads(invoke(program, [*args, "--format", "json"]).stdout)


def rows(doc, table):
    data = doc["tables"][table]
    keys = [c["id"] for c in data["columns"]]
    return [dict(zip(keys, row)) for row in data["rows"]]


for program in (spectrum, moments):
    help_text = invoke(program, ["--help"]).stdout
    for option in ("--momentum-from-pi", "--momentum-beta", "--momentum-beta-tolerance"):
        assert option in help_text

for sectors in ("2", "2,4,6"):
    doc = document(spectrum, ["--momentum-beta", "0", "--omega", "0", "--sectors", sectors])
    metadata = doc["tables"]["spectrum"]["metadata"]
    assert metadata["Momentum input kind"] == "beta"
    assert metadata["Momentum input"] == "0"
    assert metadata["Folded momentum beta"] == "0"
    assert abs(float(rows(doc, "spectrum")[0]["q"]) - math.pi / 2) < 1e-15
    for text in ("1.8", "2", "3.2"):
        ordinary = document(spectrum, ["--momentum", text, "--omega", "0", "--sectors", sectors])
        assert float(rows(ordinary, "spectrum")[0]["q"]) == float(text)

# Conversion failure must precede opening any file, including --force exports.
with tempfile.TemporaryDirectory() as tmp:
    protected = Path(tmp) / "prior.json"
    absent = Path(tmp) / "new.csv"
    for program, extra in ((spectrum, ["--omega", "0"]), (moments, [])):
        for sector in ("2", "6"):
            for query in (["--momentum-from-pi", "1e-100"], ["--momentum", "3.141592653589793"]):
                protected.write_text("preserve prior result\n")
                result = invoke(program, [*extra, *query, "--sectors", sector, "--precision", "fp64",
                                          "--json", protected, "--csv", absent, "--force"], False)
                assert "scalar-q backend" in result.stderr
                assert protected.read_text() == "preserve prior result\n"
                assert not absent.exists()

invoke(spectrum, ["--momentum-beta", "0", "--momentum-points", "3"], False)
invoke(moments, ["--momentum-from-pi", "1", "--global"], False)
for value in ("-5", "1e-8"):
    invoke(spectrum, ["--momentum-points", "3", "--points", "2", "--momentum-beta-tolerance", value], False)
    invoke(moments, ["--global", "--momentum-beta-tolerance", value], False)

# Scalar-q loss is observable separately from accurate input-to-beta conversion.
warning = invoke(spectrum, ["--momentum-from-pi", "1e-10", "--omega", "0", "--format", "json"])
assert "scalar momentum geometry" in warning.stderr
metadata = json.loads(warning.stdout)["tables"]["spectrum"]["metadata"]
assert metadata["Momentum conversion warning"] == "true"
assert float(metadata["Estimated beta conversion error"]) < 1e-12

# The native reference's own approximation of pi must not hide scalar-q loss
# merely because the user supplied beta directly, without any input subtraction.
warning = invoke(spectrum, ["--momentum-beta", "35.80001810975594200500693453648182185134",
                            "--omega", "0", "--format", "json"])
assert "scalar momentum geometry" in warning.stderr
metadata = json.loads(warning.stdout)["tables"]["spectrum"]["metadata"]
assert metadata["Momentum conversion warning"] == "true"
assert float(metadata["Estimated scalar-q beta error"]) > 8e-5

# Actual nonzero moments must be unchanged by the choice of input coordinate.
controls = ["--sectors", "2", "--two-order", "4", "--two-gap-cutoff", "1"]
a = rows(document(moments, [*controls, "--momentum-beta", "0"]), "moments")
b = rows(document(moments, [*controls, "--momentum", str(math.pi / 2)]), "moments")
for x, y in zip(a, b, strict=True):
    assert float(x["value"]) > 0
    assert math.isclose(float(x["value"]), float(y["value"]), rel_tol=1e-12)

if mpfr.upper() in ("ON", "TRUE", "1"):
    for sectors in ("2", "2,4,6"):
        doc = document(spectrum, ["--momentum-from-pi", "1e-100", "--omega", "0",
                                  "--precision", "bits-512", "--sectors", sectors])
        metadata = doc["tables"]["spectrum"]["metadata"]
        assert metadata["Momentum input"] == "1e-100"
        assert metadata["MPFR input conversion"] == "true"
        assert metadata["Momentum conversion warning"] == "false"
        assert float(metadata["Folded momentum beta"]) > 230
