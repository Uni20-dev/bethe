#!/usr/bin/env python3
"""Check default-OFF test isolation in a built experimental tree, then restore ON.

Run outside CTest: this temporarily reconfigures the supplied build and compares
its discovered tests. No tests are executed and no coverage is removed.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess


def cache_values(build):
    result = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if not line.startswith(("//", "#")) and ":" in line and "=" in line:
            key, value = line.split("=", 1)
            result[key.split(":", 1)[0]] = value
    return result


def discovered_tests(ctest, build):
    payload = subprocess.check_output(
        [ctest, "--test-dir", str(build), "--show-only=json-v1"], text=True)
    ordinary, experimental = Counter(), Counter()
    for test in json.loads(payload)["tests"]:
        if test["name"].endswith("_NOT_BUILT"):
            raise RuntimeError("build every test target before comparing discovered coverage")
        properties = {p["name"]: p["value"] for p in test.get("properties", [])}
        group = experimental if "experimental" in properties.get("LABELS", []) else ordinary
        group[test["name"]] += 1
    return ordinary, experimental


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--ctest", default="ctest")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1]
    build = args.build.resolve()
    cache = cache_values(build)
    if (Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() != source
            or cache.get("BETHE_BUILD_EXPERIMENTAL") != "ON"
            or cache.get("BETHE_BUILD_TESTS") != "ON"):
        raise SystemExit("require this source's built experimental/test configuration")
    before, experimental = discovered_tests(args.ctest, build)
    if not before or not experimental:
        raise SystemExit("build supported and experimental tests before checking isolation")
    configure = [args.cmake, "-S", str(source), "-B", str(build)]
    try:
        # No preset here: presets intentionally opt in; test the option's default.
        subprocess.run(configure + ["-U", "BETHE_BUILD_EXPERIMENTAL"], check=True)
        if cache_values(build).get("BETHE_BUILD_EXPERIMENTAL") != "OFF":
            raise RuntimeError("experimental components are not disabled by default")
        after, remaining = discovered_tests(args.ctest, build)
        if remaining or before != after:
            raise RuntimeError("default-OFF must remove only experimental tests")
        print(f"Default-OFF retains {sum(before.values())} supported tests and removes "
              f"{sum(experimental.values())} experimental tests.", flush=True)
    finally:
        subprocess.run(configure + ["-DBETHE_BUILD_EXPERIMENTAL=ON"], check=True)


if __name__ == "__main__":
    main()
