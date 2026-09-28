# Maintaining tutorials and the documentation site

Tutorials are ordinary GitHub-readable Markdown. Keep the narrative, small
reference exports, and SVG figures in `docs/tutorials/`; put executable
reproduction scripts in `scripts/`. Readers should not need a solver build
or a notebook server just to see a result.

## Numerical figures

Use a frontend's CSV/TSV/JSON exports as the data source. Include its provenance
and convergence information, check unsuccessful exit codes and missing values,
and explain energy, momentum and symmetry conventions in the tutorial.
Do not imply spectral weights when a solver supplies only kinematic bounds.

Use bare `bethe-*` command names in calculation examples, with no assumed build
directory. Keep source-relative Python commands in reproduction sections and
identify their working directory. Downloaded CSV provenance retains the actual
command used to generate it; those historical paths are not setup instructions.

### Plotting environment

From the source checkout, activate a Python environment with the plotting
dependencies. For example:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r docs/tutorials/requirements.txt
```

The environment can live elsewhere; `.venv` is just a conventional local choice.
Solvers themselves do not need Python. To regenerate data, supply a command on
`PATH` or an executable path with `--solver`; scripts that use multiple frontends
accept their containing directory. These can be build-tree or installed programs.

### Validate saved data

Validate the saved numerical examples:

```sh
python3 scripts/plot_xxz_tutorial.py --check
python3 scripts/plot_hubbard_tutorial.py --check
python3 scripts/plot_spin1_tutorial.py --check
python3 scripts/plot_lieb_liniger_tutorial.py --check
python3 scripts/plot_sutherland_tutorial.py --check
python3 scripts/plot_xxx_tutorial.py --check
python3 scripts/plot_haldane_shastry_tutorial.py --check
python3 scripts/plot_xyz_tutorial.py --check
python3 scripts/plot_q_boson_tutorial.py --check
python3 scripts/plot_biquadratic_tutorial.py --check
python3 scripts/plot_field_theory_tutorial.py --check
python3 scripts/plot_exclusion_tutorial.py --check
python3 scripts/plot_fermion_tutorial.py --check
python3 scripts/plot_multicomponent_tutorial.py --check
python3 scripts/plot_potts_qg_tutorial.py --check
python3 scripts/plot_gaudin_tutorial.py --check
python3 scripts/plot_kondo_ladder_tutorial.py --check
python3 -m unittest discover -s scripts -p 'test_tutorial*.py'
```

These checks require neither plotting packages nor a C++ build. To redraw or
regenerate, follow the [XXZ tutorial](xxz-spinons.md#5-reproduce-the-figures) or
the [Hubbard tutorial](hubbard-half-filled.md#5-reproduce-and-check).
The [ULS tutorial](su3-uls.md#5-reproduce-the-plots-and-check-the-comparison)
also describes regenerating the ULS and TB figures together.
The [Lieb–Liniger](lieb-liniger.md#5-reproduce-the-figures) and
[Sutherland](sutherland.md#6-reproduce-and-check) tutorials each have a
standalone export-and-plot script.
The [XXX](xxx.md#5-reproduce-and-validate) and
[Haldane–Shastry](haldane-shastry.md#6-reproduce-the-figures) examples also show
how to join named auxiliary tables. Their regeneration scripts share a
temporary-file capture helper and validate all exports before saving them.
The [XYZ](xyz.md#5-reproduce-and-check) tutorial checks branch/parity copies;
the [q-boson](q-boson.md#5-reproduce-the-figures) tutorial joins roots and state
tables and checks spectral moments against occupation-basis hopping.
Both [biquadratic sign tutorials](biquadratic.md#6-reproduce-both-sign-tutorials)
share one export/plot script, checking physical-spin counts, sign conventions,
Q-system completeness diagnostics and targeted finite-string branches.
The [Sine-Gordon](sine-gordon.md#5-reproduce-and-validate) and
[Lee–Yang](lee-yang.md#5-reproduce-and-validate) tutorials share field-theory
table validation, including independently converged levels, source diagnostics,
vacuum subtraction and mass/circumference scaling.
The [TASEP](tasep.md#5-reproduce-and-check) and [ASEP](asep.md#5-reproduce-and-validate)
tutorials share relaxation/root validation, including analytic endpoints,
stationary-only sectors, original Bethe equations and a small Markov-matrix check.
The [Gaudin–Yang](gaudin-yang.md#5-reproduce-and-validate) and
[t–J](tj.md#5-reproduce-and-validate) examples share nested-table checks while
keeping their equations, free-shell branches, units and momentum conventions distinct.
The [SU(n) fermion](su-fermions.md#5-reproduce-and-validate) and
[Bose–Fermi mixture](bose-fermi.md#5-reproduce-and-validate) examples share
multicomponent export validation, checking population-to-sea mappings, graded
equations and independent Gaudin–Yang/Lieb–Liniger reductions.
The [Potts](potts.md#5-reproduce-and-validate) and
[non-Hermitian XXZ](xxz-nonhermitian.md#5-reproduce-and-validate) examples
distinguish selected regular families from complete endpoint sectors. Their
shared script checks charge/momentum selection, root joins and Jordan-block
dimensions; `--oracle` adds independent small clock/spin matrices using NumPy.
The [Richardson](richardson.md#5-reproduce-and-validate) and
[central-spin](central-spin.md#5-reproduce-and-validate) examples share
Gaudin-variable validation but retain their distinct sectors and Hamiltonians.
Their `--oracle` checks every saved example against independent occupation/spin
matrices, including a central-polarization check of the energy derivative.
The [Kondo](kondo.md#5-reproduce-and-validate) and
[ladder](ladder.md#5-reproduce-and-validate) tutorials check subtracted-response
units and fixed-sector field envelopes respectively. The ladder `--oracle`
minimizes independent color matrices over all allowed populations; the
separate `reference_kondo.py` uses mpmath for the continuum response.
`scripts/tutorial_common.py` shares CSV/provenance validation and deterministic
SVG writing; each model's script and tests own its schema and physics checks.

## Build the site locally

From the source checkout, use a Python environment and an output directory of
your choice. Here both `.venv` and `.cache` are ignored by Git:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r docs/requirements-site.txt
python3 -m unittest discover -s scripts -p 'test_build_docs.py'
python3 scripts/build_docs.py --work-dir .cache/docs-site
python3 scripts/check_docs_site.py .cache/docs-site/site
python3 -m http.server 8000 --directory .cache/docs-site/site
```

Open `http://localhost:8000/`. The build stages **Git-tracked** documentation
only; `git add` new pages and assets before previewing them. Existing tracked
files use their current working-tree content. Local untracked drafts cannot
accidentally enter the site. Staging uses a fresh temporary directory per build.

`zensical.toml` controls the navigation. The build adapter keeps the repository's
protected inline math and fenced `math` blocks intact in source, converting
them only in staging to MathJax wrappers. Code examples remain code. MathJax
3.2.2 is loaded from a CDN, so equations need network access in a browser;
figures and page content are served locally. Run `scripts/check_markdown_math.py`
to check GitHub-safe source syntax (optionally `--github` for GitHub's renderer).

## Publish

The Pages workflow validates the committed tutorial data and math adapter,
then builds the site on pull requests. Pushes to `main` affecting documentation
also deploy it; manual dispatch is available. Deployment uses the GitHub Pages
environment and a Pages artifact, not a generated branch.

The site build deliberately **does not compile solvers or regenerate figures**.
Numerical regeneration is an explicit, reviewable update to data and figures.
This keeps prose-only publishing fast and independent of Uni20/MPLAPACK builds.
