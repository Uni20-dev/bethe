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

Validate the saved numerical examples:

```sh
python3 scripts/plot_xxz_tutorial.py --check
python3 scripts/plot_hubbard_tutorial.py --check
python3 -m unittest discover -s scripts -p 'test_tutorial*.py'
```

These checks require neither plotting packages nor a C++ build. To redraw or
regenerate, follow the [XXZ tutorial](xxz-spinons.md#5-reproduce-the-figures) or
the [Hubbard tutorial](hubbard-half-filled.md#5-reproduce-and-check).
`scripts/tutorial_common.py` shares CSV/provenance validation and deterministic
SVG writing; each model's script and tests own its schema and physics checks.

## Build the site locally

Use local storage for `build_codex`, following the project's build convention:

```sh
python3 -m venv build_codex/docs-venv
build_codex/docs-venv/bin/pip install -r docs/requirements-site.txt
python3 -m unittest discover -s scripts -p 'test_build_docs.py'
build_codex/docs-venv/bin/python scripts/build_docs.py
python3 scripts/check_docs_site.py
python3 -m http.server 8000 --directory build_codex/docs-site/site
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
