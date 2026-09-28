# Tutorials and model coverage

These worked examples connect actual solver exports to physical questions.
They include plots, downloadable data, reproduction scripts, normalization
checks and caveats for interpreting the result. Start with
[XXZ spinons](xxz-spinons.md) for the export-and-plot workflow, then choose a
model below.

## Running the examples

Commands use `bethe-*` names on `PATH`, not a particular build or install
directory. If you are working directly in your build directory, either prepend
`./` to those names or add that directory to `PATH` as described in
[the build guide](../building.md#run-the-programs). Run calculations wherever you
want their output files; the source checkout is not required.

The optional figure-reproduction commands are different: run `python3 scripts/…`
from the **source checkout**, with the [plotting dependencies](contributing.md#plotting-environment)
installed. A single-tool script accepts `--solver bethe-…` from `PATH` or an
explicit executable path. Multi-tool scripts take an executable directory via
`--solver-dir` (or the documented `--bin-dir`/`--build-dir` option): replace
`/path/to/bethe/bin` with your actual build or installed binary directory; it
need not be named `bin`. Without these options, scripts read the checked-in data.
Regeneration replaces the repository's example data, not files in your calculation directory.

## Available tutorials

- [XXZ spinons and continua](xxz-spinons.md): gapless/gapped bands and two-site folding.
- [Hubbard at half filling](hubbard-half-filled.md): spin–charge separation and energy conventions.
- [Doped Hubbard](hubbard-doped.md): gapless charge branches and chemical-potential subtraction.
- [SU(3)/ULS](su3-uls.md): representations, two-/four-soliton thresholds and three-site folding.
- [Takhtajan–Babujian](takhtajan-babujian.md): spin selection rules and quadrupolar continua.
- [Lieb–Liniger](lieb-liniger.md): particle/hole branches and density scaling.
- [Sutherland](sutherland.md): collision exponents, excitation labels and bounded scans.
- [XXX](xxx.md): finite-ring spinons and periodic/free-end excitation scans.
- [XXX spectral weights](xxx-structure-factor.md): two-spinon intensities, sum rules and plotting broadening.
- [Thermodynamic XXX intensity](xxx-structure-factor-thermo.md): exact two-spinon spectral density without broadening.
- [Haldane–Shastry](haldane-shastry.md): motifs, spin multiplicities and complete finite spectra.
- [XYZ](xyz.md): gapped spinons, bound branches and discrete parity copies.
- [q-bosons](q-boson.md): occupation-dependent hopping from free bosons to the phase limit.
- [Negative biquadratic chain](biquadratic.md): TL sectors, spin multiplicities and complex-root singlets.
- [Ferromagnetic biquadratic chain](biquadratic-ferromagnetic.md): one-defect bands and bound droplets.
- [Sine-Gordon](sine-gordon.md): breather thresholds and exact versus asymptotic finite-volume gaps.
- [Scaling Lee–Yang](lee-yang.md): effective central charge and a nonunitary one-particle gap.
- [TASEP](tasep.md): relaxation times, oscillation and density-dependent size scaling.
- [ASEP](asep.md): a fixed-total-rate bias scan and the symmetric limit.
- [Gaudin–Yang](gaudin-yang.md): repulsion, nested roots and finite-ring shell conventions.
- [Supersymmetric t–J](tj.md): selected dopings, projected fermions and the no-hole XXX limit.
- [SU(n) fermions](su-fermions.md): component populations, nested seas and cross-model reductions.
- [Bose–Fermi mixtures](bose-fermi.md): equal-coupling composition benchmarks and particle statistics.
- [Three-state Potts](potts.md): charged branches, discrete symmetries and finite-size CFT estimators.
- [Non-Hermitian XXZ](xxz-nonhermitian.md): regular-root families, Jordan blocks and parity-dependent Casimir terms.
- [Richardson pairing](richardson.md): blocked sectors and regular variables through pair-root collisions.
- [Central spin](central-spin.md): fixed-sector field response and energy-derived central polarization.
- [Kondo](kondo.md): host subtraction, susceptibility calibration and universal impurity response.
- [Integrable ladder](ladder.md): field envelopes, magnetization steps and highest-weight descendants.

## Coverage of implemented models

All **26 implemented model families** have worked tutorials, with 28 pages:
Hubbard has separate half-filled/doped examples and XXX has energies/weights. This inventory
includes finite systems and impurity response, not only thermodynamic
excitation tools. Proposed but unimplemented models in the
[model survey](../models.md) are outside this inventory.

| Model family | Worked tutorial |
| --- | --- |
| XXX | [Energies](xxx.md), [finite spectral weights](xxx-structure-factor.md), [thermodynamic intensity](xxx-structure-factor-thermo.md) |
| XXZ | [Tutorial](xxz-spinons.md) |
| XYZ | [Tutorial](xyz.md) |
| SU(3)/ULS | [Tutorial](su3-uls.md) |
| Spin-1 TB | [Tutorial](takhtajan-babujian.md) |
| Hubbard | [Half filling](hubbard-half-filled.md), [doping](hubbard-doped.md) |
| Lieb–Liniger | [Tutorial](lieb-liniger.md) |
| q-boson | [Tutorial](q-boson.md) |
| Three-state Potts | [Tutorial](potts.md) |
| Non-Hermitian quantum-group XXZ | [Tutorial](xxz-nonhermitian.md) |
| Kondo | [Tutorial](kondo.md) |
| Sine-Gordon | [Tutorial](sine-gordon.md) |
| Scaling Lee–Yang | [Tutorial](lee-yang.md) |
| TASEP | [Tutorial](tasep.md) |
| ASEP | [Tutorial](asep.md) |
| Gaudin–Yang | [Tutorial](gaudin-yang.md) |
| Supersymmetric t–J | [Tutorial](tj.md) |
| Negative biquadratic | [Tutorial](biquadratic.md) |
| Ferromagnetic biquadratic sign | [Tutorial](biquadratic-ferromagnetic.md) |
| Richardson pairing | [Tutorial](richardson.md) |
| Central spin | [Tutorial](central-spin.md) |
| SU(n) fermion gas | [Tutorial](su-fermions.md) |
| Bose–Fermi mixture | [Tutorial](bose-fermi.md) |
| Integrable ladder | [Tutorial](ladder.md) |
| Haldane–Shastry | [Tutorial](haldane-shastry.md) |
| Sutherland | [Tutorial](sutherland.md) |

A model tutorial introduces a useful supported calculation; it does not imply
that every frontend mode or every physical sector is implemented. Boundary
conditions, finite-size solvers, thermal modes and specialised spectrum searches
remain documented in their linked reference guides and can support additional
tutorials as the collection grows.

## Reproducibility and publication

Each completed tutorial must have runnable commands, stated energy/momentum
conventions, checked numerical exports, interpretation of the result and links
to the literature. Where figures are used, their scripts read the exports
rather than silently substituting an independent model calculation.

Major additions are committed, pushed and published to GitHub Pages, with
checks of the deployed pages. Numerical data and figures are regenerated
explicitly; ordinary site builds do not compile the solvers. See
[maintaining tutorials](contributing.md) for the workflow.
