# Tutorials and model coverage

These worked examples connect actual solver exports to physical questions.
They include plots, downloadable data, reproduction scripts, normalization
checks and caveats for interpreting the result. Start with
[XXZ spinons](xxz-spinons.md) for the export-and-plot workflow, then choose a
model below.

## Available tutorials

- [XXZ spinons and continua](xxz-spinons.md): gapless/gapped bands and two-site folding.
- [Hubbard at half filling](hubbard-half-filled.md): spin–charge separation and energy conventions.
- [Doped Hubbard](hubbard-doped.md): gapless charge branches and chemical-potential subtraction.
- [SU(3)/ULS](su3-uls.md): representations, two-/four-soliton thresholds and three-site folding.
- [Takhtajan–Babujian](takhtajan-babujian.md): spin selection rules and quadrupolar continua.
- [Lieb–Liniger](lieb-liniger.md): particle/hole branches and density scaling.
- [Sutherland](sutherland.md): collision exponents, excitation labels and bounded scans.
- [XXX](xxx.md): finite-ring spinons and periodic/free-end excitation scans.
- [Haldane–Shastry](haldane-shastry.md): motifs, spin multiplicities and complete finite spectra.
- [XYZ](xyz.md): gapped spinons, bound branches and discrete parity copies.
- [q-bosons](q-boson.md): occupation-dependent hopping from free bosons to the phase limit.

## Coverage of implemented models

The aim is a pedagogical tutorial for **every implemented model family**, not
only the thermodynamic excitation tools. This table tracks the remaining
work. A reference guide is not a substitute for a worked tutorial; rows marked
*planned* currently link to the guide only. Proposed but unimplemented models
in the [model survey](../models.md) are outside this inventory.

| Model family | Tutorial or reference guide | Worked topic still to add where planned |
| --- | --- | --- |
| XXX | [Tutorial](xxx.md) | — |
| XXZ | [Tutorial](xxz-spinons.md) | — |
| XYZ | [Tutorial](xyz.md) | — |
| SU(3)/ULS | [Tutorial](su3-uls.md) | — |
| Spin-1 TB | [Tutorial](takhtajan-babujian.md) | — |
| Hubbard | [Half filling](hubbard-half-filled.md), [doping](hubbard-doped.md) | — |
| Lieb–Liniger | [Tutorial](lieb-liniger.md) | — |
| q-boson | [Tutorial](q-boson.md) | — |
| Three-state Potts | [Guide](../potts.md); planned | Charged branch and finite-size scaling |
| Non-Hermitian quantum-group XXZ | [Guide](../xxz-nonhermitian.md); planned | Boundary fields and non-diagonalizable spectra |
| Kondo | [Guide](../kondo.md); planned | Impurity magnetization and susceptibility |
| Sine-Gordon | [Guide](../sine-gordon-excitations.md); planned | Solitons, breathers and finite-volume interpretation |
| Scaling Lee–Yang | [Guide](../lee-yang.md); planned | Vacuum scaling and a non-unitary excitation gap |
| TASEP | [Guide](../tasep.md); planned | Relaxation gaps versus ring size |
| ASEP | [Guide](../asep.md); planned | Drift and relaxation as hopping becomes symmetric |
| Gaudin–Yang | [Guide](../gaudin-yang.md); planned | Interaction dependence of continuum fermions |
| Supersymmetric t–J | [Guide](../tj.md); planned | Doping with double occupancy excluded |
| Negative biquadratic | [Guide](../biquadratic.md); planned | TL sectors, spin content and low excitations |
| Ferromagnetic biquadratic sign | [Guide](../biquadratic-ferromagnetic.md); planned | Defects and bound droplets above the vacuum |
| Richardson pairing | [Guide](../richardson.md); planned | Pair binding and blocked levels |
| Central spin | [Guide](../central-spin.md); planned | Field response in a fixed spin sector |
| SU(n) fermion gas | [Guide](../su-fermions.md); planned | Populations, repulsion and the free limit |
| Bose–Fermi mixture | [Guide](../bose-fermi.md); planned | Composition at equal masses/couplings |
| Integrable ladder | [Guide](../ladder.md); planned | Rung exchange and magnetization sectors |
| Haldane–Shastry | [Tutorial](haldane-shastry.md) | — |
| Sutherland | [Tutorial](sutherland.md) | — |

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
