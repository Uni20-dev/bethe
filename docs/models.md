# Supported models and calculations

[First calculation](getting-started.md) · [Tutorials](tutorials/index.md) ·
[Arithmetic and execution capabilities](app-capabilities.md)

Choose a physical calculation here, then follow its guide for the Hamiltonian,
normalization, allowed sectors and numerical controls. An implemented model
usually covers selected state families, not its complete spectrum.
The [model roadmap](model-roadmap.md) discusses further integrable models and
extensions; those proposals are separate from the supported calculations below.

## Current coverage

Program names below omit the `bethe-` prefix: `xxx-pbc` means
`bethe-xxx-pbc`. `pbc` means periodic boundaries; `obc` means the open-boundary
problem specified by that guide. Dispersion, continuum and finite-size tools
answer different questions. All numerical frontends support explicit MPFR in
an enabled build; inspect `bethe-info` for the formats actually configured.

## Spin chains

| Family | Programs | Supported calculation and reference | Plotted tutorial |
| --- | --- | --- | --- |
| XXX | `xxx-pbc`, `xxx-obc`, `xxx-structure-factor`, `xxx-structure-factor-thermo` | [Family guide](xxx-family.md): ground/sector minima, selected real-root excitations, finite-ring and thermodynamic two-spinon intensities | [Finite chains](tutorials/xxx.md), [spectral weights](tutorials/xxx-structure-factor.md), [thermodynamic intensity](tutorials/xxx-structure-factor-thermo.md) |
| XXZ | `xxz-pbc`, `xxz-obc`, `xxz-dispersion` | [Family guide](xxz-family.md): anisotropic finite-chain states and zero-field gapless/massive spinons and continuum edges | [Spinons and continua](tutorials/xxz-spinons.md) |
| Quantum-group XXZ | `xxz-qg-obc` | [Non-Hermitian open chain](xxz-nonhermitian.md): regular real-root families; complete free-point fixed-Sz spectra and Jordan-block sizes | [Non-Hermitian spectra](tutorials/xxz-nonhermitian.md) |
| XYZ | `xyz-pbc`, `xyz-dispersion` | [Even-ring ground branch](xyz.md); [thermodynamic spinons, bound branches and two-spinon envelopes](xyz-dispersion.md) | [Bound branches](tutorials/xyz.md) |
| SU(3)/ULS | `su3-pbc`, `su3-dispersion` | [Balanced finite singlet](su3.md); [elementary 3/bar3 lines and two-/four-soliton bounds](su3-dispersion.md) | [ULS excitations](tutorials/su3-uls.md) |
| Spin-1 TB | `tb-pbc`, `tb-dispersion` | [Even-ring singlet with string deviations](takhtajan-babujian.md); [spinons and two-/four-spinon bounds](tb-dispersion.md) | [TB spinons](tutorials/takhtajan-babujian.md) |
| Spin-1 biquadratic | `biquadratic-obc` | [Family guide](biquadratic-family.md): both signs, TL module minima, selected real/complex-root levels, ferromagnetic bands and droplets | [Negative sign](tutorials/biquadratic.md), [ferromagnetic sign](tutorials/biquadratic-ferromagnetic.md) |
| Integrable ladder | `ladder-pbc` | [Sector minima and magnetization](ladder.md), with the required integrable four-spin coupling | [Field response](tutorials/ladder.md) |
| Haldane–Shastry | `haldane-shastry-pbc` | [Exact ring energies and bounded complete motif spectra](haldane-shastry.md), with multiplicities | [Motifs and spins](tutorials/haldane-shastry.md) |
| Three-state Potts | `potts-pbc` | [Periodic vacuum and selected charged one-hole levels](potts.md), resolved by momentum/Z₃ charge | [Charged levels and CFT](tutorials/potts.md) |

## Electrons and continuum gases

| Family | Programs | Supported calculation and reference | Plotted tutorial |
| --- | --- | --- | --- |
| Hubbard | `hubbard-pbc`, `hubbard-obc`, `hubbard-dispersion`, `hubbard-continuum` | [Family guide](hubbard-family.md): finite ground states, half-filled/doped spin-charge lines and half-filled two-particle edges | [Half filling](tutorials/hubbard-half-filled.md), [doping](tutorials/hubbard-doped.md) |
| Lieb–Liniger | `lieb-liniger-pbc`, `lieb-liniger-obc`, `lieb-liniger-dispersion`, `lieb-liniger-thermal` | Repulsive [rings](lieb-liniger.md), [hard walls](lieb-liniger-open.md), [type-I/II thermodynamic curves](lieb-liniger-thermo.md) and [finite-temperature equilibrium](lieb-liniger-thermal.md) | [Particle/hole branches](tutorials/lieb-liniger.md) |
| q-bosons | `q-boson-pbc` | [Fixed-particle-number ground states and canonical excitation scans](q-boson.md), including free and phase-model limits | [Deformation and limits](tutorials/q-boson.md) |
| Gaudin–Yang | `gaudin-yang-pbc` | [Repulsive two-component gas](gaudin-yang.md), interacting odd-population shells and free/polarized limits | [Repulsion and roots](tutorials/gaudin-yang.md) |
| Supersymmetric t–J | `tj-pbc` | [Sector ground states at J=2t](tj.md), supported doped shells and no-hole/polarized sectors | [Fillings and limits](tutorials/tj.md) |
| SU(n) fermion gas | `sun-fermions-pbc` | [Repulsive multicomponent ground states](su-fermions.md), odd occupied populations and free/single-component limits | [Components and seas](tutorials/su-fermions.md) |
| Bose–Fermi mixture | `bose-fermi-pbc` | [Equal-mass/equal-coupling repulsive ground states](bose-fermi.md), supported mixed shells and pure/free limits | [Composition benchmarks](tutorials/bose-fermi.md) |
| Sutherland | `sutherland-pbc` | [Exact bosonic ground/excited energies](sutherland.md), specified labels and bounded label windows | [Collision branches](tutorials/sutherland.md) |

## Pairing and impurities

| Family | Programs | Supported calculation and reference | Plotted tutorial |
| --- | --- | --- | --- |
| Richardson pairing | `richardson` | [Attractive reduced-BCS ground states](richardson.md), distinct levels and specified paired/blocked sector | [Pairing and blocking](tutorials/richardson.md) |
| Central spin | `central-spin` | [Spin-1/2 sector minima](central-spin.md), distinct nonzero bath couplings and a central field | [Field response](tutorials/central-spin.md) |
| Kondo | `kondo-response` | [Zero-temperature scaling-limit impurity energy change, magnetization and susceptibility](kondo.md) | [Impurity response](tutorials/kondo.md) |

## Field theory

| Family | Programs | Supported calculation and reference | Plotted tutorial |
| --- | --- | --- | --- |
| Sine-Gordon | `sine-gordon-vacuum`, `sine-gordon-excited`, `sine-gordon-dispersion`, `sine-gordon-bethe-yang` | [Finite-volume vacuum](sine-gordon.md), [selected exact two-hole levels](sine-gordon-excited.md), and [particle lines, thresholds and asymptotic Bethe–Yang pairs](sine-gordon-excitations.md) | [Particles and finite volume](tutorials/sine-gordon.md) |
| Scaling Lee–Yang | `lee-yang-vacuum`, `lee-yang-excited` | [Periodic vacuum and regular zero-momentum one-particle branch](lee-yang.md) | [Vacuum and gap](tutorials/lee-yang.md) |

## Stochastic systems

| Family | Programs | Supported calculation and reference | Plotted tutorial |
| --- | --- | --- | --- |
| TASEP | `tasep-pbc` | [Leading relaxation branch](tasep.md) on a periodic ring at arbitrary filling | [Relaxation and scaling](tutorials/tasep.md) |
| ASEP | `asep-pbc` | [Bidirectional relaxation branch](asep.md), unequal rates and the symmetric limit | [Bias and symmetry](tutorials/asep.md) |

These tools return stochastic rates/decay eigenvalues, not quantum energies.

## Further models and specialized methods

The [roadmap and literature survey](model-roadmap.md) retains proposed models,
original implementation slices, reference sources and remaining extensions.
Specialized APIs and numerical investigations are linked from the family pages;
an equation being solvable does not imply a supported frontend or full spectrum.

<a id="candidate-index"></a>
<a id="lieb-liniger-the-simplest-new-interacting-family"></a>
<a id="su-n-permutation-chains-and-the-spin-1-uls-point"></a>
Earlier catalogue links now lead to the [candidate index](model-roadmap.md#candidate-index),
[Lieb–Liniger proposal](model-roadmap.md#lieb-liniger-the-simplest-new-interacting-family),
and [SU(n)/ULS proposal](model-roadmap.md#su-n-permutation-chains-and-the-spin-1-uls-point).
