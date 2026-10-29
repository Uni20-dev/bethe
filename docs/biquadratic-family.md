# Biquadratic chains and Temperley–Lieb sectors

[All models](models.md) · [Negative-sign tutorial](tutorials/biquadratic.md) ·
[Ferromagnetic-sign tutorial](tutorials/biquadratic-ferromagnetic.md)

`bethe-biquadratic-obc` covers selected states of both signs of the spin-1 pure
biquadratic chain with free ends. The auxiliary XXZ/TL equations organize the
energies; physical-spin multiplicities are additional information. A TL module
label is not automatically an SU(2) spin label.

## Choose a calculation

| Calculation | Guide |
| --- | --- |
| Negative-sign ground energy, TL module minima and real-root excitations | [Biquadratic chain](biquadratic.md) |
| Complex-root levels and bounded small-chain spectrum searches | [Open-chain Q-system](xxz-open-qsystem.md) |
| Targeted long-chain low-lying singlet | [Two-string solver](xxz-open-two-string.md) |
| Ferromagnetic ground space and exact one-defect band | [Ferromagnetic sign](biquadratic-ferromagnetic.md) |
| Two-, three- and four-defect bound droplets | [Pairs](biquadratic-bound-pairs.md), [triples](biquadratic-bound-triples.md), [quartets](biquadratic-bound-quartets.md) |
| Real-root scattering windows | [Scattering states](biquadratic-scattering.md) |
| Mixed pair/defect, triple/defect and two-pair families | [Pair plus defect](biquadratic-pair-defect.md), [triple plus defect](biquadratic-triple-defect.md), [two pairs](biquadratic-two-pairs.md) |
| Physical SU(2) spin content and multiplicities | [Spin-content guide](biquadratic-spin-content.md) |

Each guide gives its frontend options, state-family restrictions, tables and
failure diagnostics. These modes do not collectively assert a complete excited
spectrum at arbitrary length. All frontend modes support explicit MPFR when
configured, but exact multiplicity counts and labels remain integer quantities.

## Interpret excitation energies

Changing the sign reverses the ordering of auxiliary energies; it does not
turn every AF state label into the corresponding lowest ferromagnetic state.
Use the physical energy and stated gap reference. Bound-cluster deviations may
be much smaller than the root-coordinate spacing; the dedicated deviation
fields and convergence checks matter, even when rounded roots look like an
exact string. The model pages preserve the derivations and validation of each
restricted family.
