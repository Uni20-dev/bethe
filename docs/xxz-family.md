# XXZ chains: finite states and thermodynamic excitations

[All models](models.md) · [Spinon tutorial](tutorials/xxz-spinons.md)

The XXZ anisotropy changes both the physics and the root families used by the
solver. Choose the boundary condition and calculation before selecting a tool.

| What to calculate | Program | Coverage and guide |
| --- | --- | --- |
| Periodic ground states and sector minima | `bethe-xxz-pbc` | [Periodic XXZ](xxz.md): Delta>=0, and -1<Delta<0 on even rings |
| Free-end ground states and sector minima | `bethe-xxz-obc` | [Open XXZ](xxz-open.md): Delta>-1, including massive boundary roots |
| Selected real-root excitation scans | Either finite-chain tool | The respective guide restricts these scans to 0<=Delta<=1 |
| Zero-field bulk energy, spinons and two-spinon continuum edges | `bethe-xxz-dispersion` | [Thermodynamic dispersion](xxz-dispersion.md): gapless and massive regimes, Delta>-1 |
| Quantum-group open chain with imaginary end fields | `bethe-xxz-qg-obc` | [Non-Hermitian XXZ](xxz-nonhermitian.md): regular roots and the distinct free-point/Jordan-block calculation |

## MPS and iMPS comparisons

The [spinon tutorial](tutorials/xxz-spinons.md) explains the spin normalization,
physical momentum, two-site-cell folding and continuum boundaries. An envelope
locates allowed energies; it is not an operator-dependent spectral intensity.
The [non-Hermitian tutorial](tutorials/xxz-nonhermitian.md) instead concerns
different boundary fields and spectral/Jordan structure.

## Numerical methods and specialized APIs

These pages explain additional methods or restricted constructions; their
existence does not extend the finite-chain CLI domains in the table above.

- [Massive open-chain roots](xxz-open-massive.md) and
  [negative anisotropy](xxz-negative.md) explain branch choices and limits.
- [Negative odd-ring continuation](xxz-odd-continuation.md) records the restricted
  continuation problem and unresolved coverage outside the public frontend.
- [Open-chain Q-system](xxz-open-qsystem.md) and
  [targeted two-strings](xxz-open-two-string.md) underpin complex-root
  biquadratic/TL calculations.
- [Spin-helix states](xxz-spin-helix.md), [phantom roots](xxz-phantom.md),
  [phantom consistency checks](xxz-phantom-check.md) and
  [phantom wavefunctions](xxz-phantom-wave.md) concern special eigenstates.
- [Coordinate wavefunctions](xxz-coordinate-wave.md),
  [polynomial equations](xxz-polynomial.md),
  [Wronskian constraints](xxz-wronskian.md) and
  [regularity](xxz-regularity.md) discuss how candidate roots are tested.

For the physical spin-1 chains represented by some of these auxiliary equations,
continue with the [biquadratic family guide](biquadratic-family.md).
