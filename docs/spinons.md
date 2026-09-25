# Single-spinon dispersion

[Back to the overview](../README.md)

There are two related but different calculations here: an odd periodic XXX
chain supplies a finite-size one-spinon branch, while analytic functions give
the infinite-chain dispersion. Their energy references are different.

```sh
build/bethe-xxx-pbc 65 --spinons
build/bethe-xxx-pbc 15 --spinons --roots --precision long-double
```

The finite branch is not the complete Sz=1/2 spectrum. Begin with
[periodic XXX conventions](xxx.md) if the lattice momentum or quantum-number
labels are unfamiliar.

## Finite odd-chain branch

For odd N, set $M =(N -1)/2$. The one-spinon family occupies M of the M+1 slots
`-M/2, -M/2+1, ..., M/2`, leaving a hole $I_{h}$. Use
`one_spinon_state<Real>(N, hole, options)` for one member, or
`one_spinon_branch<Real>(N, options)` for all $(N +1)/2$ members in ascending k.
Each `SpinonState` contains its `RealState` in `.state`, the exact `.hole`, and:

```math
\begin{aligned}
k&=\frac{\pi}{2}-\frac{2\pi I_h}{N}
&&\text{(spinon momentum)},\\
E_{\mathrm{bulk\text{-}subtracted}}&=E_N-Ne_\infty,
&e_\infty&=\frac14-\log 2.
\end{aligned}
```

The allowed finite-size k values run from $\pi /(2N)$ to $\pi -\pi /(2N)$ in steps
of $2\,\pi /N$. This is a hole-based spinon convention, with
`P = pi*M + pi/2 - k (mod 2*pi)`, **not** momentum relative to the odd-chain
ground state. Both lattice P and spinon k are reported. The bulk-subtracted
energy retains finite-size corrections; it is not obtained by subtracting the
odd-chain ground-state energy, which itself contains a spinon.

## Thermodynamic dispersion

For the analytic zero-field thermodynamic dispersion, include
`<bethe/spinon.hpp>` (also included by `<bethe/heisenberg.hpp>`):

```cpp
long double k = 1.0L;
auto xxx = bethe::heisenberg::spinon_energy(k);          // J=1
auto xxz = bethe::xxz::spinon_energy(k, 0.5L);           // Delta=0.5, J=1
auto scaled = bethe::xxz::spinon_energy(k, 0.5L, 2.0L); // J=2
auto e_inf = bethe::heisenberg::bulk_energy_density<long double>();
```

These implement $\epsilon (k)=J \,\pi \,\sin (k)/2$ for XXX and
$\epsilon (k)=J \,\pi \,\sin (\gamma)\,\sin (k)/(2\,\gamma)$, $\gamma =\arccos (\Delta)$, for gapless
XXZ. They require finite $k$ in `[0,pi]`, positive finite J, and for XXZ
$-1<\Delta \le 1$. The isotropic endpoint uses the explicit XXX limit; k=0 and
k=pi return exactly zero. These functions are analytic, distinct from the
[finite-size XXZ solver](xxz.md).
See [References and provenance](../CITATIONS.md) for the derivations.

Related: [XXX excitation families](excitations.md) and
[reading the numerical report](command-line.md#read-and-save-the-report).
