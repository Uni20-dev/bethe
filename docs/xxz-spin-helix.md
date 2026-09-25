# Explicit periodic XXZ spin-helix states

[XXZ guide](xxz.md) | [Odd-ring continuation](xxz-odd-continuation.md) | [Regular-state test](xxz-regularity.md)

`periodic_spin_helix<Real>` in
[xxz_spin_helix.hpp](../include/bethe/xxz_spin_helix.hpp) constructs a special
family of exact eigenstates without solving Bethe equations. It supports
either ring parity, every physical magnetization, and all integer windings
`w=0,...,N-1`. It is **not** a ground-state solver or an arbitrary-coupling
excitation solver. There is no separate command-line mode yet.

## From winding to an eigenstate

For the usual $`J =1`$ spin-1/2 Hamiltonian, choose

```math
k=\frac{2\pi w}{N},\qquad \Delta_w=\cos k,\qquad E_w=\frac{N\Delta_w}{4}.
```

[Popkov, Zhang and Klümper](../CITATIONS.md#popkov-2021), equations (4), (5),
and (12), give the spin-helix construction and its fixed-magnetization
projections. Their Pauli-matrix Hamiltonian subtracts the polarized energy;
the energy above restores our spin-1/2 convention.

In the sector with $`M =N /2-S^z`$ down spins, the unnormalized amplitude is

```math
\psi(j_1,\ldots,j_M)=e^{ik(j_1+\cdots+j_M)},\qquad 0\le j_1\lt \cdots\lt j_M\lt N.
```

The state is nonzero: every amplitude has modulus one, and its squared norm
is `binomial(N,M)`. A normalized amplitude would divide by the square root
of this number. The implementation leaves that factor out to avoid
underflow and does not allocate an exponentially large wavefunction.

The local Hamiltonian action cancels around the closed ring. One can also
check this directly on a bit configuration: the numbers of clockwise and
counterclockwise domain walls agree, and each pair contributes a factor
$`\cos (k)-\Delta_{w} =0`$ to `(H-E_w)psi`. Thus every magnetization projection has
the same energy. This establishes an eigenstate, not its position in the
spectrum.

```cpp
#include <bethe/xxz_spin_helix.hpp>
#include <array>

auto const helix = bethe::xxz::periodic_spin_helix<long double>(
    7, 4, uni20::from_twice(std::int64_t{1}));  // Sz=1/2, M=3
std::array<std::size_t, 3> const occupied{0, 2, 5};
auto const amplitude = helix.unnormalized_amplitude(occupied);
// helix.delta = -cos(pi/7), helix.energy = 7*helix.delta/4
// helix.momentum_index = (4*3) mod 7 = 5
```

The momentum convention is the same as the other periodic APIs:
`P=2*pi*momentum_index/N`. Translating the *arguments* of the wavefunction
forward by one site multiplies it by $`\exp (i \,P)`$; active translation of the
ket has the inverse phase. The reflected winding $`N -w`$ reverses the
chirality. Negative Sz is stored as requested, not folded to its partner.

Integer phase sums and products are reduced modulo N before conversion to
`Real`. This avoids both integer-product overflow and large trigonometric
arguments. The arithmetic is native fp64, long double, or fp128. The
ferromagnetic $`w =0`$ and even-ring $`w =N /2`$ endpoints are supported by this
explicit construction even though other XXZ APIs have narrower domains.

## Recognizing the exceptional polynomial

For $`-1\lt \Delta_{w} \le 0`$, the helix has a repeated scaled root

```math
z_w=\cot(k/2)=\frac{\sin k}{1-\cos k}.
```

These are infinite conventional rapidities. In affine coordinates
$`z =o +s \,x`$, its monic polynomial is $`(x -(z_{w} -o)/s)^M`$. The internal
`check_helix_polynomial` compares **all** coefficients with that polynomial
and independently checks the coupling. Finding an endpoint root, matching
only the energy or momentum, or recognizing a root-of-unity anisotropy is
not enough.

Its `compatible` status means numerical agreement at the requested
tolerance, default $`128\,\epsilon`$. The coefficient error is the maximum of
`|c_j-expected_j|/max(1,|c_j|,|expected_j|)`. This is coordinate-dependent
and is not an interval certificate or a bound on extracted roots. Other
statuses distinguish a different coupling, a different polynomial, and
nonfinite arithmetic. Invalid inputs throw.

No input Delta is snapped. The result retains the helix's own coupling and
energy, the absolute coupling mismatch, and

```math
\mathrm{eigenvector\_defect\_bound}=\frac N4|\Delta_{\mathrm{input}}-\Delta_w|.
```

Apart from floating-point roundoff, this bounds
`||(H(Delta_input)-E_w)psi_w||/||psi_w||`: the difference of the two
Hamiltonians is `(Delta_input-Delta_w)*sum Sz_j*Sz_(j+1)`, whose operator
norm is at most $`N \,\lvert \Delta_{\mathrm{input}} -\Delta_{w} \rvert/4`$. This bound refers to the
explicit helix and its energy, not an arbitrary vector reconstructed from
nearby polynomial coefficients.

The odd-ring sector scan tries winding $`(N +1)/2`$ when its regular-state
test is unresolved. At $`\Delta_{w} =-\cos (\pi /N)`$, this recognizes the followed
all-phantom branch. Each entry retains the separate `regularity` and
optional `helix` results; `state_checks_complete` summarizes whether every
entry passes one of these checks or has a resolved
[mixed-phantom numerical witness](xxz-phantom-check.md). It neither changes the
continued coefficients/energies nor supplies a missing continuation result.
`regular_states_complete` retains its stricter, regular-case meaning.

## Tests and remaining coverage

Native-precision tests apply the spin Hamiltonian directly to every winding
and magnetization through seven sites, including both parity classes and
the isotropic endpoints. They also check momentum under translation,
large-integer phase arithmetic, affine polynomial matching, rejection of
perturbed polynomials, opposite chiralities and different couplings, and the
requested tolerance. Continued all-phantom polynomials are recognized in
every folded sector at N=3,5,7,9,13,17,21, in all three scalar types.

Mixed finite/infinite roots, other singular configurations, and reliable
sector-minimum tracking remain separate tasks. A
[mixed-root reduction](xxz-phantom.md) now checks the finite twisted problem
and its separate commensurability condition, without treating that as a
general nonzero-lift test. For example, an odd-ring
continued state at Delta=-1/2 need not be an all-phantom helix. It must not
be accepted just because the coupling is a root of unity. General public
negative-Delta odd-ring ground-state support remains unavailable.
