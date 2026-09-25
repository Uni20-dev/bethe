# Mixed finite and infinite XXZ roots

[Odd-ring continuation](xxz-odd-continuation.md) | [Spin helices](xxz-spin-helix.md) | [Regular-state test](xxz-regularity.md)

The internal `reduce_phantom_polynomial` in
[xxz_phantom.hpp](../include/bethe/xxz_phantom.hpp) removes a specified
one-sided cluster of infinite rapidities and checks the remaining, twisted
finite-root problem. This is a necessary-equation diagnostic, **not a
nonzero-state lifting test or a ground-state solver**. Reduction alone does
not make an unresolved entry in the sector scan acceptable; a separate
[numerical witness](xxz-phantom-check.md) checks the dressed amplitude.

## The extra condition left by infinite roots

Write $\Delta =\cos (\gamma)$ with $0<\gamma <\pi$. Suppose M roots consist of
p phantom roots at the same infinity and r=M-p finite roots.
[Popkov, Zhang and Klümper](../CITATIONS.md#popkov-2021), equation (15) and
the following reduced equations, give a commensurability condition and a
twisted finite-root problem. In our momentum and endpoint conventions:

```math
\begin{aligned}
\chi&=\pm1,\qquad z_{\mathrm{phantom}}=\chi\sqrt{\frac{1+\Delta}{1-\Delta}},\\
e^{ik_{\mathrm{phantom}}}&=e^{i\chi\gamma},\qquad e^{i(N-2r)\gamma}=1,\\
e^{iNk_j}&=e^{-2i\chi p\gamma}\prod_{l\ne j}S(k_j,k_l).
\end{aligned}
```

The integer p counts roots at **one** endpoint, not an arbitrary mixture
of the two infinities. The commensurability condition depends on the
number of finite roots r; testing only whether $\exp (i \,N \,\gamma)=1$ would
miss mixed states. At odd N and minimal $\lvert S^z \rvert =1/2$, a single phantom root
has $N -2\,r =3$, hence the familiar Delta=-1/2 exception.

For our minimal-|Sz| continuation, successive candidate clusters occur at
$\Delta =-\cos (\pi /(2\,p +1))$. This equation-level observation does not identify
every possible singular branch or establish a sector minimum.

## Deflation and the rotated polynomial equations

For input `F(x)` with $z =o +s \,x$, monic synthetic division removes p copies
of $x -(z_{\mathrm{phantom}} -o)/s$. The code reconstructs the full polynomial and
compares every coefficient, with error

```math
\max_j\frac{|\mathrm{original}_j-\mathrm{reconstructed}_j|}
{\max(1,|\mathrm{original}_j|,|\mathrm{reconstructed}_j|)}.
```

Thus a guessed multiplicity is not accepted just because F is small at an
endpoint. A different chirality or an excessive multiplicity must pass the
same complete-factorization check. The reduced coefficients retain the
original affine coordinate; no individual finite roots are extracted.

The [polynomial engine](xxz-polynomial.md) now has
`evaluate_rotated(coefficients, Delta, rotation, jacobian)`. It takes the
same parity-selected real or imaginary component as before, but of
`rotation*G`, where `G=(1+i*z)^N*K_- mod F`. For a boundary phase phi,

```math
\mathrm{rotation}=e^{i\phi/2},\qquad
\mathrm{rotation}\,G-\sigma\,\overline{\mathrm{rotation}\,G}=0,\qquad
\sigma=(-1)^{N-M-1}.
```

Phantom reduction uses `rotation=exp(-i*chi*p*gamma)`. Integer powers of
$\Delta +i \,\chi \,\sqrt{1-\Delta ^{2}}$ generate the required phases without an
inverse trigonometric function. The separate commensurability residual is
checked before normalizing radial roundoff in the rotation.

The analytic Jacobian includes this fixed rotation; it still differentiates
the coefficient-dependent polynomial reduction. `evaluate` remains exactly
the zero-twist call. Energy and total-momentum formulas are unchanged, but
`momentum_defect` is specifically a **zero-twist** check and must not be
applied as a quantization criterion to the reduced twisted state.

The regular-state diagnostic accepts the same optional rotation. A
diagonal twist multiplies the algebraic vacuum eigenvalues by constant
nonzero phases; the logarithmic derivative and the regular norm argument
are otherwise unchanged. This applies to the reduced state, not to the
singular map that adds phantom roots back.

## Tolerances and returned evidence

`regular_reduced_equations` requires a successful factorization, the
commensurability condition, and a `regular_on_shell` reduced state.
Other statuses distinguish `no_endpoint_factor`, `phase_mismatch`,
`unresolved_reduced_state`, and `nonfinite`. Invalid arguments throw.
The function supports $0<p <M$ and $-1<\Delta \le 0$; the all-phantom construction
already has its own [helix implementation](xxz-spin-helix.md).

The output retains the finite coefficients, rotation, reconstruction and
commensurability errors, finite-root regularity result, and the difference
between full and reduced energies. Infinite roots contribute zero energy
on the exact limiting solution. That energy difference is an additional
diagnostic, not an independently certified energy-error bound.

The factor/phase tolerance defaults to $128\,\epsilon$; the reduced
coefficient-residual tolerance defaults to $32\,\epsilon$. Neither is
loosened automatically, and the requested coupling is never replaced.
Deflation and the changed residual normalization can amplify roundoff:
for example, the N=21, p=9 fp64 audit has a reduced residual about
2.65e-14, above its default 7.11e-15 threshold. It correctly stays
unresolved at that tolerance. A separate regression requests $512\,\epsilon$
for the larger-cluster audit and checks that the coefficients do not change.
This is not a hidden fallback or evidence of exact commensurability.

## An independent mixed-state wavefunction

For r=1, our coordinate-ansatz deduction gives an explicit useful check.
Put `k0=chi*gamma`, $M =p +1$, and choose the finite momentum k satisfying
`exp(i*N*k)=exp(-2*i*p*k0)`. On ordered occupied sites, with rank a starting
at zero, the unnormalized wavefunction is

```math
\psi(j_0,\ldots,j_p)=e^{ik_0\sum_a j_a}
\sum_{a=0}^{p}e^{2ik_0a+i(k-k_0)j_a}.
```

The rank factor is the amplitude ratio for moving the finite magnon past
a phantom. Native-precision tests check this formula by applying the spin
Hamiltonian directly, without reusing the reduced Bethe residual. Both
chiralities and the regular finite-momentum modes through N=9 are covered,
with multiple phantom roots. Endpoint finite momenta and the infinite
polynomial coordinate at k=0 are explicitly excluded from this regular
finite-root check.

Additional tests compare the rotated coefficient equations with direct
rootwise products (including conjugate roots), their Jacobian with finite
differences, and the one-magnon twist with exact momentum quantization.
Continuation/deflation audits cover all mixed multiplicities in the
minimal-|Sz| branch at N=5,7,9,13,17,21; p=1 passes the default tolerance.
Deliberately off-shell finite roots and incorrect commensurability are
rejected independently of successful endpoint deflation.

A [general coordinate-space dressing map](xxz-phantom-wave.md) now
constructs the lifted amplitudes from a finite-state callback. Full
Hamiltonian-identity tests and nonzero lifts of continued two-finite-root
states provide additional independent checks. The map also has a genuine
kernel, exhibited explicitly in the tests. A bounded numerical witness
search now checks amplitudes against uncertainty estimates; it is not a
proof of general nonzero lifting. The
remaining singular-root families, and reliable sector-minimum tracking
remain work to do before public odd-ring ground-state integration.
