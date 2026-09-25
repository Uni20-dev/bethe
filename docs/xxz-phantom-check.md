# Numerical witnesses for mixed phantom states

[Phantom reduction](xxz-phantom.md) | [Dressing map](xxz-phantom-wave.md) | [Odd-ring scan](xxz-odd-continuation.md)

`check_phantom_lift<Real>` in
[xxz_phantom_check.hpp](../include/bethe/xxz_phantom_check.hpp) combines
endpoint deflation, finite-root regularity, root recovery, and coordinate
dressing. Its purpose is to find a **numerically resolved nonzero amplitude**,
not to assume that the phantom map has no kernel. It does not allocate a
spin basis or diagonalize a Hamiltonian.

The input includes an explicit phantom count p and chirality. The
[reduction](xxz-phantom.md) must first pass its factorization,
commensurability, and finite twisted-equation checks. The reduced
[root recovery](polynomial-roots.md) must also return `resolved`.
Neither the input coefficients nor Delta are changed.

## Resolution, not merely a nonzero floating-point value

For a recovered root x in affine coordinates `z=center+scale*x`, use its
recovery diagnostic's largest root-uncertainty estimate, transformed to z,
plus a native rounding allowance. For a radius rho around z, the momentum
factor $`v =-(1-i \,z)/(1+i \,z)`$ has the variation estimate

```math
\mathrm{radius}_v=\frac{2\rho}{|1+iz|(|1+iz|-\rho)}.
```

The denominator must be positive. Conversion roundoff is added separately.
This formula is an exact perturbation inequality for a supplied disk in
exact arithmetic; the input disks here are **numerical estimates**, not
certified root enclosures.

The coordinate evaluator now accepts optional momentum radii. It propagates
variation through the plane powers, pair factors, and subset recurrence.
For two uncertain factors, the product radius is
`abs(a)*radius_b + abs(b)*radius_a + radius_a*radius_b`. For the pair factor,

```math
\mathrm{radius}_F(i,j)=|v_i|\mathrm{radius}_j+|v_j|\mathrm{radius}_i
+\mathrm{radius}_i\mathrm{radius}_j+2|\Delta|\mathrm{radius}_i.
```

Divide by the **nominal** pair normalization, holding that common global
scale fixed. Powers use repeated uncertain multiplication, avoiding the
cancellation in `(abs(v)+radius)^x-abs(v)^x`. This propagates nonlinear
variation, not just the first derivative. Exact arithmetic would give an
envelope for all momentum factors in the supplied disks; the implementation
uses ordinary rounding, not outward-rounded interval arithmetic.

The dressing sum accumulates the finite-amplitude variations and all
absolute permutation terms. A separate heuristic arithmetic allowance is

```math
\mathrm{amplitude\_tolerance}
\left[1+NM+r^2+\binom Mp\right]\mathrm{absolute\_term\_sum}.
```

The default `amplitude_tolerance` is $`128\,\epsilon`$. This allowance accounts
for the extent of phase powers and finite sums at a conservative numerical
scale; it is **not a derived interval error bound**. Delta and the nominal
phantom phase are fixed in the input-variation calculation. None of these
diagnostics rigorously bounds continuation error, deflation error, or
uncertainty in the requested coupling.

A witness is recorded only if its magnitude exceeds the sum of propagated
input variation and arithmetic allowance. The absolute-term scale must be
normal and the total allowance positive, so underflow to a zero allowance
cannot produce a spurious infinite-resolution success. This is stronger
than testing `amplitude != 0`, but remains a numerical state check, not a
proof about an exact Bethe vector or a sector minimum.

## Work budgets and outcomes

Ordered configurations are tried lexicographically until a witness is found,
the basis is exhausted, or a budget prevents further work. Defaults are
16 configurations and 1000000 total subset updates. Each configuration
costs `binomial(M,p)*2^r` updates, with O(r²) arithmetic per update; the
complete cost is checked before starting that configuration. A partial
amplitude is never accepted. The finite-degree and iteration limits in
`root_options` apply separately to root recovery.
Counters include attempted configurations and charge a complete subset
table for each attempted finite-amplitude evaluation, even if arithmetic
fails before that evaluation completes.

Statuses distinguish:

- `nonzero_witness`: a tested amplitude passes the resolution comparison.
- `reduction_unresolved`: the specified endpoint/phase/finite-state
  hypothesis did not pass; its detailed reduction result is retained.
- `roots_unresolved`: root recovery did not resolve distinct finite roots.
- `resolution_unresolved`: conversion or amplitude resolution is
  insufficient, including exhaustive searches with no resolved witness.
- `work_limit`: a degree, amplitude, or subset budget prevents completion.
- `nonfinite`: the reduction or amplitude arithmetic is nonfinite.

In particular, **none of the failure statuses means the vector is zero**.
The result retains the best tested configuration and its complex amplitude,
absolute-term sum, variation and arithmetic allowances, resolution ratio,
and work counters. `all_configurations_tested` distinguishes actual
exhaustion from a truncated search, without needing to form `binomial(N,M)`.

## Use in the internal sector scan

The odd-ring candidate scan tries this diagnostic only after the regular
Gaudin criterion and explicit helix match are unresolved. It tries mixed
multiplicities and both chiralities, retaining every attempted result in
`sector.phantom_lifts`. Configuration and subset budgets are shared across
these hypotheses **per sector**; root-iteration limits apply per hypothesis.
`phantom_work_limited` records whether any attempt was prevented by a budget.
The fifth scan argument supplies `PhantomLiftOptions`; a zero amplitude
budget disables this extra work.

`state_checks_complete` now includes numerical phantom witnesses alongside
regular states and helix matches. This flag is explicitly not a rigorous
state certificate. It also remains separate from equation convergence,
Wronskian consistency, and candidate-energy comparison. Disabling witnesses
does not change continuation coefficients, energies, Newton iteration
counts, the selected index, or the near-degeneracy report.

Tests cover continued mixed branches through N=13, both chiralities,
precisely exhausted budgets, deliberately unattainable resolution demands,
and rejected off-shell hypotheses. The variation recurrence is checked
against independently enumerated, perturbed permutation sums. Existing
direct-Hamiltonian tests remain independent checks of the constructed
vectors. All run in fp64, long double, and fp128. Public odd-ring
ground-state support still needs reliable sector-minimum tracking and
coverage of the remaining singular families.
