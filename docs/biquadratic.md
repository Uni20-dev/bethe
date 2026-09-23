# Free-end spin-1 biquadratic chains and Temperley–Lieb models

[Overview](../README.md) · [Model catalogue](models.md) · [Precision and CLI controls](command-line.md)

`bethe-biquadratic-obc` calculates the singlet ground state, TL module minima,
and restricted real-root excitations of the **even-length, free-end** spin-1
pure biquadratic chain, N>=2:

```text
H_b = -sum_(i=1)^(N-1) (S_i.S_(i+1))^2.
```

The coefficient is -1, the S operators are spin 1, and there is no added
constant. A bond has energies -4, -1, -1 in total spin 0, 1, 2 respectively.
This is neither the [Takhtajan–Babujian point](takhtajan-babujian.md) nor
the [ULS point](su3.md): both of those also have a bilinear term.

```sh
build/bethe-biquadratic-obc 64
build/bethe-biquadratic-obc 4 --roots --precision long-double
# With binary128 enabled:
build/bethe-biquadratic-obc 64 --precision fp128
```

Two useful checks are `E(2)=-4` and
`E(4)=-(15+sqrt(17))/2=-9.561552812808830274910704927...`.
The report includes physical, TL, and reference XXZ energies, convergence
diagnostics, and CPU time. `--roots` prints the auxiliary XXZ rapidities
and integer labels. There is no lattice momentum for these open chains.

## Excited levels and their multiplicities

```sh
# Lowest energy in each TL module (not each physical-spin sector):
build/bethe-biquadratic-obc 16 --sectors
build/bethe-biquadratic-obc 16 --through-lines 2
# Lowest 10 supported real-root levels in ell=2, including its minimum:
build/bethe-biquadratic-obc 16 --excitations 10
build/bethe-biquadratic-obc 16 --through-lines 4 --excitations all
# A specific Bethe level (its root count determines ell):
build/bethe-biquadratic-obc 8 --quantum-numbers 1,2,5 --roots
```

Each numerical level carries a TL through-line label ell and a **multiplicity
per TL eigenvector**. For even ell=0,2,4,6,... these multiplicities are
1,8,55,377,... in the spin-1 chain. In particular, the ell=2 multiplicity
space is one spin-1 triplet plus one spin-2 quintuplet, not eight multiplets
and not a single spin-1 multiplet. The code reports the total multiplicity;
it does not yet decompose general ell into physical SU(2) spins.
Coincident levels are not numerically merged: sum their weights if a genuine
degeneracy has been established. A decimal energy tolerance alone is not a
reliable criterion for doing so.

`--excitations` defaults to ell=2; use `--through-lines` to select any even
ell in [0,N]. Gaps are relative to the **global singlet ground state**, not
the selected module's minimum. COUNT limits retained levels, not work: every
candidate in the selected family is solved before energy ordering. The
`--max-candidates` guard (default 10000) checks the combinatorial count before
allocating roots or solving a reference state. Failed candidates are excluded,
reported, and cause exit status 2; a failed ground reference makes gaps
unavailable even if some candidate levels converge.

**`all` does not mean the complete spectrum.** For M=(N-ell)/2 roots the
family consists of increasing integer labels selected from 1,...,N-M,
giving `choose(N-M,M)` candidates. Complex-root levels are missing. For
example, N=4, ell=0 has two TL eigenvalues, but only its ground state is in
this real-root family. The other singlet energy `(-15+sqrt(17))/2` is not
returned. By contrast, the zero- and one-root modules are complete here:
for N=4, ell=2 the three energies are `-6-sqrt(2), -6, -6+sqrt(2)`, each
with multiplicity 8. The first has gap `2.147339250435735226109016203...`.

`--sectors`, `--excitations`, and `--quantum-numbers` are separate modes.
Explicit labels determine ell and cannot be combined with `--through-lines`;
use `--quantum-numbers none` for the zero-root ell=N level. Multiplicities
use checked uint64 arithmetic: ell>=46 at d=3 is reported as
`overflow (>uint64)`, never wrapped, rounded, or silently replaced by 1.
Energy calculations can still converge when this integer count is unavailable.

## The algebra connects spectra, not physical spin labels

Define `e_i=(S_i.S_(i+1))^2-1=3*P_(singlet,i,i+1)`. These operators satisfy
the Temperley–Lieb relations with loop weight lambda=3:

```text
e_i^2 = lambda*e_i,
e_i*e_(i+/-1)*e_i = e_i,
[e_i,e_j] = 0                         for |i-j|>1.
H_TL = -sum_i e_i,
H_b = H_TL - (N-1).
```

Another representation of the same open-chain algebra uses spin-1/2 XXZ
operators. For general lambda=2*Delta>2, our reference Hamiltonian is

```text
H_ref = sum_i [sx_i*sx_(i+1) + sy_i*sy_(i+1) + Delta*sz_i*sz_(i+1)]
        + sqrt(Delta^2-1)/2 * (sz_1-sz_N),
E_TL = 2*E_ref - (N-1)*lambda/4.
```

At lambda=3 this gives **Delta=3/2**, end-field coefficient sqrt(5)/4,
and `E_b=2*E_ref-7*(N-1)/4`. The end fields belong to the auxiliary
spin-1/2 representation; the physical spin-1 chain still has free ends.
The existing `bethe-xxz-obc` solves the **zero-end-field** XXZ model and
cannot be used with just an energy shift to obtain these results.

The equivalence is within TL modules. A module labelled by ell through-lines
has different multiplicities in different spin-chain representations.
For the generic singlet-projector spin-chain representation of local dimension d,

```text
m_0=1, m_1=d, m_(ell+1)=d*m_ell-m_(ell-1).
```

For spin 1, d=3 and the sequence is 1, 3, 8, 21, 55, ...; only ell with
the same parity as N occurs. For spin 1/2, d=2 and `m_ell=ell+1`.
Consequently an XXZ spin multiplet's degeneracy is not the biquadratic
degeneracy. Nor should ell/2 be relabelled as the physical spin-1 total spin:
the multiplicity space can contain several SU(2) multiplets.
An accidental degeneracy shared by different TL modules adds their contributions.
These formulas do not claim to cover arbitrary TL representations such as RSOS models.

The present ground state lies in the ell=0 module and is a single physical
singlet. In particular, **do not attach a factor of two for dimerization**
to a finite even open chain. Its ends select the dimer pattern.

## Bethe equations and numerical branch

Following [Albertini](../CITATIONS.md#albertini-2000), Eqs. (6)–(10), set
`Delta=cosh(eta)`. A module with ell through-lines uses M=(N-ell)/2 real
roots in (0,pi). Its lowest state has consecutive labels `I_i=i`,
i=1,...,M; other supported states leave holes in the label window:

```text
Theta(alpha;w) = 2*atan2(sin(alpha/2), tanh(w)*cos(alpha/2)),
2*N*Theta(alpha_i;eta/2)
  - sum_(j!=i) [Theta(alpha_i-alpha_j;eta)+Theta(alpha_i+alpha_j;eta)]
  = 2*pi*I_i,
E_ref = (N-1)*Delta/4 - sum_i (Delta^2-1)/(Delta-cos(alpha_i)).
```

Sending the largest root to pi gives `I=N-M+1`, but that endpoint has a
vanishing Bethe wavefunction and is excluded. Thus `1<=I_1<...<I_M<=N-M`.
The finite regular roots select quantum-group highest weights, so the
corresponding TL module is ell=N-2M: counting every auxiliary XXZ Sz sector
again would duplicate descendants. Its full module dimension is
`choose(N,M)-choose(N,M-1)`, generally larger than the real-family count.

The reflected sum must retain its phase branch when alpha_i+alpha_j>pi.
The solver uses `x_i=Theta(alpha_i;eta/2)/2` in (0,pi/2), with
`alpha_i=2*atan2(tanh(eta/2)*sin(x_i),cos(x_i))`. In this coordinate the
energy contribution is simply `-(Delta+cos(2*x_i))`, avoiding cancellation
in the rapidity energy formula near Delta=1. No complex boundary root is
needed for this reference model's even ground state.

Newton iteration uses an analytic Jacobian and a line search preserving
strictly ordered interior roots. The seed is the bare driving phase
`x_i=pi*I_i/(2*N)`, not an approximate thermodynamic result. Dense storage
is O(N^2), and the partial-pivot solve costs O(N^3) per update. It uses the
shared native-precision Newton helper with recoverable small-pivot failure,
without changing Uni20's process-wide error policy.

`residual_norm` is the maximum absolute logarithmic equation residual
divided by 2N. The default tolerance is 32 epsilon of the selected type;
it is **not an energy-error bound**. Failed solves retain mutually consistent
roots, energy, residual, and accepted-update count. A zero iteration budget
evaluates the seed. The CLI returns 2 for an unconverged result and 1 for
invalid input; it never silently changes precision or relaxes the tolerance.

## Library layers

```cpp
#include <bethe/biquadratic.hpp>
auto state = bethe::biquadratic::ground_state<long double>(64);
if (!state.reference.converged) { /* state.energy is an incomplete estimate */ }
auto sector = bethe::biquadratic::sector_ground_state<long double>(64, 2);
auto levels = bethe::biquadratic::real_excitations<long double>(16, 2, {.count=10});
if (!levels.converged()) { /* failed candidates and/or reference; inspect diagnostics */ }
// Each level has state.energy, state.through_lines, state.multiplicity and optional gap.

// The generic even, zero-through-line TL ground state, lambda>2:
auto tl = bethe::temperley_lieb::open_ground_state(64, 3.0L);
// Exact representation multiplicity, nullopt if uint64_t would overflow:
auto copies = bethe::temperley_lieb::spin_chain_multiplicity(3, 4); // 55
```

The three layers have separate responsibilities:

- [xxz_quantum_group.hpp](../include/bethe/xxz_quantum_group.hpp): the
  reference equations, roots, energy and numerical diagnostics, for Delta>1.
- [temperley_lieb.hpp](../include/bethe/temperley_lieb.hpp): the TL energy
  normalization and spin-chain representation multiplicities, for lambda>2.
- [biquadratic.hpp](../include/bethe/biquadratic.hpp): the physical spin-1
  Hamiltonian and its additive energy shift.

All numerical layers retain fp64, native long double, or enabled fp128.
The multiplicity helper uses checked integer arithmetic, not floating point.
`sector_ground_state`, `solve_real`, and `real_excitations` are available in
both the TL and physical biquadratic layers. TL calls additionally take a
loop weight; their energies do not include the physical spin-1 shift. The
generic TL layer leaves representation multiplicities to its caller.
The `GroundState`/`OpenGroundState` type aliases retain the original ground API.

## Validation and future slices

Tests build the physical spin-1 bond matrix and square it independently.
For N=2,4,6 they reconstruct its **entire ED spectrum** from independently
diagonalized XXZ sectors and TL multiplicities. They also verify the unique
ground state and explicitly detect the wrong answer when end fields are
omitted. This validates the spectral mapping; it does not make a full-spectrum
Bethe solver available. Further tests check the original multiplicative Bethe
equations, the analytic Jacobian, exact two-/four-site energies in each native
precision, longer chains through N=128, and CLI/failure contracts.
Excitation tests isolate TL modules by subtracting adjacent auxiliary Sz
spectra, and match every returned real-root level for N=2,4,...,10 without
reusing an ED eigenvalue. They verify module minima, eightfold physical
degeneracies, exact one-root energies and gaps at native precision, candidate
limits, failed-reference behavior, and the missing complex-root singlet.

Odd free-end chains need their own one-domain-wall/spinon branch. Their
low-lying states describe motion of that defect, not just a factor-two
choice of dimer pattern. Periodic chains require sector-dependent XXZ twists
and periodic representation bookkeeping. Neither is enabled here, nor are
complex-root excitations, physical-spin sector scans, general boundary fields,
or thermodynamics.

The original spectral mapping is due to
[Barber–Batchelor](../CITATIONS.md#barber-batchelor-1989). For the open-chain
normalization and real-root equations see Albertini above; for TL modules
and multiplicities see [Aufgebauer–Klümper](../CITATIONS.md#aufgebauer-klumper-2010),
Secs. 2.3 and 3. All three references also appear with `--references`.
