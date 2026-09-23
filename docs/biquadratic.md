# Free-end spin-1 biquadratic chains and Temperley–Lieb models

[Overview](../README.md) · [Model catalogue](models.md) · [Precision and CLI controls](command-line.md)

`bethe-biquadratic-obc` calculates the unique singlet ground state of the
**even-length, free-end** spin-1 pure biquadratic chain, N>=2:

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
`Delta=cosh(eta)`. The even ground state has M=N/2 real roots in (0,pi)
and consecutive labels `I_i=i`, i=1,...,M:

```text
Theta(alpha;w) = 2*atan2(sin(alpha/2), tanh(w)*cos(alpha/2)),
2*N*Theta(alpha_i;eta/2)
  - sum_(j!=i) [Theta(alpha_i-alpha_j;eta)+Theta(alpha_i+alpha_j;eta)]
  = 2*pi*I_i,
E_ref = (N-1)*Delta/4 - sum_i (Delta^2-1)/(Delta-cos(alpha_i)).
```

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
It supports general ell, but the numerical solver currently supplies only
the **even ell=0 ground state**, not the other modules' energies.

## Validation and future slices

Tests build the physical spin-1 bond matrix and square it independently.
For N=2,4,6 they reconstruct its **entire ED spectrum** from independently
diagonalized XXZ sectors and TL multiplicities. They also verify the unique
ground state and explicitly detect the wrong answer when end fields are
omitted. This validates the spectral mapping; it does not make a full-spectrum
Bethe solver available. Further tests check the original multiplicative Bethe
equations, the analytic Jacobian, exact two-/four-site energies in each native
precision, longer chains through N=128, and CLI/failure contracts.

Odd free-end chains need their own one-domain-wall/spinon branch. Their
low-lying states describe motion of that defect, not just a factor-two
choice of dimer pattern. Periodic chains require sector-dependent XXZ twists
and periodic representation bookkeeping. Neither is enabled here, nor are
excitations, physical-spin sector scans, general boundary fields, or thermodynamics.

The original spectral mapping is due to
[Barber–Batchelor](../CITATIONS.md#barber-batchelor-1989). For the open-chain
normalization and real-root equations see Albertini above; for TL modules
and multiplicities see [Aufgebauer–Klümper](../CITATIONS.md#aufgebauer-klumper-2010),
Secs. 2.3 and 3. All three references also appear in executable help.
