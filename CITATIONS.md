# References and provenance

Use any executable's `--references` option for its relevant literature and
applicability. Ordinary help and no-argument usage omit the bibliography.
Bibliographic metadata and tool selections live in
[data/citations.json](data/citations.json); the [bibliography](#bibliography)
below and the C++ registry are generated from it. The provenance notes remain
hand-written. See [maintaining citations](docs/citations.md).

The [paper archive](papers/README.md) provides selected reference PDFs with
edition, source, license, and checksum records. Other papers remain linked;
being cited here does not grant permission to redistribute a copy.

The bibliography also contains literature for the
[model catalogue and proposals](docs/models.md). Those references do not
imply implemented methods and are not selected by executable help unless
the corresponding tool actually supports the cited calculation.

The Heisenberg ground-state solver succeeds Ian McCulloch's
`misc/heisenberg-energy.cpp` in the Matrix Product Toolkit (copyright
2015-2023). The original is maintained there, rather than duplicated here.
The successor retains its GPL-3.0-or-later notice.

The logarithmic equations, quantum numbers, energy convention, and real-root
iteration follow [karbach-1998](#karbach-1998), Eqs. (6)-(9).
Table I supplies the N=16 regression reference, including its N/4 energy shift.

The sector-state selection extends the consecutive quantum-number construction
in that paper (Eq. (16)). For odd-chain one-spinon states and hole labeling, see
[groha-2017](#groha-2017), Sec. 3. The present
solver concerns the unperturbed integrable chain, not the decay calculation.

The zero-field thermodynamic XXX and gapless XXZ dispersions follow
the notes [caux-xxx-spinons](#caux-xxx-spinons) and
[caux-xxz-spinons](#caux-xxz-spinons). The thermodynamic XXZ curve is currently
a library facility, not a mode of the finite-chain XXZ executable. We use positive spinon
momentum k in [0,pi], reversing the sign of the notes' convention. At finite
odd N our explicit convention is k=pi/2-2*pi*I_h/N, related to lattice momentum
by P=pi*M+pi/2-k modulo 2*pi. The bulk reference energy is J*(1/4-log(2)).

Research using these tools should acknowledge the Matrix Product Toolkit and
Uni20 where relevant, and cite the methods used in the calculation.

The free-end spin-1 pure biquadratic solver implements the spectral mapping
of [barber-batchelor-1989](#barber-batchelor-1989). Its reference XXZ equations
and normalization follow [albertini-2000](#albertini-2000), Eqs. (2)–(10): our
reference Hamiltonian is half Eq. (3), after a staggered rotation of the xy
exchange, so `E_b=2*E_ref-7*(N-1)/4`. The opposite XXZ end fields are essential;
the physical spin-1 model has free ends. We implement even-chain TL module
minima and real-root label scans in `1<=I<=N-M`, including the global singlet,
not complex-root levels or the odd-chain spinon band. Finite regular Bethe
states are assigned ell=N-2M without counting auxiliary descendants again.
The TL module/multiplicity
separation follows [aufgebauer-klumper-2010](#aufgebauer-klumper-2010),
Secs. 2.3 and 3, specialized to generic singlet-projector spin-chain
representations. It does not implement periodic twists, full-spectrum Bethe
enumeration, physical-spin decomposition, or thermodynamics. See the
[biquadratic/TL guide](docs/biquadratic.md) for the explicit conventions.

The ferromagnetic option reverses the physical Hamiltonian, not the auxiliary
XXZ equations: `H=+sum(S.S)^2=(N-1)+sum e_i`. The one-defect module has
`E-E0=3+2*cos(pi*j/N)`, `j=1,...,N-1`. The exact OBC gap above the entire
ground space, `3-2*cos(pi/N)`, follows from
[koma-nachtergaele-1997](#koma-nachtergaele-1997), Proposition 2 and Eq. (3.30),
using TL equivalence and rescaling their normalized XXZ Hamiltonian by
`2*Delta=3`. The sign convention and ground-space degeneracies are also
discussed in [zhou-2025-biquadratic](#zhou-2025-biquadratic).
Real-family selection reverses energy ordering before truncation. Q-system
discoveries are reordered in physical energy without strengthening their
completeness claim. See the [ferromagnetic excitation guide](docs/biquadratic-ferromagnetic.md).

The [targeted bound-pair solver](docs/biquadratic-bound-pairs.md) uses an
empty-sea two-string regularization of [bajnok-2020](#bajnok-2020), Eq. (5.12),
with both signs of the finite deviation. Its phase label and deviation sign
must satisfy N-J odd for positive deviation, even for negative deviation.
The thermodynamic two-defect threshold 5/3 follows from
[nachtergaele-spitzer-starr-2007](#nachtergaele-spitzer-starr-2007), Theorem 2.1,
after the same factor-three rescaling. These are targeted finite-chain modes,
not an ideal-string substitution, full two-defect spectrum or form-factor calculation.

The free-end XXX equations follow the rational limit of Eqs. (11)-(12) and
footnote 2 in [mei-2017](#mei-2017).
Our rapidity is `z=2*lambda`; our spin-1/2 Hamiltonian is the paper's XXX
Hamiltonian divided by four and shifted by `N/4`. This gives the ferromagnetic
reference `(N-1)/4`. Taking the positive-root logarithmic branch gives the
integer quantum-number convention documented in the
[open-chain guide](docs/open-chains.md#allowed-states-and-convergence).

For the highest-weight multiplet and string classification relevant to
real-root excitation scans, see [caux-xxx-states](#caux-xxx-states).
The scans enumerate only the supported all-real quantum-number window;
they do not implement the string sectors described in those references.

The finite-size periodic XXZ sector solver follows the logarithmic equations,
energy and momentum in [vlijm-2016](#vlijm-2016), Eqs. (1), (3)-(5).
We set J=1, h=0 and add N*Delta/4 to the paper's shifted Hamiltonian.
Our coordinate z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma), rewrites the
equations algebraically into the rational form documented in the
[XXZ guide](docs/xxz.md#scaled-rapidities-and-bethe-equations);
it has the smooth XXX limit z=2*lambda_XXX. The excitation scan's one-rapidity
infinity bound is derived by taking lambda_j to infinity in Eq. (3):
theta_1 tends to pi-gamma and theta_2 to pi-2*gamma, giving
I_infinity=[N-M+1-(N-2*M+2)*gamma/pi]/2. We intersect the strict bound with
the conventional XXX window and exclude a precision-dependent rounding band
at the infinity threshold. This deliberately restricted family is not a
complete classification of XXZ solutions. We do not implement the paper's
correlation functions.

For periodic XXZ at Delta>1, the trigonometric equations follow
[dugave-2015](#dugave-2015), Eqs. (1.1)-(1.2), at zero twist and field.
Our spin operators are sigma/2, so our J=1 energy is one quarter of the
paper's J=1 Pauli energy. We rewrite real rapidities as
`z=tan(lambda)/tanh(eta/2)`, Delta=cosh(eta), retaining the scattering phase's
winding with atan2. The consecutive sector-minimum labels are checked against
finite-chain energy and momentum ED, including odd lengths. This does not
extend the existing excitation window or implement the paper's form factors.

The internal negative-Delta ground-root engine uses the periodic equations
in [kozlowski-2017](#kozlowski-2017), Eqs. (0.4), (0.7), and the open reflection
equations of [mei-2017](#mei-2017). Kozlowski explicitly assumes even rings;
we do not extend that ground-state identification to odd periodic chains.
The [implementation note](docs/xxz-negative.md) derives the rank-subtracted,
scaled equations and records the unresolved odd-ring and public-API work.
This internal checkpoint does not extend the executable help's scope.

The internal [polynomial XXZ equations](docs/xxz-polynomial.md) start from
[caux-xxz-coordinate](#caux-xxz-coordinate), Eqs. `xxz.be` and `xxz.e`,
at zero twist and J=1, restoring the ferromagnetic energy N*Delta/4.
The divided-Horner recurrence, analytic coefficient derivatives, and collision
regularization are our algebraic implementation. This is a tested equation
building block, not a complete odd-ring ground-state or complex-string solver.

The internal [polynomial root recovery](docs/polynomial-roots.md) uses
Ehrlich-Aberth iteration as reviewed in [bini-robol-2013](#bini-robol-2013),
equation (3), with Gauss-Seidel updates. Native-real error-free transforms
and compensated complex Horner evaluation follow the numerical strategy
discussed in [cameron-graillat-2022](#cameron-graillat-2022), Sections 2-4.
Our fixed-precision implementation and diagnostic stopping rules are not
MPSolve or a reproduction of the latter paper's validated error bounds.
It does not perform arbitrary-precision escalation or interval certification.

The massive free-end XXZ ground-state module also follows the reflection
equations in [mei-2017](#mei-2017). Its distinguished boundary root and
finite-size deviation are motivated by [grijalva-2019](#grijalva-2019),
Sec. 4.3.2; the [module guide](docs/xxz-open-massive.md) derives the implemented
regularization. This is a ground-state solver, not a boundary correlation
function, general string solver, or massive excitation scan. The existing
`bethe-xxz-obc` frontend selects it for Delta>1 ground states and sector minima.

The free-end XXZ solver uses [mei-2017](#mei-2017), Eqs. (11)-(12),
with the paper's Hamiltonian divided by four and shifted by `N*Delta/4`.
The resulting ferromagnetic reference is `(N-1)*Delta/4`. In the massless
regime, the cosh ratio in Eq. (11) supplies the boundary reflection phase;
in our scaled coordinate it contributes `4*atan(c*z)`,
`c=(1-Delta)/(1+Delta)`. This term remains nonzero at the XX point, producing
the open-chain standing-wave denominator N+1. Taking one rapidity to infinity
gives `I_infinity=N-M+1-(N-2*M+1)*gamma/pi`. As in the periodic case, we intersect
this strict bound with the conventional XXX window and exclude a rounding
band at the threshold; this is a restricted positive finite-root family.
See the [free-end XXZ guide](docs/xxz-open.md#boundary-equations-and-limiting-cases)
for the transformed equations and supported excitation window.

The periodic Hubbard solver uses [lieb-wu-2003](#lieb-wu-2003),
Eqs. (1), (11), and (14)-(18): unshifted Hamiltonian, energy, logarithmic
equations, quantum-number parity, and the consecutive ground-state labels.
We set t=1, use u=U/4, and write the scattering phases as positive `2*atan`
terms with the signs shown in the [Hubbard guide](docs/hubbard.md).
The paper's Eq. (18) supplies the centered labels for even N and odd M,
including our supported doped sectors. At half filling with even M we use
its parity rules and a full charge Brillouin zone with the +pi endpoint,
together with centered spin labels. This extension is checked against small-ring fermionic
exact diagonalization, including momentum; Eq. (18) alone is not a statement
of both parity branches. The reflection reduction, analytic Jacobian, and
damped-Newton continuation are our numerical implementation. We implement
neither the full excited-state string classification nor correlation functions.

The Hubbard full and partial particle-hole transformations are reviewed in
[rylands-2022](#rylands-2022), Sec. II, Eqs. (4)-(7). We apply them to fixed
particle/spin sectors on even rings and open chains of either length parity.
The unshifted-Hamiltonian energy offsets
are derived from `n_down -> 1-n_down`; momentum offsets follow from the
complementary occupied momenta under `k -> pi-k`. Both are independently
tested against fermionic exact diagonalization. See the
[sector guide](docs/hubbard-sectors.md) for the explicit conventions and the
distinction between physical observables and auxiliary repulsive roots.
We do not implement that paper's quench dynamics or overlap formulas.

The free-end Hubbard equations follow [deguchi-yue-1997](#deguchi-yue-1997),
Eqs. (2.5)-(2.8) and (3.1)-(3.2), with all boundary potentials and fields set
to zero. The charge boundary factor supplies the length L+1; both direct
and reflected spin self-scattering terms are excluded. We use positive
roots and consecutive integer labels, as documented in the
[open Hubbard guide](docs/hubbard-open.md#reflected-scattering-and-ground-state-labels).
The one-dimensional open-chain spin ordering of
[lieb-mattis-1962](#lieb-mattis-1962) underlies selecting the minimum-spin
branch in each fixed-Sz sector, without the periodic shell-parity restrictions.
The Jacobian, continuation scheme, and symmetry-sector bookkeeping are our
implementation, tested against independent fermionic exact diagonalization.
We implement neither boundary fields nor the papers' conformal-dimension analysis.

The repulsive periodic Bose gas follows [lieb-liniger-1963](#lieb-liniger-1963).
The explicit Hamiltonian, logarithmic equations, label parity, energy and
momentum were cross-checked against [essler-de-klerk-2023](#essler-de-klerk-2023),
Eqs. (4), (24)-(27). We set hbar^2/(2m)=1 and retain the paper's coupling 2c.
Its root-density relation (32)-(33), restricted to a filled symmetric interval,
provides an independent test-only bulk benchmark. The complementary weak-coupling
equation, scaled residual, Newton implementation, and finite-window scan are
our numerical choices; see the [model guide](docs/lieb-liniger.md).
[lieb-1963-excitations](#lieb-1963-excitations) supplies excitation background,
not an implemented thermodynamic dispersion mode. No matrix elements or
finite-temperature thermodynamics from the cited papers are implemented.

The periodic SU(3) permutation chain follows
[sutherland-1975](#sutherland-1975), with explicit nested equations, labels,
and energy checked against [doikou-nepomechie-1998](#doikou-nepomechie-1998),
Eqs. (2.17)-(2.19) and (2.24)-(2.29). We implement only the all-real,
filled-sea balanced singlet of Sec. 2.3 for L divisible by three. The reference
uses `sum(P-1)/2`; our `H=sum P` gives `E=2*E_paper+L`.
The density-inversion seed, reflection reduction, Jacobian, and damped Newton
solve are our implementation. The [SU(3) guide](docs/su3.md) derives the
six-site singlet oracle and the spin-1 ULS energy shift, independently checked
in color/spin space. No excited strings, scattering matrices, boundary fields,
or general-rank SU(n) solver are implied by citing this paper.

The periodic Gaudin–Yang gas follows [gaudin-1967](#gaudin-1967) and
[yang-1967](#yang-1967), with the explicit Hamiltonian and nested equations
from [oelkers-2006](#oelkers-2006), Eqs. (1), (2), (7), (8), (25).
We implement the repulsive even-N, odd-minority filled-sea branch, plus exact
free gases, using the paper's kinetic units and coupling 2c. Its weak/strong
limits (12), (15)-(16) supplement the independent two-body jump condition
and dilute-Hubbard discretization tests. Reflection reduction, scaled spin
coordinates, continuation, complementary phases, and residual normalization
are our numerical choices; see the [Gaudin–Yang guide](docs/gaudin-yang.md).
Other periodic shell branches, attraction, hard walls, excited states, and
thermodynamic integral equations are not implemented by this tool.

<!-- BEGIN GENERATED BIBLIOGRAPHY -->

## Bibliography

Generated from [data/citations.json](data/citations.json); edit the registry, not this section.

### barber-batchelor-1989

Michael N. Barber and Murray T. Batchelor. *Spectrum of the biquadratic spin-1 antiferromagnetic chain*.
Phys. Rev. B 40, 4621-4626 (1989).

[DOI](<https://doi.org/10.1103/PhysRevB.40.4621>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Original free-end spin-1 biquadratic/TL spectral correspondence. Even-chain ground states, TL module minima and restricted real-root excitations of H=-sum (S.S)^2, with representation multiplicities.

### albertini-2000

Giuseppe Albertini. *Is the purely biquadratic spin 1 chain always massive?*.
arXiv:cond-mat/0012439 (2000).

[arXiv](<https://arxiv.org/abs/cond-mat/0012439>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Hamiltonians and energy shift (2)-(5), real-root Bethe equations (6)-(10), and integer-label window 1\<=I\<=N-M. Our spin-half reference is half of (3), after a staggered rotation. The odd-chain spinon band is not implemented.

### aufgebauer-klumper-2010

Britta Aufgebauer and Andreas Klümper. *Quantum spin chains of Temperley-Lieb type: periodic boundary conditions, spectral multiplicities and finite temperature*.
J. Stat. Mech. 2010, P05018 (2010).

[arXiv](<https://arxiv.org/abs/1003.1932>).

Relevant tool modes:

- `bethe-biquadratic-obc`: TL loop weight and quantum-group XXZ end fields, Sec. 2.3; open spin-chain module multiplicities, Sec. 3, especially (48)-(55). Periodic twists, physical-spin decomposition and thermodynamics are not implemented.

### essler-korepin-1999

Fabian H. L. Essler and Vladimir E. Korepin. *Form factors in the half-filled Hubbard model*.
Phys. Rev. B 59, 1734-1738 (1999).

[DOI](<https://doi.org/10.1103/PhysRevB.59.1734>), [arXiv](<https://arxiv.org/abs/cond-mat/9808018>).

Relevant tool modes:

- `bethe-hubbard-dispersion`: Half-filled spinon and holon dispersions, Eqs. (4)-(5), with our U equal to four times the paper's U; no form factors or spectral weights are implemented.

### melzer-1995

Ezer Melzer. *On the scaling limit of the 1D Hubbard model at half filling*.
Nucl. Phys. B 443, 553-564 (1995).

[DOI](<https://doi.org/10.1016/0550-3213(95)00105-2>), [arXiv](<https://arxiv.org/abs/cond-mat/9410043>).

Relevant tool modes:

- `bethe-hubbard-dispersion`: Exact massive-branch K-Bessel series, resummed here into nonoscillatory positive integrals to retain the weak-coupling charge gap; hopping and interaction conventions are converted as documented.

### essler-2010

Fabian H. L. Essler. *Threshold singularities in the one-dimensional Hubbard model*.
Phys. Rev. B 81, 205120 (2010).

[DOI](<https://doi.org/10.1103/PhysRevB.81.205120>), [arXiv](<https://arxiv.org/abs/1002.1671>).

Relevant tool modes:

- `bethe-hubbard-dispersion`: Doped zero-field density and dressed-energy equations (103)-(106), with dressed momenta from (28)-(29); our mu\_unshifted equals mu+2u in (105), u=U/4. Elementary lines only, not threshold exponents or continuum minimization.

### luo-pu-guan-2024

Jia-Jia Luo, Han Pu, and Xi-Wen Guan. *Exact results of the one-dimensional repulsive Hubbard model*.
Rep. Prog. Phys. 87, 117601 (2024).

[DOI](<https://doi.org/10.1088/1361-6633/ad7b70>), [arXiv](<https://arxiv.org/abs/2307.00890>).

Relevant tool modes:

- `bethe-hubbard-dispersion`: Real charge particle/hole and spinon excitation interpretation, Sec. II; doped charge-particle is a real-root addition, not the gapped half-filled antiholon or a k-Lambda string.

### karbach-1998

Michael Karbach, Kun Hu, and Gerhard Müller. *Introduction to the Bethe ansatz II*.
Computers in Physics 12, 565 (1998).

[DOI](<https://doi.org/10.1063/1.168740>), [arXiv](<https://arxiv.org/abs/cond-mat/9809163>).

Relevant tool modes:

- `bethe-xxx-pbc`: Periodic XXX equations, energy normalization, and sector quantum numbers; Eqs. (6)-(9), (16), Table I.
- `bethe-xxz-pbc`: XXX limit at Delta=1 and the conventional real-root quantum-number window.
- `bethe-tj-pbc`: No-hole reduction to the periodic XXX sector solver: H\_tJ=2\*H\_XXX-L/2. Fermionic translation adds the filled-reference phase (-1)^(L-1).

### groha-2017

Stefan Groha and Fabian H. L. Essler. *Spinon decay in the spin-1/2 Heisenberg chain with weak next nearest neighbour exchange*.
J. Phys. A 50, 334002 (2017).

[arXiv](<https://arxiv.org/abs/1702.06550>).

Relevant tool modes:

- `bethe-xxx-pbc`: Odd-chain one-spinon states and hole labels, Sec. 3; only the unperturbed integrable chain is used.

### caux-xxx-spinons

Jean-Sébastien Caux. *The Bethe Ansatz: XXX spinons*.
Online notes.

[XXX spinons](<https://integrability.org/g_h_e.html>).

Relevant tool modes:

- `bethe-xxx-pbc`: Thermodynamic XXX spinon dispersion and the bulk reference for --spinons.

### caux-xxz-spinons

Jean-Sébastien Caux. *The Bethe Ansatz: XXZ spinons*.
Online notes.

[XXZ spinons](<https://integrability.org/g_sc_p_e.html>).

### mei-2017

Zhongtao Mei and C. J. Bolech. *Derivation of matrix product states for the Heisenberg spin chain with open boundary conditions*.
Phys. Rev. E 95, 032127 (2017).

[arXiv](<https://arxiv.org/abs/1609.08045>).

Relevant tool modes:

- `bethe-xxx-obc`: Free-end XXX equations: rational limit of Eqs. (11)-(12) and footnote 2; our spin-1/2 normalization differs.
- `bethe-xxz-obc`: Free-end XXZ equations and boundary reflection phase, Eqs. (11)-(12); our Hamiltonian is divided by four and shifted. Negative-Delta ground states use an algebraically rank-subtracted, rescaled form of these equations.

### caux-xxx-states

Jean-Sébastien Caux. *The Bethe Ansatz: real rapidities, SU(2) descendants, and two-string states*.
Online notes.

[Real rapidities](<https://integrability.org/c_h_e_rr.html>), [SU(2) descendants](<https://integrability.org/c_h_e_rr_10.html>), [Two-string states](<https://integrability.org/c_h_e_s_2.html>).

Relevant tool modes:

- `bethe-xxx-pbc`: Highest-weight and string-classification background; strings and infinite-root descendants are not implemented.
- `bethe-xxx-obc`: SU(2) and string-classification background, not the boundary equations; scans include only the supported finite-real family.
- `bethe-xxz-pbc`: Background for the XXX window used to restrict the XXZ scan; not a complete XXZ state classification.

### vlijm-2016

R. Vlijm, I. S. Eliëns, and J.-S. Caux. *Correlations of zero-entropy critical states in the XXZ model: integrability and Luttinger theory far from the ground state*.
SciPost Phys. 1, 008 (2016).

[DOI](<https://doi.org/10.21468/SciPostPhys.1.1.008>), [arXiv](<https://arxiv.org/abs/1606.09516>).

Relevant tool modes:

- `bethe-xxz-pbc`: Massless periodic XXZ equations, energy, and momentum, Eqs. (1), (3)-(5); our energy includes the N\*Delta/4 shift.

### bajnok-2020

Zoltán Bajnok, Etienne Granet, Jesper Lykke Jacobsen, and Rafael I. Nepomechie. *On Generalized Q-systems*.
JHEP 03 (2020) 177 (2020).

[arXiv v2](<https://arxiv.org/abs/1910.07805v2>), [HTML](<https://arxiv.org/html/1910.07805v2>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Open quantum-group-invariant XXZ Q-system, Sec. 5: Wronskian (5.17), Bethe equations (5.12), and admissibility conditions. Used for Q-system module searches (validated through N=8; larger sizes experimental), the regularized AF two-string singlet, and empty-sea two-/three-string ferro branches on odd/even long chains. Original equations also underlie odd/even selected real roots and our restricted high-label scattering windows. Two-strings retain signed real deviations; three-strings retain complex deviations. Neither numerical completeness nor long-chain energy ordering is rigorously certified.

### popkov-2021

Vladislav Popkov, Xin Zhang, and Andreas Klümper. *Phantom Bethe excitations and spin helix eigenstates in integrable periodic and open spin chains*.
Phys. Rev. B 104, L081410 (2021).

[DOI](<https://doi.org/10.1103/PhysRevB.104.L081410>), [arXiv v2](<https://arxiv.org/abs/2102.03295v2>), [HTML](<https://arxiv.org/html/2102.03295v2>).

### kozlowski-2017

Karol K. Kozlowski. *On condensation properties of Bethe roots associated with the XXZ chain*.
arXiv:1508.05741v2 (2017).

[arXiv v2](<https://arxiv.org/abs/1508.05741v2>), [HTML](<https://arxiv.org/html/1508.05741v2>).

Relevant tool modes:

- `bethe-xxz-pbc`: Even-ring negative-anisotropy sector ground states: Eqs. (0.4), (0.7), with Pauli exchange divided by four. We use rank-subtracted equations scaled near Delta=-1; the reference explicitly assumes even length and does not justify odd-ring ground-state selection.

### lieb-wu-2003

Elliott H. Lieb and F. Y. Wu. *The one-dimensional Hubbard model: A reminiscence*.
Physica A 321, 1-27 (2003).

[arXiv v2](<https://arxiv.org/abs/cond-mat/0207529v2>).

Relevant tool modes:

- `bethe-hubbard-pbc`: Lieb-Wu ground-state equations and quantum-number parity, Eqs. (1), (11), (14)-(18); t=1 and unshifted U\*n\_up\*n\_down.

### rylands-2022

Colin Rylands, Bruno Bertini, and Pasquale Calabrese. *Integrable quenches in the Hubbard model*.
J. Stat. Mech. 2022, 103103 (2022).

[DOI](<https://doi.org/10.1088/1742-5468/ac98be>), [arXiv](<https://arxiv.org/abs/2206.07985>).

Relevant tool modes:

- `bethe-hubbard-pbc`: Full particle-hole and partial particle-hole (Shiba) transformations, Sec. II, Eqs. (4)-(7); used for sector mappings, not quench dynamics.
- `bethe-hubbard-obc`: Bipartite particle-hole and Shiba transformations, Sec. II; used for physical-to-auxiliary sector mappings, not quench dynamics.

### deguchi-yue-1997

Tetsuo Deguchi and Ruihong Yue. *Exact solutions of 1-D Hubbard model with open boundary conditions and the conformal dimensions under boundary magnetic fields*.
arXiv:cond-mat/9704138 (1997).

[arXiv v2](<https://arxiv.org/abs/cond-mat/9704138v2>).

Relevant tool modes:

- `bethe-hubbard-obc`: Free-end Hubbard energy and nested reflection equations, Eqs. (2.5)-(2.8), (3.1)-(3.2), with all boundary potentials zero; t=1 and unshifted interaction.

### lieb-mattis-1962

Elliott Lieb and Daniel Mattis. *Theory of Ferromagnetism and the Ordering of Electronic Energy Levels*.
Phys. Rev. 125, 164-172 (1962).

[DOI](<https://doi.org/10.1103/PhysRev.125.164>).

Relevant tool modes:

- `bethe-hubbard-obc`: Spin ordering for the nearest-neighbour open chain: the sector minimum has S=|Sz|. No periodic-shell spin-branch selection is required.

### lieb-liniger-1963

Elliott H. Lieb and Werner Liniger. *Exact Analysis of an Interacting Bose Gas. I. The General Solution and the Ground State*.
Phys. Rev. 130, 1605 (1963).

[DOI](<https://doi.org/10.1103/PhysRev.130.1605>).

Relevant tool modes:

- `bethe-lieb-liniger-dispersion`: Repulsive zero-temperature root-density equation and bulk ground-state energy, in units hbar=2m=1.
- `bethe-lieb-liniger-pbc`: Repulsive Bose gas on a ring, H=-sum d\_j^2+2c sum delta, and its ground state; c\>0.

### lieb-1963-excitations

Elliott H. Lieb. *Exact Analysis of an Interacting Bose Gas. II. The Excitation Spectrum*.
Phys. Rev. 130, 1616 (1963).

[DOI](<https://doi.org/10.1103/PhysRev.130.1616>).

Relevant tool modes:

- `bethe-lieb-liniger-dispersion`: Thermodynamic type-I particle and type-II hole excitation branches; energies are fixed-N gaps relative to the ground state. No spectral weights or finite-temperature TBA.
- `bethe-lieb-liniger-pbc`: Excited-state background; we solve finite-volume states in explicit label windows, not a thermodynamic dispersion calculation.

### gaudin-1971

M. Gaudin. *Boundary Energy of a Bose Gas in One Dimension*.
Phys. Rev. A 4, 386 (1971).

[DOI](<https://doi.org/10.1103/PhysRevA.4.386>).

Relevant tool modes:

- `bethe-lieb-liniger-obc`: Repulsive Bose gas with Dirichlet walls, positive real roots and reflected scattering. We implement finite-volume states and bounded label windows, not general boundary potentials.

### gaudin-1967

M. Gaudin. *Un système à une dimension de fermions en interaction*.
Phys. Lett. A 24, 55-56 (1967).

[DOI](<https://doi.org/10.1016/0375-9601(67)90193-4>).

Relevant tool modes:

- `bethe-gaudin-yang-pbc`: Original spin-1/2 continuum fermion solution; we implement repulsive periodic ground states in selected sectors, not attraction.

### yang-1967

C. N. Yang. *Some Exact Results for the Many-Body Problem in One Dimension with Repulsive Delta-Function Interaction*.
Phys. Rev. Lett. 19, 1312 (1967).

[DOI](<https://doi.org/10.1103/PhysRevLett.19.1312>).

Relevant tool modes:

- `bethe-gaudin-yang-pbc`: Original multicomponent delta-gas solution; our implementation has two spin components only.

### yang-1968

C. N. Yang. *S Matrix for the One-Dimensional N-Body Problem with Repulsive or Attractive delta-Function Interaction*.
Phys. Rev. 168, 1920 (1968).

[DOI](<https://doi.org/10.1103/PhysRev.168.1920>).

### sutherland-1975

Bill Sutherland. *Model for a multicomponent quantum system*.
Phys. Rev. B 12, 3795 (1975).

[DOI](<https://doi.org/10.1103/PhysRevB.12.3795>).

Relevant tool modes:

- `bethe-su3-pbc`: Original multicomponent permutation-chain solution; the implementation selects only the fundamental SU(3) periodic balanced ground state.

### essler-korepin-1992

Fabian H. L. Essler and Vladimir E. Korepin. *Higher conservation laws and algebraic Bethe Ansätze for the supersymmetric t-J model*.
Phys. Rev. B 46, 9147 (1992).

[DOI](<https://doi.org/10.1103/PhysRevB.46.9147>), [arXiv](<https://arxiv.org/abs/hep-th/9207007>).

Relevant tool modes:

- `bethe-tj-pbc`: Projected t-J Hamiltonian (1.3)-(1.5), Sutherland BFF equations (3.73), and energy (3.75). We remove the shift 2\*N\_e-L: E=2\*N\_h-sum 1/(lambda^2+1/4). Doped mixed-spin coverage is restricted to odd N\_up and N\_down; J=2t=2.

### babujian-1982

H. M. Babujian. *Exact solution of the one-dimensional isotropic Heisenberg chain with arbitrary spins S*.
Phys. Lett. A 90, 479-482 (1982).

[DOI](<https://doi.org/10.1016/0375-9601(82)90403-0>).

Relevant tool modes:

- `bethe-tb-pbc`: Original integrable higher-spin chain family; only the periodic even-length spin-1 singlet ground state is implemented.

### vlijm-caux-2014

Rogier Vlijm and Jean-Sébastien Caux. *Computation of dynamical correlation functions of the spin-1 Babujan-Takhtajan chain*.
J. Stat. Mech. 2014, P05009 (2014).

[DOI](<https://doi.org/10.1088/1742-5468/2014/05/P05009>), [arXiv](<https://arxiv.org/abs/1401.4450>).

Relevant tool modes:

- `bethe-tb-pbc`: Spin-1 Hamiltonian and complex Bethe equations (1.2)-(1.4), filled two-string sea, and finite-deviation equations (3.8)-(3.11). Our bilinear coefficient is 1, i.e. J=4 in this paper. Equation (3.18) supplies only an initial guess; finite-size deviations are solved, not dropped. No dynamical correlations or excitations are implemented.

### dukelsky-2004

J. Dukelsky, S. Pittel, and G. Sierra. *Colloquium: Exactly solvable Richardson-Gaudin models for many-body quantum systems*.
Rev. Mod. Phys. 76, 643-662 (2004).

[DOI](<https://doi.org/10.1103/RevModPhys.76.643>), [arXiv](<https://arxiv.org/abs/nucl-th/0405011>).

Relevant tool modes:

- `bethe-richardson`: Reduced BCS pairing model and blocked-level sectors. We support attractive uniform pairing with distinct doubly degenerate single-particle levels, not arbitrary degeneracies or the full Gaudin family.

### faribault-2011

Alexandre Faribault, Omar El Araby, Christoph Sträter, and Vladimir Gritsev. *Gaudin models solver based on the correspondence between Bethe ansatz and ordinary differential equations*.
Phys. Rev. B 83, 235124 (2011).

[DOI](<https://doi.org/10.1103/PhysRevB.83.235124>), [arXiv](<https://arxiv.org/abs/1103.0472>).

Relevant tool modes:

- `bethe-richardson`: Eigenvalue variables (7), quadratic equations in Sec. II, and Richardson equations (32). We use e\_i=2\*epsilon\_i and y\_i=g\*Lambda(e\_i), continue the filled lowest levels from g=0, and include diagonal pair scattering. Energy is sum e\_i\*y\_i-g\*M\*(L-M+1), plus blocked single-particle energies. No pair-root reconstruction or form factors are implemented.
- `bethe-central-spin`: Quadratic eigenvalue-variable method and central-spin realization, Secs. II and IV.3. We compactify the inverse-field continuation to reach B=0; the seed minimizes the central-spin Hamiltonian, not the Richardson energy.

### faribault-schuricht-2013

Alexandre Faribault and Dirk Schuricht. *Spin decoherence due to a randomly fluctuating spin bath*.
Phys. Rev. B 88, 085323 (2013).

[DOI](<https://doi.org/10.1103/PhysRevB.88.085323>), [arXiv](<https://arxiv.org/abs/1306.2541>).

Relevant tool modes:

- `bethe-central-spin`: Spin-half central-spin Hamiltonian and all-down-reference Bethe equations, Eqs. (2)-(6), with no bath field: H=B\*S0^z+sum A\_j\*S0.Sj. Only fixed-magnetization ground energies for distinct nonzero couplings are implemented; no dynamics or form factors.

### sutherland-1968

Bill Sutherland. *Further Results for the Many-Body Problem in One Dimension*.
Phys. Rev. Lett. 20, 98-100 (1968).

[DOI](<https://doi.org/10.1103/PhysRevLett.20.98>).

Relevant tool modes:

- `bethe-sun-fermions-pbc`: Original multicomponent delta-gas solution. We implement equal-mass repulsive fermions on a periodic ring, not arbitrary mixtures or statistics.

### lee-2011

J. Y. Lee, X. W. Guan, and M. T. Batchelor. *Yang-Yang method for the thermodynamics of one-dimensional multi-component interacting fermions*.
J. Phys. A: Math. Theor. 44, 165002 (2011).

[arXiv](<https://arxiv.org/abs/1011.0128>).

Relevant tool modes:

- `bethe-sun-fermions-pbc`: Nested finite-size Bethe equations (2)-(3), with c'=c/2, and logarithmic signs from (13)-(19). Our Hamiltonian is -sum d\_j^2+2c sum delta with c\>=0 and no Zeeman term. Only the centered all-odd occupied-population ground branch is implemented; no strings, TBA or attraction.

### imambekov-demler-2006

Adilet Imambekov and Eugene Demler. *Exactly solvable case of a one-dimensional Bose-Fermi mixture*.
Phys. Rev. A 73, 021602(R) (2006).

[DOI](<https://doi.org/10.1103/PhysRevA.73.021602>), [arXiv](<https://arxiv.org/abs/cond-mat/0505632>).

### wang-1999

Yupeng Wang. *Exact solution of a spin-ladder model*.
Phys. Rev. B 60, 9236 (1999).

[DOI](<https://doi.org/10.1103/PhysRevB.60.9236>), [arXiv](<https://arxiv.org/abs/cond-mat/9901168>).

Relevant tool modes:

- `bethe-ladder-pbc`: Permutation-form ladder Hamiltonian (2), rung basis (3), chemical-potential form (4), and three nested rational equations (5). We use spin-1/2 operators, leg coefficient 1, four-spin coefficient 4, and J\_r=2\*J in (2)/(4), restoring E=E\_perm-L/4+J\_r\*(L/4-N\_s). Only zero-field periodic ground states and singlet-count sector minima are implemented.

### hakobyan-2004

Tigran Hakobyan. *The ordering of energy levels for SU(n) symmetric antiferromagnetic chains*.
Nucl. Phys. B 699, 575-594 (2004).

[DOI](<https://doi.org/10.1016/j.nuclphysb.2004.07.032>), [arXiv](<https://arxiv.org/abs/cond-mat/0403587>).

Relevant tool modes:

- `bethe-ladder-pbc`: SU(n) multiplets, Young diagrams and weight spaces. The periodic energy-ordering theorem has a row-parity restriction; we do not assume each population's own highest weight is lowest. Compatible dominant multiplets and displaced packed real seas are compared, not the full excited-state spectrum.

### bogoliubov-1997

N. M. Bogoliubov, A. G. Izergin, and N. A. Kitanine. *Correlation functions for a strongly correlated boson system*.
arXiv:solv-int/9710002 (1997).

[arXiv](<https://arxiv.org/abs/solv-int/9710002>).

Relevant tool modes:

- `bethe-q-boson-pbc`: Deformed Fock algebra (1.2)-(1.8), free/phase limits, and exact phase-model momenta (3.4). Our Hamiltonian is twice (1.1), with the +2N shift included.

### baxter-1973

R. J. Baxter. *Eight-vertex model in lattice statistics and one-dimensional anisotropic Heisenberg chain. III. Eigenvectors of the transfer matrix and Hamiltonian*.
Ann. Phys. 76, 48-71 (1973).

[DOI](<https://doi.org/10.1016/0003-4916(73)90441-7>).

### haldane-1988

F. D. M. Haldane. *Exact Jastrow-Gutzwiller resonating-valence-bond ground state of the spin-1/2 antiferromagnetic Heisenberg chain with 1/r^2 exchange*.
Phys. Rev. Lett. 60, 635 (1988).

[DOI](<https://doi.org/10.1103/PhysRevLett.60.635>).

Relevant tool modes:

- `bethe-haldane-shastry-pbc`: Original periodic spin-1/2 inverse-square chain; we use H=(pi/N)^2 sum\_{i\<j} S\_i.S\_j/sin^2(pi\*(i-j)/N), J=1.

### shastry-1988

B. Sriram Shastry. *Exact solution of an S=1/2 Heisenberg antiferromagnetic chain with long-ranged interactions*.
Phys. Rev. Lett. 60, 639 (1988).

[DOI](<https://doi.org/10.1103/PhysRevLett.60.639>).

Relevant tool modes:

- `bethe-haldane-shastry-pbc`: Independent exact solution of the long-range antiferromagnet; finite even and odd rings, not nearest-neighbor XXX.

### jiang-lamers-miao-2026

Yunfeng Jiang, Jules Lamers, and Yuan Miao. *Norms, overlaps and Yangian descendants for the Haldane-Shastry spin chain*.
arXiv:2606.20168v2 (2026).

[arXiv](<https://arxiv.org/abs/2606.20168v2>).

Relevant tool modes:

- `bethe-haldane-shastry-pbc`: Motif spectral rules and Yangian multiplicities, Secs. 2.2-2.3. Convert their H=sum(1-P)/(4\*sin^2) by H\_ours=E\_ferro-2\*(pi/N)^2\*H\_theirs. No wavefunctions, norms or overlaps are implemented.

### sutherland-1971

Bill Sutherland. *Exact Results for a Quantum Many-Body Problem in One Dimension*.
Phys. Rev. A 4, 2019 (1971).

[DOI](<https://doi.org/10.1103/PhysRevA.4.2019>).

Relevant tool modes:

- `bethe-sutherland-pbc`: Original trigonometric inverse-square gas. We select periodic scalar bosons with collision behavior |x\_i-x\_j|^lambda, lambda\>=0; the coupling coefficient alone does not fix this domain.

### gurappa-panigrahi-1999

N. Gurappa and Prasanta K. Panigrahi. *Equivalence of the Sutherland Model to Free Particles on a Circle*.
arXiv:hep-th/9908127v2 (1999).

[arXiv](<https://arxiv.org/abs/hep-th/9908127v2>).

Relevant tool modes:

- `bethe-sutherland-pbc`: Hamiltonian (18), ground energy (19), and the partition spectrum after (23), with beta=lambda and hbar=2m=1. Ascending integer labels include common boosts. Exact energies and momenta only; finite label windows do not claim global spectral completeness.

### andrei-1980

N. Andrei. *Diagonalization of the Kondo Hamiltonian*.
Phys. Rev. Lett. 45, 379 (1980).

[DOI](<https://doi.org/10.1103/PhysRevLett.45.379>).

### destri-de-vega-1992

C. Destri and H. J. de Vega. *New thermodynamic Bethe ansatz equations without strings*.
Phys. Rev. Lett. 69, 2313-2317 (1992).

[DOI](<https://doi.org/10.1103/PhysRevLett.69.2313>).

### gwa-spohn-1992

Leh-Hun Gwa and Herbert Spohn. *Bethe solution for the dynamical-scaling exponent of the noisy Burgers equation*.
Phys. Rev. A 46, 844 (1992).

[DOI](<https://doi.org/10.1103/PhysRevA.46.844>).

### yang-yang-1966

C. N. Yang and C. P. Yang. *One-Dimensional Chain of Anisotropic Spin-Spin Interactions. I. Proof of Bethe's Hypothesis for Ground State in a Finite System*.
Phys. Rev. 150, 321 (1966).

[DOI](<https://doi.org/10.1103/PhysRev.150.321>).

### shastry-sutherland-1990

B. Sriram Shastry and Bill Sutherland. *Twisted boundary conditions and effective mass in Heisenberg-Ising and Hubbard rings*.
Phys. Rev. Lett. 65, 243 (1990).

[DOI](<https://doi.org/10.1103/PhysRevLett.65.243>).

### de-vega-gonzalez-ruiz-1994

H. J. de Vega and A. González-Ruiz. *Boundary K-matrices for the XYZ, XXZ and XXX spin chains*.
J. Phys. A: Math. Gen. 27, 6129-6138 (1994).

[DOI](<https://doi.org/10.1088/0305-4470/27/18/021>), [arXiv](<https://arxiv.org/abs/hep-th/9306089>).

### essler-de-klerk-2023

F. H. L. Essler and A. J. J. M. de Klerk. *Statistics of matrix elements of local operators in integrable models*.
arXiv:2307.12410v1 (2023).

[arXiv](<https://arxiv.org/abs/2307.12410v1>), [Equations (HTML)](<https://arxiv.org/html/2307.12410v1>).

Relevant tool modes:

- `bethe-lieb-liniger-pbc`: Explicit normalization and finite-ring equations, Eqs. (4), (24)-(27); root-density equation (32)-(33) used for validation, not a thermodynamics API. No matrix elements are implemented.

### yang-yang-1969

C. N. Yang and C. P. Yang. *Thermodynamics of a One-Dimensional System of Bosons with Repulsive Delta-Function Interaction*.
J. Math. Phys. 10, 1115 (1969).

[DOI](<https://doi.org/10.1063/1.1664947>).

Relevant tool modes:

- `bethe-lieb-liniger-thermal`: Finite-temperature thermodynamics of repulsive continuum bosons, at fixed chemical potential or by inversion at fixed density; hbar=2m=k\_B=1.

### doikou-nepomechie-1998

Anastasia Doikou and Rafael I. Nepomechie. *Bulk and Boundary S Matrices for the SU(N) Chain*.
Nucl. Phys. B 521, 547-572 (1998).

[DOI](<https://doi.org/10.1016/S0550-3213(98)00239-9>), [arXiv](<https://arxiv.org/abs/hep-th/9803118>).

Relevant tool modes:

- `bethe-su3-pbc`: Nested equations and energy, Eqs. (2.17)-(2.19), logarithmic labels (2.24)-(2.29), and the filled-sea singlet in Sec. 2.3. Our H=sum P gives E=2\*E\_paper+L. No strings, S matrices or boundary fields are implemented.

### oelkers-2006

N. Oelkers, M. T. Batchelor, M. Bortz, and X.-W. Guan. *Bethe Ansatz study of one-dimensional Bose and Fermi gases with periodic and hard wall boundary conditions*.
J. Phys. A 39, 1073-1098 (2006).

[DOI](<https://doi.org/10.1088/0305-4470/39/5/005>), [arXiv v2](<https://arxiv.org/abs/cond-mat/0511694v2>).

Relevant tool modes:

- `bethe-gaudin-yang-pbc`: Hamiltonian and periodic fermion equations (1), (2), (7), (8), (25); even N, odd minority population as in Sec. 5. Weak and strong limits (12), (15)-(16) provide checks. Hard walls, attraction and excitations are not implemented.

### grijalva-2019

Sebastian Grijalva, Jacopo De Nardis, and Veronique Terras. *Open XXZ chain and boundary modes at zero temperature*.
SciPost Phys. 7, 023 (2019).

[DOI](<https://doi.org/10.21468/SciPostPhys.7.2.023>), [arXiv v4](<https://arxiv.org/abs/1901.10932v4>).

Relevant tool modes:

- `bethe-xxz-obc`: Massive open XXZ boundary root and finite-size deviation, Sec. 4.3.2, specialized to zero boundary fields; used for Delta\>1 ground states, not boundary correlations or excitation scans.

### dugave-2015

M. Dugave, F. Göhmann, K. K. Kozlowski, and J. Suzuki. *On form-factor expansions for the XXZ chain in the massive regime*.
J. Stat. Mech. (2015) P05037 (2015).

[DOI](<https://doi.org/10.1088/1742-5468/2015/05/P05037>), [arXiv v2](<https://arxiv.org/abs/1412.8217v2>).

Relevant tool modes:

- `bethe-xxz-pbc`: Massive periodic XXZ Hamiltonian and trigonometric Bethe equations, Eqs. (1.1)-(1.2), with zero twist and field; Pauli exchange is divided by four. Used for Delta\>1 sector ground states, not form factors, strings or excitations.

### caux-xxz-coordinate

Jean-Sébastien Caux. *The Bethe Ansatz: coordinate wavefunctions and periodic XXZ equations*.
Online notes.

[General M](<https://integrability.org/c_h_s_m.html>), [Two-magnon scattering](<https://integrability.org/c_h_s_2.html>).

### korepin-2009

Vladimir E. Korepin. *Norm of Bethe Wave Function as a Determinant*.
arXiv:0911.1881 (historical account of the 1982 norm formula) (2009).

[arXiv](<https://arxiv.org/abs/0911.1881>), [HTML](<https://arxiv.org/html/0911.1881v1>), [Original 1982 paper](<https://doi.org/10.1007/BF01212176>).

### bini-robol-2013

Dario A. Bini and Leonardo Robol. *Solving secular and polynomial equations: a multiprecision algorithm*.
Author manuscript, May 10, 2013 (2013).

[Author manuscript](<https://web.dm.unipi.it/robol/assets/pdf/secular-paper.pdf>).

### cameron-graillat-2022

Thomas R. Cameron and Stef Graillat. *On a compensated Ehrlich-Aberth method for the accurate computation of all polynomial roots*.
Electronic Transactions on Numerical Analysis 55, 401-423 (2022).

[DOI](<https://doi.org/10.1553/etna_vol55s401>), [Open-access article](<https://etna.ricam.oeaw.ac.at/vol.55.2022/pp401-423.dir/pp401-423.pdf>).

### koma-nachtergaele-1997

Tohru Koma and Bruno Nachtergaele. *The spectral gap of the ferromagnetic XXZ chain*.
Lett. Math. Phys. 40, 1-16 (1997).

[arXiv](<https://arxiv.org/abs/cond-mat/9512120>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Ferromagnetic OBC gap and one-defect band: Proposition 2 and Eq. (3.30), transferred through TL equivalence at Delta=3/2 and multiplied by 2\*Delta=3. Gap above the entire ground space is 3-2\*cos(pi/N), not a zero-mode splitting.

### zhou-2025-biquadratic

Huan-Qiang Zhou, Qian-Qian Shi, Ian P. McCulloch, and Murray T. Batchelor. *Goldstone modes and the golden spiral in the ferromagnetic spin-1 biquadratic model*.
J. Phys. A: Math. Theor. 58, 39LT01 (2025).

[DOI](<https://doi.org/10.1088/1751-8121/ae05d9>), [arXiv](<https://arxiv.org/abs/2302.13126>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Ferromagnetic convention H=+sum (S.S)^2, Eq. (3), and Fibonacci ground-space multiplicities under free ends. Ground-space zero modes are distinct from the positive-energy TL defect band; entanglement and ground-state wavefunctions are not implemented.

### nachtergaele-spitzer-starr-2007

Bruno Nachtergaele, Wolfgang Spitzer, and Shannon Starr. *Droplet Excitations for the Spin-1/2 XXZ Chain with Kink Boundary Conditions*.
Ann. Henri Poincare 8, 165-201 (2007).

[DOI](<https://doi.org/10.1007/s00023-006-0304-6>), [arXiv](<https://arxiv.org/abs/math-ph/0508049>).

Relevant tool modes:

- `bethe-biquadratic-obc`: Droplet interpretation and thermodynamic module-edge limit, Theorem 2.1; rescaling by 2\*Delta=3 gives thresholds 5/3 for two defects and 2 for three. Finite-chain pair/triple modes solve the original Bethe equations with finite deviations, not a bulk/ideal-string substitution. No form factors or physical-spin decomposition are implemented.

### reichert-2019

Benjamin Reichert, Grigori E. Astrakharchik, Aleksandra Petkovic, and Zoran Ristivojevic. *Exact Results for the Boundary Energy of One-Dimensional Bosons*.
Phys. Rev. Lett. 123, 250602 (2019).

[DOI](<https://doi.org/10.1103/PhysRevLett.123.250602>), [arXiv](<https://arxiv.org/abs/1908.08172>).

Relevant tool modes:

- `bethe-lieb-liniger-obc`: Equations (7)-(8) fix the hard-wall ground-state convention and exclusion of self-image scattering; the finite-volume solver extends the logarithmic labels to selected excitations. No boundary-energy integral-equation or thermodynamics API is implemented.

### franchini-2011

Fabio Franchini. *Notes on Bethe Ansatz Techniques*.
SISSA lecture notes, May 15, 2011 (2011).

[Open lecture notes](<https://people.sissa.it/~ffranchi/BAnotes.pdf>).

Relevant tool modes:

- `bethe-lieb-liniger-thermal`: Section 2.11: thermal filling and entropy of Lieb-Liniger Bethe states. No attractive strings, trapped gases or dynamical correlations.
- `bethe-lieb-liniger-dispersion`: Dressed-energy equation and Fermi boundary condition (2.87)-(2.88), related particle-hole energies (2.95); our phase is +2 atan(k/c).

### caux-lieb-liniger

Jean-Sébastien Caux. *The Bethe Ansatz: The ground state and the Lieb equation*.
Online lecture notes, integrability.org.

[Lieb equation](<https://integrability.org/g_l_Le.html>).

Relevant tool modes:

- `bethe-lieb-liniger-dispersion`: Pedagogical ground-state integral equation and normalization at fixed density.

### pozsgay-2014-q-boson

B. Pozsgay. *Quantum quenches and Generalized Gibbs Ensemble in a Bethe Ansatz solvable lattice model of interacting bosons*.
J. Stat. Mech. (2014) P10045 (2014).

[DOI](<https://doi.org/10.1088/1742-5468/2014/10/P10045>), [arXiv](<https://arxiv.org/abs/1407.8344>).

Relevant tool modes:

- `bethe-q-boson-pbc`: Periodic q-boson Hamiltonian (2.1), Bethe equations (2.13), and shifted energy 4 sum sin^2(k/2). Fixed-N ground states and canonical real-root excitation scans; no quenches or correlation functions.

### guan-batchelor-2011

X.-W. Guan and M. T. Batchelor. *Polylogs, thermodynamics and scaling functions of one-dimensional quantum many-body systems*.
J. Phys. A: Math. Theor. 44, 102001 (2011).

[DOI](<https://doi.org/10.1088/1751-8113/44/10/102001>), [arXiv](<https://arxiv.org/abs/1010.4842>).

Relevant tool modes:

- `bethe-lieb-liniger-thermal`: Yang-Yang pseudoenergy, root density and pressure equations (2)-(5). We solve the full integral equations, not the finite-c polylog expansion.

<!-- END GENERATED BIBLIOGRAPHY -->
