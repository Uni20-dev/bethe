# References and provenance

The Heisenberg ground-state solver succeeds Ian McCulloch's
`misc/heisenberg-energy.cpp` in the Matrix Product Toolkit (copyright
2015-2023). The original is maintained there, rather than duplicated here.
The successor retains its GPL-3.0-or-later notice.

The logarithmic equations, quantum numbers, energy convention, and real-root
iteration follow Michael Karbach, Kun Hu, and Gerhard Müller,
*Introduction to the Bethe ansatz II*, Computers in Physics **12**, 565 (1998),
[doi:10.1063/1.168740](https://doi.org/10.1063/1.168740),
[arXiv:cond-mat/9809163](https://arxiv.org/abs/cond-mat/9809163), Eqs. (6)-(9).
Table I supplies the N=16 regression reference, including its N/4 energy shift.

The sector-state selection extends the consecutive quantum-number construction
in that paper (Eq. (16)). For odd-chain one-spinon states and hole labeling, see
Stefan Groha and Fabian H. L. Essler, *Spinon decay in the spin-1/2 Heisenberg
chain with weak next nearest neighbour exchange*, J. Phys. A **50**, 334002
(2017), [arXiv:1702.06550](https://arxiv.org/abs/1702.06550), Sec. 3. The present
solver concerns the unperturbed integrable chain, not the decay calculation.

The zero-field thermodynamic XXX and gapless XXZ dispersions follow
Jean-Sébastien Caux's *The Bethe Ansatz* notes:
[XXX spinons](https://integrability.org/g_h_e.html) and
[XXZ spinons](https://integrability.org/g_sc_p_e.html). We use positive spinon
momentum k in [0,pi], reversing the sign of the notes' convention. At finite
odd N our explicit convention is k=pi/2-2*pi*I_h/N, related to lattice momentum
by P=pi*M+pi/2-k modulo 2*pi. The bulk reference energy is J*(1/4-log(2)).

Research using these tools should acknowledge the Matrix Product Toolkit and
Uni20 where relevant, and cite the methods used in the calculation.

The free-end XXX equations follow the rational limit of Eqs. (11)-(12) and
footnote 2 in Zhongtao Mei and C. J. Bolech, *Derivation of matrix product states
for the Heisenberg spin chain with open boundary conditions*, Phys. Rev. E
**95**, 032127 (2017), [arXiv:1609.08045](https://arxiv.org/abs/1609.08045).
Our rapidity is `z=2*lambda`; our spin-1/2 Hamiltonian is the paper's XXX
Hamiltonian divided by four and shifted by `N/4`. This gives the ferromagnetic
reference `(N-1)/4`. Taking the positive-root logarithmic branch gives the
integer quantum-number convention documented in the
[open-chain guide](docs/open-chains.md#allowed-states-and-convergence).

For the highest-weight multiplet and string classification relevant to
real-root excitation scans, see Caux's notes on
[real rapidities](https://integrability.org/c_h_e_rr.html),
[SU(2) descendants](https://integrability.org/c_h_e_rr_10.html), and
[two-string states](https://integrability.org/c_h_e_s_2.html).
The scans enumerate only the supported all-real quantum-number window;
they do not implement the string sectors described in those references.

The finite-size periodic XXZ sector solver follows the logarithmic equations,
energy and momentum in R. Vlijm, I. S. Eliëns, and J.-S. Caux,
*Correlations of zero-entropy critical states in the XXZ model: integrability
and Luttinger theory far from the ground state*, SciPost Phys. **1**, 008 (2016),
[doi:10.21468/SciPostPhys.1.1.008](https://doi.org/10.21468/SciPostPhys.1.1.008),
[arXiv:1606.09516](https://arxiv.org/abs/1606.09516), Eqs. (1), (3)-(5).
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

The free-end XXZ solver uses Mei and Bolech's Eqs. (11)-(12), cited above,
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
