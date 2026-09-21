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
integer quantum-number convention documented in README.md.
