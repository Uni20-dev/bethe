# Spin-1 Takhtajan–Babujian chains

[Overview](../README.md) | [Model catalogue](models.md) | [Precision and CLI controls](command-line.md)

`bethe-tb-pbc` calculates the zero-field singlet ground state of the periodic
spin-1 Takhtajan–Babujian (TB, also Babujan–Takhtajan) chain for **even L>=4**:

```math
H=\sum_j\left[\mathbf S_j\cdot\mathbf S_{j+1}-(\mathbf S_j\cdot\mathbf S_{j+1})^2\right].
```

The S operators are spin 1. The bilinear coefficient is 1, with no additive
constant. A polarized bond has energy 0; its total-spin-0, 1, 2 eigenvalues
are -6, -2, 0. This normalization is **J=4** in
[Vlijm–Caux](../CITATIONS.md#vlijm-caux-2014), Eq. (1.2), whose Hamiltonian
has a J/4 prefactor. It differs from the often-used coefficients
`cos(theta), sin(theta)` at theta=-pi/4 by a factor sqrt(2).

This is not the generic bilinear spin-1 Heisenberg chain, and not the
[SU(3) ULS point](su3.md), whose biquadratic term has the opposite sign.
Odd lengths, other spin representations, magnetization sectors, excitations,
fields and open boundaries are not yet supported here.

## First calculations

```sh
build/bethe-tb-pbc 16
build/bethe-tb-pbc 6 --roots --precision long-double
# With binary128 enabled:
build/bethe-tb-pbc 64 --precision fp128
```

The four-site ground energy is $-11-\sqrt{41}$, approximately
`-17.40312423743284868648821767`. Ground momentum is zero for this even-ring
family. The report includes both phase and modulus residuals, CPU time, and
accepted Newton updates. `--roots` prints string centres and deviations,
then the real and imaginary parts of every rapidity.

## Why complex roots are necessary

Relative to the all-up spin-1 reference, the singlet has L spin lowerings,
and hence L rapidities. They form L/2 conjugate pairs:

```math
\lambda_{j,\pm}=x_j\pm i(1/2+\delta_j).
```

An ideal two-string sets delta_j=0. That is not the finite-chain solution:
the original equations then contain a singular within-pair scattering
factor. Here both x_j and delta_j are solved, retaining positive finite
deviations. Even the four-site energy changes appreciably if the imaginary
parts are replaced by +/-1/2 after solving for the centres.

The original complex equations and energy, in our normalization, are

```math
\left(\frac{\lambda_j+i}{\lambda_j-i}\right)^L
=\prod_{k\ne j}\frac{\lambda_j-\lambda_k+i}{\lambda_j-\lambda_k-i},
\qquad E=-4\sum_j\frac1{1+\lambda_j^2}.
```

Conjugate contributions make E real. Translation obeys
`exp(iP)=product_j (lambda_j+i)/(lambda_j-i)`; reflection and conjugation
give P=0 for this filled even-ring sea. Tests verify that product independently
as well as translation in the spin basis.

## Finite-deviation equations and branch selection

The filled two-string sea has consecutive centred string labels
$I_{j} =j -(L /2-1)/2$, j=0,...,L/2-1. These are not the quantum numbers of the
individual complex roots. In the relation (3.9) of Vlijm–Caux, the ordered
string-sign sum cancels I_j, giving $J _++J _-=0$. The half-odd-integer root
labels and Eq. (3.11) select positive deviations, $J _+=-1/2, J _-=1/2$.
Coinciding root quantum numbers across different strings do not imply
coinciding rapidities.

For clarity, the actual equations are recorded here. Define
$A (a,b)=\operatorname{atan2} (a,b)$, $B (a,b)=\log (a \,a +b \,b)/2$, $d_{\mathrm{jk}} =x_{j} -x_{k}$, and widths

```math
w_{jk}=(2+\delta_j+\delta_k,-\delta_j-\delta_k,1+\delta_j-\delta_k,1-\delta_j+\delta_k),
\qquad \mathrm{signs}=(+1,-1,+1,-1).
```

Our normalized residuals are

```math
\begin{aligned}
F_j&=A(x_j,3/2+\delta_j)+A(x_j,1/2-\delta_j)
-\frac1L\sum_{k\ne j}\sum_{a=1}^{4}A(d_{jk},w_{jk}[a]),\\
G_j&=B(x_j,3/2+\delta_j)-B(x_j,1/2-\delta_j)
-\frac{\log(1+\delta_j)-\log\delta_j}{L}\\
&\quad-\frac1L\sum_{k\ne j}\sum_{a=1}^{4}\mathrm{signs}[a]B(d_{jk},w_{jk}[a]).
\end{aligned}
```

These specialize Eqs. (3.8) and (3.10) to the ground-state two-string sea.
The atan2 branch matters when the second width is negative. Its first
argument is nonzero for distinct centres, so no branch value at (0,0) is
needed. The self-pair modulus is evaluated from delta directly, avoiding
subtraction of two nearly equal imaginary parts. Positive deviations below
1/2 keep this branch away from the driving poles at lambda=+/-i.

Reflection symmetry removes redundant equations: centres occur as +/-x,
deviations are equal in each reflected pair, and a central string at x=0
has an identically zero phase equation. There are L/2 real Newton variables.
The seed uses the thermodynamic density $1/(2\,\cosh (\pi \,x))$ and the asymptotic
deviation estimate in Eq. (3.18). **Only the seed is approximate:** convergence
is assessed using the finite-size equations above.

## Library and numerical contract

```cpp
#include <bethe/takhtajan_babujian.hpp>
auto state = bethe::takhtajan_babujian::ground_state<long double>(16);
if (!state.converged) { /* state.energy is an incomplete estimate */ }
```

The result retains centres, deviations, exact half-integer string labels,
adjacent conjugate rapidities, energy, momentum, residuals, status and update
count. Computation and reporting retain fp64, long double or optional fp128.
There is no precision fallback or tolerance relaxation.

Newton uses an analytic real Jacobian, recoverable partial-pivot linear
solves, and a residual-decreasing line search preserving ordered centres and
positive, representable deviations. Dense storage is O(L^2), factorization
O(L^3). Allocation-overflow checks are not a promise that very large chains
fit in available memory. Singular/tiny pivots return an unresolved status
without changing Uni20's error policy.

`residual_norm=max(max|F_j|,max|G_j|)` defaults to a tolerance of 32 epsilon.
It is an equation-residual target, **not an energy-error bound**. A zero
iteration budget evaluates the seed only. Incomplete solves retain the
final iterate and report iteration limit, stalled, or ill-conditioned.
The CLI exits 0 on convergence, 2 for an incomplete solve, 1 for invalid or
unsupported input. Convergence identifies the prescribed ground-state
branch; it is not a scan or certificate of the entire spectrum.

Validation includes the original multiplicative complex equations through
L=128 in every scalar type; independent spin-1 matrix construction and
ground-state energy/translation checks at L=4,6,8; the native-precision
four-site analytic energy; finite-difference Jacobian checks; symmetry,
budget, precision-I/O and invalid-input tests. The spin-basis ED oracle
deliberately uses fp64, separately from native-precision analytic/root checks.

The next extension is state selection with broken two-strings. Real roots,
three-strings, singular solutions and their quantum-number branches need
explicit treatment; changing I_j in this ground-only solver is insufficient.
