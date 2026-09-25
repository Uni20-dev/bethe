# Periodic SU(3) permutation chain

[Overview](../README.md) · [CLI controls](command-line.md) · [Model catalogue](models.md#su-n-permutation-chains-and-the-spin-1-uls-point)

Each site has three states, called colors. The nearest-neighbor operator P
exchanges them: $`P \lvert a,b \gt  = \rvert b,a \gt `$. We solve the antiferromagnetic ring

```math
H=\sum_{j=1}^{L}P_{j,j+1},\qquad P_{L,L+1}=P_{L,1},\qquad J=1.
```

This is a useful non-Abelian lattice benchmark without a hopping or coupling
parameter to tune. It is also the spin-1 Uimin-Lai-Sutherland (ULS) point,
with an energy shift explained below, **not** generic spin-1 bilinear exchange.

## First calculation and supported scope

```sh
build/bethe-su3-pbc 6
build/bethe-su3-pbc 9 --roots --precision long-double
build/bethe-su3-pbc 48 --precision fp128 --format plain
```

The last command requires an [MPLAPACK-enabled build](building.md#enable-binary128).
fp64 is the default; long-double and fp128 retain their native arithmetic,
tolerance parsing, and output precision.

Currently supported: the balanced singlet ground state for **L>=3 divisible
by 3**, with L/3 sites of each color, ordinary periodic boundaries, and no
fields. Lengths 3, 6, 9, ... include both parities. There are two coupled
families of finite real roots, and no excited-state enumeration in this
frontend. Other lengths, unequal populations, descendants, complex strings,
twists, open ends, and general SU(n) are separate extensions.

The six-site result is $`E =-1-\sqrt{13}`$, approximately -4.605551275463989.
The report includes populations, both root counts, total and per-site energy,
momentum, residuals for each equation family, Newton updates, and solver CPU
time. `--roots` prints two tables with exact integer/half-integer labels.
Run `--help` for controls and `--references` for relevant literature; normal calculation
output does not repeat the bibliography.

## Why two sets of roots?

Starting from a reference color, the first nesting level describes sites
which are not that color. The second level distinguishes the remaining two
colors. With root counts M1 and M2, the populations are
$`(L -M_1, M_1 -M_2, M_2)`$. The balanced singlet therefore has

```math
\begin{aligned}
M_1&=\frac{2L}{3},\qquad M_2=\frac L3,\\{}
I_j&=j-\frac{M_1-1}{2},\quad j=0,\ldots,M_1-1,\\{}
J_a&=a-\frac{M_2-1}{2},\quad a=0,\ldots,M_2-1.
\end{aligned}
```

These are consecutive, centered labels with no holes. The first-level labels
are always half-odd integers; the second-level labels are integers when M2
is odd, and half-odd integers when M2 is even. At L=3, for example,
$`I ={-1/2,1/2}`$, $`J ={0}`$. At L=6 they are
$`I ={-3/2,-1/2,1/2,3/2}`$, $`J ={-1/2,1/2}`$.

The state selection follows the filled-sea singlet in Sec. 2.3 of
[Doikou and Nepomechie](../CITATIONS.md#doikou-nepomechie-1998), building on
[Sutherland's multicomponent solution](../CITATIONS.md#sutherland-1975).
Selecting this state uses more information than just solving equations to a
small residual; arbitrary centered counts are not a public ground-state API.

## Equations and normalization

Use conventional rapidities lambda and mu and $`\theta_{n} (x)=2 \arctan (2x /n)`$:

```math
\begin{aligned}
L\theta_1(\lambda_j)-\sum_{k\ne j}\theta_2(\lambda_j-\lambda_k)
+\sum_a\theta_1(\lambda_j-\mu_a)&=2\pi I_j,\\{}
\sum_j\theta_1(\mu_a-\lambda_j)-\sum_{b\ne a}\theta_2(\mu_a-\mu_b)&=2\pi J_a,\\{}
E&=L-\sum_j\frac1{\lambda_j^2+1/4},\\{}
P&=\sum_j[\pi-\theta_1(\lambda_j)]\pmod{2\pi}.
\end{aligned}
```

Only the first-level roots enter energy and momentum directly, but they are
determined by **both** equations. The mu roots are not physical momenta.
Our roots are not the XXX tool's $`z =2\,\lambda`$.

These are the all-real specialization of the reference's Eqs. (2.17)-(2.19)
and (2.24)-(2.29), using L for length and 3 for the number of colors.
Its Hamiltonian is `sum(P-1)/2`, so **$`E_{\mathrm{here}} =2\,E_{\mathrm{paper}} +L`$**.
There are L bonds, including the closing bond.

Both root sets are reflection symmetric. A lambda pair contributes 2*pi
to the momentum sum; M1 is even, so the supported state has P=0 exactly.
We report `momentum_index=0`, with the lattice convention
`P=2*pi*momentum_index/L (mod 2*pi)`.

## Numerical method and failure reporting

The implementation solves only the positive roots. Negative roots are their
exact reflections, and an odd-sized sea has a fixed zero root. The Newton
matrix has order `M1/2 + floor(M2/2)`, about L/2, rather than L.
The full arrays, including negative and zero roots, are returned to callers.

The starting guess inverts the integrated thermodynamic filled-sea density.
For nesting level a=1,2 and a positive label q, it is

```math
x_a(q)=\frac3\pi\operatorname{artanh}\!\left[\tan(\pi q/L)\tan(\pi a/6)\right].
```

This is only a seed: the solver refines the **finite-size equations**, not a
bulk approximation. We use an analytic reflection-reduced Jacobian, Uni20's
dense linear solve, and a damped Newton step preserving positive, strictly
ordered roots within each sea. Phase sums and energy use compensated summation.
Dense factorization costs O(L^3) per update and O(L^2) storage.

Each residual is its logarithmic equation mismatch divided by L. The report
takes the maximum absolute residual of the independent positive-root
equations in each family, then the maximum of those two values. Reflection
supplies the negative equations; the central-zero equations vanish by
symmetry. Tests separately evaluate **all** equations from the returned
roots. A family with no positive roots has reported residual zero.

The default tolerance is 32 times the selected scalar type's epsilon.
It is an equation criterion, **not an energy-error bound**. `--max-iterations`
counts accepted Newton updates (default 10000); zero returns the evaluated
density seed. Unattainable tolerances can stall at representable precision.

Exit status is 0 for convergence, 2 for a budget-limited or stalled solve, and
1 for invalid input or another error. Nonconverged energies and roots remain
explicitly labeled estimates, with diagnostics for the actual returned roots.
There is no precision fallback. Increasing L can amplify root conditioning
and energy cancellation: check convergence and compare precision/tolerance
before interpreting the last printed digits.

## C++ library

```cpp
#include <bethe/su3.hpp>

auto state = bethe::su3::ground_state<long double>(12);
if (!state.converged) {
  // Inspect state.status, state.residual_norm, and state.iterations.
}
auto const& lambda = state.rapidities[0];
auto const& mu = state.rapidities[1];
auto const& I = state.quantum_numbers[0];
auto const& J = state.quantum_numbers[1];
```

`State<Real>` keeps `populations`, both root arrays, both label arrays,
`energy`, `momentum`, `momentum_index`, and convergence diagnostics.
`ground_quantum_numbers(L)` exposes the label construction. Labels use
`uni20::half_int`; all real-valued state data use `Real`. As in the other
modules, `SolverOptions<Real>` supplies `residual_tolerance` and
`max_iterations`. The internal equation class is not a public arbitrary-sector
solver, and the Hubbard API is unchanged.

## Independent checks and the spin-1 mapping

For L=3, the totally antisymmetric color singlet has energy -3. Its roots are
$`\lambda ={-1/\sqrt{12},1/\sqrt{12}}`$, $`\mu ={0}`$. This checks an irrational root
in every precision, not just an exactly representable energy.

For L=6, form singlets by multiplying two three-color Levi-Civita tensors:
`B_A(c_1,...,c_6)=epsilon(c_A)*epsilon(c_complement)`, with each site's indices
in increasing order. Choose A=(123),(124),(125),(134),(135). These five
independent states span the singlet space. Acting with the six swaps gives
$`H B = B R`$, where our exact color-space calculation gives

```math
R=\begin{pmatrix}
-5&1&1&1&-3\\{}
2&-2&0&0&2\\{}
-1&1&-1&1&-1\\{}
-1&1&1&-1&-1\\{}
1&1&1&1&3
\end{pmatrix}.
```

R is not symmetric because this basis is not orthonormal. Its characteristic
polynomial is $`E \,(E +2)^{2}\,(E ^{2}+2\times 10^{-12})`$, giving the singlet minimum
$`-1-\sqrt{13}`$. Independent full-Hilbert-space diagonalization confirms it is
the six-site ground energy. The analytic value is evaluated in each test's
scalar type; fp128 and extended long-double checks reject a double-rounded
answer. The ED oracle itself is intentionally fp64 and is not claimed as a
high-precision reference.

The tests also compare L=3,6,9 against direct color-space diagonalization
and its translation eigenvalue, verify the original multiplicative equations
(including self factors and logarithmic parity), check the Jacobian by finite differences,
and delete the auxiliary sea in the same equation implementation to recover
the SU(2) chain: $`E_{\mathrm{permutation}} =2\,E_{\mathrm{XXX}} +L /2`$, $`z_{\mathrm{XXX}} =2\,\lambda`$.
Finite sizes through L=192 approach the bulk reference
$`e_{\mathrm{infinity}} =1-\log (3)-\pi /(3\,\sqrt{3})`$, obtained by rescaling the reference's
Eq. (2.49). This bulk value is a validation target, not a thermodynamics API.

Finally, for two spin-1 sites the swap identity is

```math
P=\mathbf S_i\cdot\mathbf S_j+(\mathbf S_i\cdot\mathbf S_j)^2-1.
```

It is verified directly as a 9-by-9 matrix identity. Thus for
`H_ULS=sum[S_i.S_j+(S_i.S_j)^2]`, the energy is **$`E_{\mathrm{ULS}} =E +L`$**.
If instead the bilinear-biquadratic Hamiltonian uses coefficients
`cos(theta), sin(theta)` at theta=pi/4, its energy is $`(E +L)/\sqrt{2}`$.
Neither convention changes the roots. Do not apply this mapping to another
angle or assume that SU(3) singlet coverage includes all spin-1 sectors.
