# Supersymmetric t–J chains

[Overview](../README.md) | [Model catalogue](models.md) | [Precision and CLI controls](command-line.md)

`bethe-tj-pbc` solves selected fixed-population ground states of the
periodic t–J chain at **t=1, J=2**. Each site is empty, spin up, or spin down:
double occupancy is excluded. This is not a solver at arbitrary J/t.

```math
H=-\sum_{j,\sigma}\mathcal P
\left(c^\dagger_{j,\sigma}c_{j+1,\sigma}+\mathrm{h.c.}\right)\mathcal P
+2\sum_j\left[\mathbf S_j\cdot\mathbf S_{j+1}-\frac{n_jn_{j+1}}4\right].
```

P projects out double occupancy. Fermions have periodic boundary conditions,
including the boundary hopping sign in a Fock basis. We require L>=3, so
there is no two-site bond-counting ambiguity. There is no magnetic field or
chemical-potential term. Energies include the displayed density interaction.

## First calculations and supported populations

```sh
# No holes: the Heisenberg limit, with its exchange factor and energy shift.
build/bethe-tj-pbc 6
# Two holes, three electrons of each spin.
build/bethe-tj-pbc 8 --particles 6 --roots
# An imbalanced doped sector: N_up=3, N_down=1.
build/bethe-tj-pbc 5 --particles 4 --sz 1 --precision long-double
# Exact polarized free fermions.
build/bethe-tj-pbc 7 --particles 4 --sz 2 --roots
```

`--particles N` defaults to L. `--sz` defaults to 0 for even N or 1/2 for
odd N, and fixes $`N_{\mathrm{up}} =N /2+S^z`$, $`N_{\mathrm{down}} =N /2-S^z`$. Every returned energy is
a minimum **within those populations**, not a minimization over filling.
Three branches are supported:

- **Doped, mixed spin:** both N_up and N_down must be odd. Either parity
  of L is allowed. Spin reversal is supported and reuses the same equations.
- **No holes:** every physical spin population, using the periodic XXX
  sector solver. This includes odd L and populations excluded by the doped
  branch's parity restriction.
- **Vacuum or full polarization:** any particle count from 0 to L, using
  exact periodic free-fermion occupations.

Thus `8 --particles 4` is currently rejected: its default N_up=N_down=2
lies outside the doped real-root family. Do not silently shift its quantum
numbers or spin to a supported state. Excitations, other doped shell
parities, strings, open boundaries, and J!=2t are not implemented.

## Why the Sutherland formulation?

[Essler–Korepin](../CITATIONS.md#essler-korepin-1992) present three nested
gradings. The Sutherland/BFF form, their Eq. (3.73), gives two real-root
families for the branch used here. With spin reversal folded so N_up>=N_down,
their counts are $`M_1 =N_{h} +N_{\mathrm{down}}`$, $`M_2 =N_{h}`$, not N_e and N_down. At no holes
the second family disappears and the first equations reduce to XXX. The
Lai/FFB and alternative FBF gradings are not additional implemented solvers.

Define $`e_{a} (x)=(x +i \,a /2)/(x -i \,a /2)`$. The original equations are

```math
\begin{aligned}
e_1(\lambda_j)^L&=\frac{\prod_{k\ne j}e_2(\lambda_j-\lambda_k)}
{\prod_\alpha e_1(\lambda_j-\mu_\alpha)},\\{}
1&=\prod_j e_1(\mu_\alpha-\lambda_j).
\end{aligned}
```

There is **no mu–mu self-scattering**. It would be incorrect to copy the
SU(3) equations unchanged. In logarithmic form, with $`\theta_{a} (x)=2\,\arctan (2\,x /a)`$,

```math
\begin{aligned}
L\theta_1(\lambda_j)-\sum_{k\ne j}\theta_2(\lambda_j-\lambda_k)
+\sum_\alpha\theta_1(\lambda_j-\mu_\alpha)&=2\pi I_j,\\{}
\sum_j\theta_1(\mu_\alpha-\lambda_j)&=2\pi J_\alpha.
\end{aligned}
```

We use consecutive centered labels at each level:
$`I_{j} =j -(M_1 -1)/2`$, $`J_{\alpha} =\alpha -(M_2 -1)/2`$, with zero-based indices.
Taking the signs of the original rational equations into account requires
`2I = L-M1+1-M2 (mod 2)` and `2J = M1 (mod 2)`. For these centered seas,
these conditions require both spin populations to be odd. This is a branch
restriction, not a prohibition on other physical populations. Both the
logarithmic and original multiplicative equations are tested independently.

## Energy and momentum conventions

The paper's supersymmetric Hamiltonian includes a chemical-potential shift:
$`H_{\mathrm{susy}} =H +2\,N_{e} -L`$, Eq. (1.5). Removing it from Eq. (3.75) gives our energy

```math
E=2N_h-\sum_j\frac1{\lambda_j^2+1/4}.
```

The primary rapidities are conventional lambda, not the XXX executable's
z=2*lambda. For no holes,

```math
H_{tJ}=2H_{\mathrm{XXX}}-\frac L2.
```

Numerically we sum the root contributions directly, rather than subtracting
these two extensive terms. This preserves an O(1) nearly polarized energy
even when L is large enough for that subtraction to lose it in fp64.

This energy reduction does not permit copying spin-chain momentum blindly.
Translation of the fully occupied fermionic reference has sign
$`(-1)^{L -1}`$. Add pi to the XXX spin momentum on even L, and zero on odd L.
The centered doped family has total fermionic momentum zero. These conventions
are checked against translation in the projected Fock-space ground eigenspace.

Fully polarized states fill integer modes closest to zero, $`k =2\,\pi \,j /L`$.
For even particle count we choose the sea with positive total momentum;
its reflected partner is degenerate (and coincides modulo L at full filling).
The root report prints those free modes instead of manufacturing nested roots.

## Library and numerical contract

```cpp
#include <bethe/tj.hpp>
auto state = bethe::tj::ground_state<long double>(8, 3, 3);
// Arguments: sites, N_up, N_down, optional SolverOptions<Real>.
if (!state.converged) { /* report an incomplete solve */ }
```

All finite calculations retain fp64, long double, or optional fp128. The
doped solver works with the positive halves of reflection-paired roots,
inserting exact zero roots for odd counts. Newton uses an analytic Jacobian
and a line search preserving positive, ordered roots. A small native-precision
pivoted linear solve returns an unresolved status for singular/numerically
tiny pivots; it does not invoke or modify a process-global abort policy.
Dense Newton storage is O((M1+M2)^2); factorization is cubic. The tangent seed
is only an initial guess, not a thermodynamic approximation to the result.
This grading is economical near unit filling: a dilute gas on a long ring
can still require many hole roots. A different grading may be preferable
for that regime, but is not yet implemented.

`residual_norm` is the largest logarithmic residual divided by L, with each
level also retained in `level_residuals`. The default is 32 epsilon, not an
energy-error bound. A zero budget evaluates the seed only. Accepted updates
consume `max_iterations`; no-hole states use the existing XXX fixed-point
updates, while free states need none. Failure returns the final iterate,
energy estimate, residual, count, and a status: iteration limit, stalled,
or ill-conditioned. The CLI exits 2 for incomplete solves, 1 for invalid or
unsupported input, and 0 for convergence. Tolerances are never loosened and
arithmetic precision is never silently changed.

Tests cover native-precision analytic energies (including E=-4 for one
electron of each spin and E=-3-sqrt(5) for L=5, N_up=3, N_down=1), the XXX
reduction, free limits, original equations, Jacobian differences, budgets,
spin reversal, and singular Newton systems. Independent projected-Fock-space
energies and momenta cover all supported populations through L=8. Root and
equation checks also reach L=96. The ED oracle intentionally uses fp64;
separate analytic and equation checks run in all supported scalar types.

Next extensions require an explicit finite-ring branch audit for the missing
population parities; open boundaries and excited-state selection need their
own equations and labels. Neither shifted centered labels nor copying
Hubbard's finite-ring selection is a validated extension.
