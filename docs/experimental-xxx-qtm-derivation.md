# Deriving the isotropic QTM reference formulas

[Usage, measures and limitations](experimental-xxx-qtm.md) ·
[Reference source](../experimental/reference/qtm/reference_xxx_qtm.py)

This page records the algebra implemented by the two pointwise references.
Equation numbers refer to the published papers: JMP Theorem 4 is Eq. (282),
whereas arXiv v3 calls it (284). JMP means
[Babenko et al., J. Math. Phys. 62, 041901 (2021)](../CITATIONS.md#babenko-gohmann-kozlowski-suzuki-2021);
PRL means [Babenko et al., Phys. Rev. Lett. 126, 210602 (2021)](../CITATIONS.md#babenko-gohmann-kozlowski-sirker-suzuki-2021).

Both return the symmetric Pauli-zz amplitude, with its integration measure
exposed separately. The independently normalized real-rapidity vertex control
uses the operator/Jacobian conversion in the usage guide. The computations
below do not establish a global contour deformation or a general proof of
higher-sector equivalence.

## Shared ordinary special functions

Writing $`G`$ for the Barnes function, JMP (275)–(276) give

```math
\widehat p(z)=\frac{\pi}{2}
-i\log\frac{\cosh[\pi(z+i/2)/2]}{\cosh[\pi(z-i/2)/2]},
```

```math
\Psi(z)=
\Gamma\!\left(\frac12-\frac{iz}{2}\right)
\Gamma\!\left(1-\frac{iz}{2}\right)
\frac{G(1-iz/2)^2}{G(1/2-iz/2)^2}.
```

Define

```math
Q=\frac12\sum_{j=1}^{\ell}[\widehat p(v_j)-\widehat p(u_j)],
\qquad
S(z)=i\,
\frac{\Gamma(1+iz/2)\Gamma(1/2-iz/2)}
     {\Gamma(1-iz/2)\Gamma(1/2+iz/2)}.
```

In products over $`U\ominus V`$, a hole has sign $`+1`$ and a particle
sign $`-1`$. A double product includes both signs, **including diagonal
self-pairs**. Logarithms are used to form Gamma/Barnes products; the final
amplitudes, rather than individual log branches, are compared.

## JMP: direct quadrature and the distinguished-column trap

JMP (279)–(280) use

```math
g(x,y)=\frac{\pi}{\sinh[\pi(y-x)]},\qquad
e(z)=\frac1z-\frac1{z-i},
```

```math
\phi^\sigma(t)=
\prod_{w\in U\ominus V}
\Gamma\!\left(\frac12-\frac{\sigma i(t-w)}2\right)
\Gamma\!\left(1+\frac{\sigma i(t-w)}2\right).
```

The two matrices are

```math
\Omega_{jk}=g(u_j,v_k)
-\int_{\mathbb R}\frac{dt}{2\pi i}
\frac{\phi^-(v_k)}{\phi^+(t)}g(u_j,t)e(t-v_k),
```

```math
\overline\Omega_{jk}=g(v_j,u_k)
-\int_{\mathbb R}\frac{dt}{2\pi i}
\frac{\phi^-(t)}{\phi^+(u_k)}g(v_j,t)e(u_k-t).
```

Construct all possible modified columns, not just the first:

```math
D_{jk}=\int_{\mathbb R}\frac{dt}{\pi}
\frac{\phi^-(v_k)}{\phi^+(t)}\frac{g(u_j,t)}{v_k-t},
\qquad
\overline D_{jk}=\int_{\mathbb R}\frac{dt}{\pi}
\frac{\phi^-(t)}{\phi^+(u_k)}\frac{g(v_j,t)}{u_k-t}.
```

For matrices $`B,C`$, let

```math
T(B,C)=\sum_{k=1}^{\ell}
\det(B\text{ with column }k\text{ replaced by column }k\text{ of }C).
```

The symmetric amplitude used for the pointwise comparison is

```math
A_\ell^{\rm JMP}
=(-1)^{\ell-1}\sin^2 Q
\left[\prod_{x,y\in U\ominus V}\Psi(x-y)\right]
T(\Omega,D)T(\overline\Omega,\overline D).
```

**Why not use Eq. (282) literally?** After Eq. (243), the authors use
symmetry of the outer integral to single out the first column of each
determinant. Eq. (282) therefore has measure $`1/[(\ell-1)!]^2`$ and only
those first-column replacements. Expressed with our common measure, that
literal integrand is

```math
A_{\ell,\mathrm{first}}^{\rm JMP}
=(-1)^{\ell-1}\ell^2\sin^2 Q
\left[\prod_{x,y\in U\ominus V}\Psi(x-y)\right]
\det\Omega^{[1]}\det\overline\Omega^{[1]}.
```

It is **not generally equal pointwise** to the symmetric PRL amplitude.
Averaging independently over the $`\ell`$ choices in each determinant
restores $`A_\ell^{\rm JMP}`$. This is an integration-symmetry convention,
not a defect in the paper. The code returns both versions, and a regression
test checks that they differ before symmetrization and agree afterwards.

Quadrature uses the real line, split at zero and all real rapidity parts.
`split=False` and a larger `maxdegree` provide alternative convergence
checks. Returned quadrature errors are mpmath estimates for the *matrix-entry
integrals*, not certified errors on the final determinant product.

## PRL: the isotropic limit cannot be taken entrywise alone

Set $`q=e^{-\gamma}`$, rescale the paper's rapidities by $`\gamma`$, and put
$`b=q^4=e^{-4\gamma}`$. A multiplicative parameter becomes

```math
q^c e^{2i\gamma z}=b^{c/4-iz/2}.
```

Thus the two basic hypergeometric functions in PRL (6) become ordinary
$`{}_{2\ell}F_{2\ell-1}(1)`$. Their parameters are implemented directly in
`hyper_parameters`. In both cases their excess is

```math
\sum_r b_r-\sum_r a_r=1.
```

However, $`\det M`$ and $`\det\widehat M`$ both vanish as $`\gamma\to0`$.
Just substituting $`q=1`$ gives zero and loses the leading answer.

### The first correction comes from the series tail

For either of these excess-one series, define

```math
\Phi={}_pF_{p-1}(\{a\};\{b\};1),
\qquad K=\frac{\prod_r\Gamma(b_r)}{\prod_r\Gamma(a_r)}.
```

The ordinary summand has tail $`K/n^2`$. At any **fixed** summation index,
the first-order change of its basic-hypergeometric counterpart cancels,
including the argument $`b=e^{-4\gamma}`$. Differentiating a finite
truncation therefore misses the actual leading correction.

With $`h=4\gamma`$, the transition region $`n\sim1/h`$ instead gives

```math
\Phi_q-\Phi
=Kh\int_0^\infty
\left[\frac{e^{-x}}{(1-e^{-x})^2}-\frac1{x^2}\right]dx+o(h)
=-\frac{Kh}{2}+o(h).
```

The integrand is finite at zero and integrable at infinity. Its primitive
is $`1/x-1/(e^x-1)`$, with endpoint values $`1/2`$ and zero.
Consequently the analytic first coefficient used here is

```math
\left.\frac{d\Phi_q}{d\gamma}\right|_{0^+}=-2K.
```

This matched-tail derivation applies at generic parameters away from poles;
it is not a uniform error bound near all degenerate rapidity configurations.
The test suite checks it directly against the unexpanded basic series.

### Singular determinants without inverses or finite differences

Let $`a(z)=-iz/2`$, $`a_{ij}=a(v_i-v_j)`$. The limit of PRL (7) is

```math
r_{ji}=\frac1{4a_{ij}(1+a_{ij})}
\prod_{u\in U}\frac{a(v_i-u)}{a(v_i-u)+1/2}
\prod_{r\ne i,j}\frac{a(v_i-v_r)+1/2}{a(v_i-v_r)}.
```

Its first derivative in $`\gamma`$ vanishes. The scattering factor also
has zero first derivative: the $`e^{i\gamma z}`$ phase in the finite-$`q`$
definition must be retained to cancel the first correction from q-Gamma.
Thus $`M_0,M_1`$ follow by using $`\Phi,-2K`$, respectively, in PRL (8)–(9).

The bar in that paper means **invert all $`P,H`$**, equivalently negate
all rapidities after rescaling; it does **not** mean complex conjugation.
The second matrix is obtained by replacing $`(U,V)`$ with $`(-V,-U)`$.

Set

```math
d_M=T(M_0,M_1),\qquad
d_{\widehat M}=T(\widehat M_0,\widehat M_1).
```

These are the coefficients of $`\gamma`$ in the respective determinants.
No inverse of a singular matrix and no numerical differentiation is needed.
At $`\ell=1`$, Gauss summation provides a separate closed-form check:

```math
\Phi_1=\frac{\Gamma(1/2+a)}{\sqrt{\pi}\Gamma(1+a)},\qquad
\Phi'_1=a\Phi_1,\qquad
(M_1)_{11}=i(v-u)\overline\Phi_1,\quad a=-i(v-u)/2.
```

### Limiting scalar prefactor

Define

```math
F(z)=\sinh^2(\pi z/2)\Psi(z)\Psi(-z),
\qquad s=\sum_j(v_j-u_j).
```

The resulting symmetric isotropic amplitude is

```math
A_\ell^{\rm PRL}
=d_Md_{\widehat M}
\left[\frac{\pi\Psi(0)}{2}\right]^{2\ell}
\frac{\pi^2\sin^2 Q}{\sinh^2(\pi s/2)}
\frac{\prod_{i\lt j}F(u_i-u_j)F(v_i-v_j)}
     {\prod_{i,j}F(u_i-v_j)}.
```

For checking constants, modular transformation gives

```math
\vartheta_1(\gamma z\mid e^{-a\gamma})
\sim B_a\sinh(\pi z/a),\qquad
B_a=2\sqrt{\frac{\pi}{a\gamma}}e^{-\pi^2/(4a\gamma)}.
```

The exponential factors cancel in PRL (10). Before the determinant product,
the remaining scalar scales as $`\gamma^{-2\ell-2}`$. The rescaled rapidity
measure contributes $`\gamma^{2\ell}`$, and the two determinant zeros
contribute $`\gamma^2`$. This accounts for the finite amplitude above.
Here $`k`$ is the paper's discrete branch label, not the physical momentum.
The staggered $`k=1`$ part vanishes in this isotropic limit; only $`k=0`$
is evaluated here, consistently with JMP (277).


## Removing the common contour offset analytically

The complex-partition audit makes the fixed-rapidity boundary
step simpler: **no offset extrapolation is needed for this purpose**.
This follows directly from the structure of the isotropic PRL formula,
not from a fit to the vertex result.

Write the symmetric amplitude as

```math
A_H(\eta)=K(\boldsymbol\beta,\eta)B_H(\boldsymbol\beta).
```

All factors in the PRL expression other than its explicit momentum sine
depend on rapidity differences: the hypergeometric parameters, scattering
ratios, determinant coefficients, pair products, and the sum
$`\sum(v-u)`$. They are unchanged by a common translation $`i\eta`$.
This is visible already in the structure of
[PRL Eqs. (6)–(11)](https://arxiv.org/html/2012.07378).
Here we use the explicit isotropic determinant coefficients derived in the
[preceding limiting algebra](#prl-the-isotropic-limit-cannot-be-taken-entrywise-alone).

The elementary identity

```math
e^{i\widehat p(z-i)}=e^{-i\widehat p(z)}
```

follows from the cosh ratio in $`\widehat p`$. Consequently the sine-squared
factor is the same for **every** assignment:

```math
K(\boldsymbol\beta,\eta)
=\sin^2\!\left[\frac12\sum_{j=1}^n
\widehat p\left(\frac{\beta_j}{\pi}+i\eta\right)\right].
```

Any $`2\pi`$ log-branch changes disappear in the sine squared. At the real
boundary this is simply

```math
K_0(\boldsymbol\beta)
=\sin^2\!\left[\sum_{j=1}^n\arctan(e^{\beta_j})\right],
\qquad A_H(0)=K_0B_H.
```

The audit computes $`B_H`$ by **omitting** the explicit sine factor from
the product, not by dividing $`A_H`$ by a possibly small $`K`$. This also
avoids a numerical $`0/0`$ at zero physical momentum, although other
coincident-rapidity/removable-singularity limits are not thereby solved.
It evaluates special functions at well-separated interior contour
offsets, then inserts the exact physical $`K_0`$.

This explains the offset dependence tested by the vertex audit. Near
$`\eta=0`$, $`K(-\eta)=K(\eta)^*`$; when the summed reduced coefficient is
real, its real correction starts at $`O(\eta^2)`$ and its imaginary part
at $`O(\eta)`$. There is no residual offset dependence to fit in $`B_H`$.
For nonzero $`K`$ and total, the cancellation condition is also independent
of the common offset:

```math
\frac{\sum_H|A_H(\eta)|}{|\sum_H A_H(\eta)|}
=\frac{\sum_H|B_H|}{|\sum_H B_H|}.
```

For the independent JMP route the corresponding coefficient is

```math
B_H^{\rm JMP}
=(-1)^{\ell-1}\!
\left[\prod_{x,y\in U\ominus V}\Psi(x-y)\right]
T(\Omega,D)T(\bar\Omega,\bar D).
```

Here $`T`$ is the sum of the distinguished-column determinants, not one
column chosen without symmetrization. The fixed-real-line JMP quadrature
and the PRL series are evaluated at **identical** $`U,V`$ for each
partition; their products $`B_H`$ are compared before any sum over $`H`$.
Their individual determinant coefficients need not be equal.

The [partition audit](../experimental/reference/qtm/audit_xxx_qtm_partitions.py) retains
the full complex amplitudes, reduced coefficients, explicit coordinates
and route-specific factors. It repeats the calculation at 30/42 decimal
digits and offsets $`0.25/0.1`$, changing the JMP quadrature splitting and
resolution as well. It also reconstructs **each** $`\eta=10^{-10}`$ complex
amplitude from $`K(\eta)B_H`$. These comparison inputs are generated afresh
by default, or read from an explicitly supplied, source-verified vertex audit.
The [usage guide](experimental-xxx-qtm.md#choose-an-audit) gives both commands.

This factorization is a fixed-rapidity algebraic statement within the
implemented isotropic formula. It is not a proof that an outer integration
contour may be shifted without residues, or that tails and collisions are
uniformly controlled. Nor does an individual complex $`B_H`$ represent a
positive Lehmann weight.
