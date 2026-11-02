# Free-end XXX boundary dynamics

This experimental calculation produces local and end-to-end time correlations
for an even spin-half chain with free ends, zero field and
$`H=\sum_{j=1}^{N-1}\mathbf S_j\cdot\mathbf S_{j+1}`$. The exchange is one.
The intermediate states are the supported positive finite-real-root triplets.
Their weights are physically normalized; missing states are not compensated
by rescaling the retained spectrum.

This is a finite-chain calculation. Its real-root family is not a complete
spinon decomposition, and it does not include general complex strings or
arbitrary bulk-site matrix elements. Same-sector magnetization requires a
separate signed matrix-element calculation.

## Run a fresh calculation

Configure with `BETHE_BUILD_EXPERIMENTAL=ON` and `BETHE_BUILD_APPS=ON`.
From your build directory:

```sh
cmake --build . --target bethe-experimental-xxx-boundary-dynamics
./experimental/bethe-experimental-xxx-boundary-dynamics 16 --threads 2 --csv time.csv
./experimental/bethe-experimental-xxx-boundary-dynamics 8 --precision bits-160 --json result.json
```

The default grid has 101 points from time zero to 100. `--time-max` and
`--time-points` select another uniform grid. Real-valued controls are parsed
directly at the requested precision. Available arithmetic is fp64, native
long double and, when configured, binary128 or explicit MPFR precision.
`--references` prints the literature; `--help` lists every exportable table.

| Table | Content | Screen option |
| --- | --- | --- |
| `time` | Real, imaginary and absolute local/end-to-end correlations | Default primary table |
| `moments` | Retained moments, full local rules and omitted amounts | Default |
| `spectrum` | Accepted gaps, local weights and signed endpoint coefficients | `--spectrum` |
| `diagnostics` | Every requested label, including failed or unattempted states | `--diagnostics` |
| `roots` | Ground and attempted excited roots, in solver coordinates $`z=2\lambda`$ | `--roots` |
| `density` | Full ground-state boundary two-site density matrix | `--density` |

For example, `--csv-table spectrum=lines.csv` writes the weights without
also displaying them. `--json result.json` writes the declared tables together
with conventions, provenance, status and CPU time. Unknown table names are
errors. Exit status 2 denotes an incomplete or inconsistent numerical result;
failed requests remain in `diagnostics`, with unavailable values left absent.
A failed ground state produces no time curve. An incomplete excited family
retains its computed contribution with explicit failed records.

## Weights, endpoint signs and time convention

The [boundary-weight construction](experimental-xxx-boundary-form-factors.md)
gives the local raising weight between the singlet and highest-weight triplet.
SU(2) supplies the longitudinal weight. With $`\omega_n=E_n-E_0`$,

```math
\begin{aligned}
w_n&=\tfrac12|\langle n,S^z=1|S_1^+|0\rangle|^2,\\{}
C_{11}(t)&=\sum_{n\in\mathcal R}w_n e^{-i\omega_n t},\\{}
C_{1N}(t)&=\sum_{n\in\mathcal R}r_0r_nw_n e^{-i\omega_n t}.
\end{aligned}
```

Here $`\mathcal R`$ is the retained real-root family, and $`r_n=\pm1`$ is
the spatial-reflection eigenvalue. The observable is
$`C_{1N}(t)=\langle0|S_1^z(t)S_N^z(0)|0\rangle`$.
No periodic Fourier factor N or spectral factor 2 pi enters these finite
time sums. The full local norm is 1/4, but the end-to-end equal-time value is
a signed correlation, not 1/4.

For the positive integer labels $`I_j`$ and M down spins, the parity is

```math
r=(-1)^{NM+\sum_j I_j+M(M-1)/2}.
```

SU(2) descendants inherit this parity. The implementation computes the exponent
modulo two, avoiding integer products that might overflow. The ground parity
is $`(-1)^{N/2}`$. The tests compare the rule directly with reflected ABA
wavefunctions. Squared weights and this reflection identity determine the
endpoint cross coefficient; they do not supply arbitrary local operator phases.

## Full sum rules from the ground state

The complete local norm and first moment are

```math
m_0=\tfrac14,\qquad
m_1=\langle0|S_1^z(H-E_0)S_1^z|0\rangle
=-\langle S_1^zS_2^z\rangle=-\tfrac13\langle\mathbf S_1\cdot\mathbf S_2\rangle.
```

The first equality for $`m_1`$ fixes the sign and operator convention.
Only the adjacent bond enters the double commutator; SU(2) gives the last
equality. The boundary bond is computed from the ground roots alone, so the
full first-moment reference is independent of the restricted excitation scan.

Use standard roots $`u=z/2`$ and the open-chain Gaudin matrix G in the
[boundary-weight guide](experimental-xxx-boundary-form-factors.md). Then

```math
b_j=\frac{2u_j}{(u_j^2+1/4)^2},\qquad
Gx=b,\qquad
B=\langle\mathbf S_1\cdot\mathbf S_2\rangle=\tfrac14-\tfrac12b^Tx.
```

This follows by differentiating the integrable boundary-impurity coupling
$`J_b=1/(1+\theta^2)`$ at the homogeneous point and applying
Hellmann–Feynman. The reference is Yupeng Wang, *Exact solution of the open
Heisenberg chain with two impurities*, Phys. Rev. B **56**, 14045 (1997):
[published article](https://doi.org/10.1103/PhysRevB.56.14045),
[open preprint](https://arxiv.org/abs/cond-mat/9805253).
The code uses one linear solve through ordinary Uni20 dispatch. It checks the
original root residual and the normwise backward error of this solve, without
finite differencing or selecting a backend manually.

SU(2) also fixes the full two-site density matrix. Put $`c=B/3`$; in the
ordered basis up-up, up-down, down-up, down-down,

```math
\rho_{12}=\begin{pmatrix}
1/4+c&0&0&0\\{}
0&1/4-c&2c&0\\{}
0&2c&1/4-c&0\\{}
0&0&0&1/4+c
\end{pmatrix}.
```

Its singlet eigenvalue is $`1/4-B`$ and each triplet eigenvalue is
$`1/4+B/3`$. The reported matrix and probabilities are not clamped.

For exact retained weights and energies, $`1/4-m_0^{\mathcal R}`$ bounds the
absolute omitted local or endpoint correlation at every time. Numerical root
and form-factor errors remain separate. `omitted_weight()` and
`omitted_first_moment()` return no value when a computed retained moment exceeds
its full reference, even by roundoff; they never manufacture a zero bound.
The overall acceptance check allows a working-precision margin of
4096 N epsilon for those comparisons. This margin is not a certified error bar.

## Library and Async use

Link the opt-in `bethe::experimental` target. The synchronous wrapper and the
composable operation use the same owned calculation:

```cpp
#include <bethe/experimental/xxx_boundary_dynamics.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace boundary = bethe::experimental::xxx::open;
uni20::async::TbbScheduler scheduler(2);
auto pending = boundary::boundary_real_triplet_spectrum_async<double>(
    8, {}, 10000, {&scheduler});
auto const& result = pending.get_wait(scheduler);
if (!result.converged()) throw std::runtime_error("inspect boundary diagnostics");
auto local = bethe::local_time_correlation<double>(result.lines, 10.0);
auto endpoint = bethe::cross_time_correlation<double>(result.end_to_end_lines, 10.0);
```

Each worker solves and evaluates one state using shared immutable ground
preparation. Results are collected in lexicographic label order, independent
of scheduling. The finite family has M=N/2−1 roots chosen from integers 1
through N−M, hence $`\binom{N/2+1}{2}`$ candidates. `max_candidates` is a
preflight bound, not an early termination with silently missing states.

The operation owns its solver options and numerical inputs. The borrowed
scheduler must outlive it. Cancellation or an exception stops further admission,
drains admitted kernels and reports an error through the final result; dropping
the handle alone does not cancel work. No `run_all()` is needed.

For MPFR, construct `SolverOptions<uni20::mpreal>` with an explicit
`uni20::Precision`, and construct time values at that same precision. Neither
submission nor evaluation changes an ambient precision default. The ground-only
`boundary_two_site_density(N, ground)` API also works without an excitation scan.

## Validation and limits

Small-chain ABA and ED controls check individual weights, reflection parity,
the full two-site density, and local/end-to-end time correlations through
time 100. Independent perturbed-boundary Hamiltonians check the bond response
by Hellmann–Feynman. The dimer gives
$`C_{11}(t)=e^{-it}/4`$ and $`C_{1N}(t)=-e^{-it}/4`$.
Explicit MPFR 128/256-bit controls check precision propagation, while serial,
concurrent, held-work, failed-ground and cancellation cases check ownership
and failure accounting. CLI tests cover normalization and export contracts.

Long-time phase accuracy still depends on energy accuracy. Increasing the
number of time samples does not improve the root solutions or restore omitted
complex-root states. The complete two-site ground density and full local
sum rules should not be confused with the retained dynamical spectrum.

The ground-root density and full first-moment reference are evaluated independently
of spectral form-factor preparation. If spectral preparation fails but the density
solve succeeds, `moments` still exports the full reference; retained moments and
omissions are absent, not zero. The time table remains empty and the exit code is 2.
