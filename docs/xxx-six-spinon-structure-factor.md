# Six-spinon spectral points

`bethe/xxx_six_spinon_structure_factor.hpp` evaluates the longitudinal
six-spinon response of the infinite antiferromagnetic XXX chain at zero field.
It combines the [physical W6 density](xxx-spinon-density.md),
[sampling geometry](xxx-spinon-sampling.md) and
[bounded density recovery](xxx-spinon-density-recovery.md). It generates fresh
samples; no research archive or Python process is required.

This library interface is an initial integration checkpoint. It has a finite
rapidity cutoff and ordinary Monte Carlo errors. A successful result certifies
completion of its requested numerical checks, not convergence of the cutoff,
proposal resolution or sampling distribution. Installed command-line spectrum
and moment front ends are a subsequent integration step.

## Physical measure

The returned value is S6zz(k,omega), with J positive. W6 already contains the
longitudinal normalization and its single state-counting factorial. Every
sample includes all three momentum sections and both labelled solutions of the
energy constraint, weighted by the complete balanced proposal. No extra
factorial or 2 pi multiplier is applied to the spectrum.

The proposal cutoff means max(abs(beta)) is at most R for **all six** rapidities,
including the solved pair. With p=2 atan(exp(beta)), increasing R includes
progressively softer spinons near p=0 and pi. Cutting off these configurations
can remove weight within the bright continuum; it is not a high-energy cutoff.
At fixed k, changing exchange obeys S6zz(k,omega;J)=S6zz(k,omega/J;1)/J.

## A small calculation

This example uses MPFR, an explicit 384-bit geometry context, and density
controls that can recover from 128/192 to 224/320 bits. Include MPFR support
when [building Bethe](building.md). The settings demonstrate the interface;
R=8 and 512 draws are deliberately modest and do not give a converged physical
spectrum.

```cpp
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <uni20/async/tbb_scheduler.hpp>
#include <iostream>

int main()
{
    using namespace bethe::heisenberg;
    using R = uni20::mpreal;
    using Kind = detail::SpinonProposalKind;
    R const one(1, uni20::Precision::bits(384));
    R const pi = 4 * uni20::math::atan(one);
    detail::SpinonJointProposal<R> proposal(
        8 * one, one / 2,
        {{Kind::background, one / 2, one, one},
         {Kind::one_pair, one / 4, 4 * one, one / 2},
         {Kind::two_pairs, one / 4, 8 * one, one}});
    SpinonRecoveryPolicy recovery{
        "1e-7", "0",
        {{128, 192, {{{32, 12}, {64, 20}}}},
         {224, 320, {{{64, 20}, {128, 24}}}}},
        60};
    ThermodynamicSixSpinonStructureFactor<R, R> spectrum(proposal, recovery, one);
    uni20::async::TbbScheduler scheduler(2);
    auto pending = spectrum.at(pi / 2, 2 * one).evaluate_async(
        {.replicates = 8, .samples_per_replicate = 64, .seed = 42},
        {.scheduler = &scheduler, .max_in_flight = 2});
    auto result = pending.get_wait(scheduler);
    if (result.status != SixSpinonSpectrumStatus::estimated) {
        std::cerr << result.reason << '\n';
        return 1;
    }
    std::cout << uni20::format_scalar(*result.value, {.precision = 16}) << '\n';
    if (result.sampling_error)
        std::cout << uni20::format_scalar(*result.sampling_error, {.precision = 16}) << '\n';
}
```

The same `Point` offers synchronous `evaluate`, individual `replicate` calls,
and `observation` with an exact binary midpoint identity. Numerical state is
owned until completion; an originating temporary evaluator may be destroyed.
The supplied scheduler must outlive all admitted work, even if its Async handle
is dropped. Cancellation stops admission, drains work already admitted, and
retains its terminal records. No process-wide scheduler or precision is changed.

## Results and uncertainty

Each result records k, omega, J, R, proposal parameters, seed, replicate ids,
uniform precision, recovery policy, completed samples, supported roots, density
attempts and wall time. All complete replicate means are retained. This permits
later covariance analysis when frequencies share the same replicate streams.

The uncertainty fields have distinct meanings:

- `sampling_error` is the standard error estimated from equal-size independent
  replicate means. One replicate has no error estimate. A small number of
  replicates can underrepresent a rare, high-weight region.
- `contour_diagnostic` is the sampled, weighted contour-tail diagnostic from
  the density evaluations.
- `arithmetic_diagnostic` is the sampled, weighted density arithmetic screen,
  including the low/guard precision and control comparisons.

The latter two are numerical diagnostics, not rigorous confidence bounds.
They exclude sampling error, rapidity-tail bias and geometry discretization.
Geometry precision is independent of density recovery: retrying the density
cannot repair an unresolved energy surface. Compare whole calculations at
stronger geometry precision using the same sample identities. Geometry must
also retain the largest density guard precision. R=1040 requires much higher
geometry precision and stronger contour budgets than this R=8 example.

`estimated` requires every requested replicate to complete. Failed samples
are retained with their index, failure kind and density status, and prevent a
mean over the surviving replicates. Cancellation similarly yields no estimate.
`no_support_sampled` means that the draws missed the energy surface; it does
not report a zero spectral value. This is especially relevant near the very
narrow upper-edge region. Resolved exclusions outside the unrestricted continuum
and the upper-edge limit return zero; the lower threshold has no finite value.

## Checks and measured cost

Tests compare complete physical observations with an independent Python
construction using direct positive proposal PDFs, an acos energy root and
high-precision contour sums. They exercise every momentum section, soft
configurations and a frequency above the two-spinon edge. Other checks cover
whole-geometry precision growth from 384 to 512 bits, exchange scaling, identical
serial/concurrent replicate summaries, failures, cancellation and owned state.
The optional `tests/reference_xxx_six_spinon_spectrum.py` generates the compact
physical fixtures; mpmath is only needed to regenerate these references.

A fresh pilot with the example settings completed 512 draws, including 92
supported labelled roots and 184 density-control attempts. It gave
S6zz=0.00065061 with sampling standard error 0.00029002 (about 45%). Its sampled
contour and arithmetic diagnostics were 5.7e-37 and 2.6e-36. On the development
machine it used 43.23 seconds wall time and 78.12 seconds CPU time on two
workers. These numbers document cost and the separation of numerical checks
from sampling accuracy; they are not a reference thermodynamic spectrum.
