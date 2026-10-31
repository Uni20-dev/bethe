// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/solver.hpp>

namespace bethe::heisenberg::detail
{
/// Importance proposal on the fixed-momentum six-spinon surface. Choose one
/// of 30 labelled charts: one eliminated momentum, one native rapidity, four
/// joint-proposal rapidities. Unsupported draws are zero contributions, never
/// redrawn. The reciprocal COMPLETE chart intensity includes the momentum
/// constraint Jacobian and dp/d beta, but no W6, factorial or 1/(2*pi).
template <uni20::Real R> class SpinonFixedMomentumProposal {
  public:
    using Proposal = SpinonJointProposal<R>;
    struct Sample
    {
        bool inside = false;
        std::array<R, 6> beta;
        R log_momentum_measure, log_energy; // energy at J=1
        std::size_t component = 0, pairing = 0, eliminated = 0, native = 0;
        int section = -1;
    };

    SpinonFixedMomentumProposal(Proposal proposal, R momentum)
        : proposal_(std::move(proposal)), momentum_(std::move(momentum)),
          zero_(uni20::scalar_like(proposal_.cutoff(), 0)), one_(uni20::scalar_like(zero_, 1)),
          pi_(R{4} * uni20::math::atan(one_)), lower_(R{2} * uni20::math::atan(uni20::math::exp(-proposal_.cutoff()))),
          upper_(pi_ - lower_), tolerance_(R{64} * uni20::numeric_limits<R>::epsilon(one_) * pi_)
    {
      bethe::detail::require_same_precision(momentum_, one_);
      if (!uni20::isfinite(momentum_) || momentum_ < zero_ || momentum_ >= R{2} * pi_)
        throw std::invalid_argument("require momentum in [0,2*pi)");
      if (!(zero_ < lower_ && lower_ < upper_ && upper_ < pi_))
        throw std::runtime_error("moment cutoff unresolved; increase geometry precision");
    }

    R const& momentum() const { return momentum_; }
    Proposal const& proposal() const { return proposal_; }

    std::array<R, 30> chart_log_intensities(std::array<R, 6> const& beta) const
    {
      auto const global = proposal_.global_chart_log_densities(beta);
      std::array<R, 6> native, cosh;
      R total_cosh = zero_;
      for (int i = 0; i < 6; ++i)
      {
        native[i] = proposal_.native_log_density(beta[i]);
        cosh[i] = log_cosh(beta[i]);
        total_cosh += cosh[i];
      }
      std::array<R, 30> result;
      std::size_t a = 0;
      for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j, ++a)
        {
          // global[a] = log(q4(other four) g_i g_j). Remove the
          // eliminated coordinate's g, and convert the five free PDFs to dp.
          result[2 * a] = global[a] - native[i] + total_cosh - cosh[i];
          result[2 * a + 1] = global[a] - native[j] + total_cosh - cosh[j];
        }
      return result;
    }

    R log_momentum_measure(std::array<R, 6> const& beta) const
    {
      auto const h = chart_log_intensities(beta);
      return uni20::math::log(R{30} * one_) - log_sum(h);
    }

    Sample draw(std::array<R, 6> const& u) const
    {
      for (auto const& x : u)
      {
        bethe::detail::require_same_precision(x, one_);
        if (!uni20::isfinite(x) || !(zero_ < x && x < one_))
          throw std::invalid_argument("require open uniform coordinate");
      }
      R const position = R{30} * u[0];
      std::size_t chart = 0;
      while (chart + 1 < 30 && position >= R(chart + 1))
        ++chart;
      R const remainder = position - R(chart);
      if (!(zero_ < remainder && remainder < one_))
        throw std::runtime_error("moment proposal selector lost open subinterval resolution");
      auto const free = proposal_.draw({remainder, u[1], u[2], u[3], u[4]});
      Sample result;
      result.beta.fill(zero_);
      result.log_momentum_measure = result.log_energy = zero_;
      result.component = free.component;
      result.pairing = free.pairing;
      result.eliminated = chart / 5;
      result.native = chart % 5;
      if (result.native >= result.eliminated) ++result.native;
      result.beta[result.native] = proposal_.native_inverse(u[5]);
      std::size_t next = 0;
      for (std::size_t i = 0; i < 6; ++i)
        if (i != result.eliminated && i != result.native) result.beta[i] = free.beta[next++];
      bethe::detail::CompensatedSum<R> total(zero_);
      for (std::size_t i = 0; i < 6; ++i)
        if (i != result.eliminated)
        {
          R const distance = R{2} * uni20::math::atan(uni20::math::exp(-uni20::math::abs(result.beta[i])));
          R const p = result.beta[i] <= zero_ ? distance : pi_ - distance;
          if (!(lower_ < p && p < upper_)) throw std::runtime_error("free momentum lost cutoff resolution");
          total.add(p);
        }
      // The open interval has width < pi, so at most one section is supported
      // for these five free momenta. Examine all three rather than folding
      // momenta with a discontinuous floating-point remainder operation.
      for (int section = 0; section < 3; ++section)
      {
        R const p = momentum_ + R(2 * section) * pi_ - total.value();
        if (p < lower_ - tolerance_ || p > upper_ + tolerance_) continue;
        if (p <= lower_ + tolerance_ || p >= upper_ - tolerance_)
          throw std::runtime_error("fixed-momentum boundary unresolved; increase geometry precision");
        if (result.inside) throw std::runtime_error("ambiguous fixed-momentum section");
        R const b = p <= pi_ / R{2} ? uni20::math::log(uni20::math::tan(p / R{2}))
                                    : -uni20::math::log(uni20::math::tan((pi_ - p) / R{2}));
        if (!uni20::isfinite(b) || !(uni20::math::abs(b) < proposal_.cutoff()))
          throw std::runtime_error("solved rapidity lost cutoff resolution");
        result.beta[result.eliminated] = b;
        result.inside = true;
        result.section = section;
      }
      if (result.inside)
      {
        result.log_momentum_measure = log_momentum_measure(result.beta);
        std::array<R, 6> log_sines;
        for (int i = 0; i < 6; ++i)
          log_sines[i] = -log_cosh(result.beta[i]);
        result.log_energy = uni20::math::log(pi_ / R{2}) + log_sum(log_sines);
      }
      return result;
    }

  private:
    R log_cosh(R const& x) const
    {
      R const a = uni20::math::abs(x);
      return a + uni20::math::log1p(uni20::math::exp(-R{2} * a)) - uni20::math::log(R{2} * one_);
    }
    R log_sum(std::span<R const> terms) const
    {
      R const largest = *std::max_element(terms.begin(), terms.end());
      R sum = zero_;
      for (auto const& v : terms)
        sum += uni20::math::exp(v - largest);
      R const result = largest + uni20::math::log(sum);
      if (!uni20::isfinite(result)) throw std::runtime_error("fixed-momentum intensity exceeds exponent range");
      return result;
    }
    Proposal proposal_;
    R momentum_, zero_, one_, pi_, lower_, upper_, tolerance_;
};
} // namespace bethe::heisenberg::detail
