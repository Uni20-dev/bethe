// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_weighted_value.hpp>
#include <bethe/detail/xxx_transition_rate.hpp>
#include <bethe/quadrature_accumulation.hpp>

namespace bethe::heisenberg
{
enum class TwoSpinonMomentStatus
{
  estimated,
  physical_zero,
  numerical_failure,
  cancelled
};

template <uni20::Real R> struct TwoSpinonMoments
{
    // m0,m1 (or I0,I1), then their integrated transition-kernel diagnostics.
    bethe::QuadratureReduction<R, 4> integral;
    TwoSpinonMomentStatus status = TwoSpinonMomentStatus::numerical_failure;
    std::optional<R> momentum;
    R exchange, gap_cutoff, kernel_tolerance;
    std::array<std::optional<R>, 2> sum_rule_denominators, sum_rule_fractions;
    std::size_t order = 0, max_nodes = 0, max_kernel_evaluations = 0, max_kernel_order = 0;
    bool unrestricted_static_divergence = false;
};

/// Direct longitudinal moments, integrated in the relative rapidity g>0.
/// The finite domain is g<gap_cutoff, with common shift unrestricted. Kernel
/// diagnostics exclude outer quadrature and the omitted large-gap tail.
template <uni20::Real R> class ThermodynamicTwoSpinonMoments {
    struct Kernel
    {
        R exchange, zero, one, pi;
        StructureFactorOptions<R> options;
        detail::XXXTransitionRate<R> rate;
        Kernel(R j, StructureFactorOptions<R> o)
            : exchange(std::move(j)), zero(bethe::detail::precision_zero(exchange)), one(zero + R{1}),
              pi(R{4} * uni20::math::atan(one)), options(std::move(o)), rate(options)
        {
          bethe::detail::require_same_precision(exchange, options.tolerance);
          if (!uni20::isfinite(exchange) || exchange <= zero) throw std::invalid_argument("require finite J>0");
        }
    };
    std::shared_ptr<Kernel const> kernel_;

  public:
    ThermodynamicTwoSpinonMoments(R exchange, StructureFactorOptions<R> options)
        : kernel_(std::make_shared<Kernel>(std::move(exchange), std::move(options)))
    {}

    class Query {
      public:
        using Result = TwoSpinonMoments<R>;
        std::size_t size() const { return panels_ * rule_->x.size(); }
        [[nodiscard]] uni20::async::Async<Result> evaluate_async(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          auto work = bethe::accumulate_quadrature_async<R, 4>(
              size(), [self = *this](std::size_t i, std::stop_token stop) { return self.observe(i, stop); },
              kernel_->one, execution);
          uni20::async::Async<Result> output;
          auto task = finish(work.read(), initial_, output.write());
          bethe::detail::schedule_owned(scheduler, std::move(task));
          return output;
        }
        Result evaluate(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          return evaluate_async(execution).get_wait(scheduler);
        }

      private:
        friend class ThermodynamicTwoSpinonMoments;
        std::shared_ptr<Kernel const> kernel_;
        std::shared_ptr<bethe::detail::GaussLegendreRule<R> const> rule_;
        Result initial_;
        std::size_t panels_ = 0;
        std::optional<R> log_sine_, log_cosine_;

        Query(std::shared_ptr<Kernel const> kernel, std::optional<R> momentum, std::size_t order, R cutoff,
              std::size_t max_nodes)
            : kernel_(std::move(kernel))
        {
          auto const& k = *kernel_;
          bethe::detail::require_same_precision(cutoff, k.one);
          if (order < 2 || order > 256 || !max_nodes || !uni20::isfinite(cutoff) || cutoff <= k.zero)
            throw std::invalid_argument("require order 2..256, positive gap cutoff and node budget");
          if (momentum)
          {
            bethe::detail::require_same_precision(*momentum, k.one);
            if (!uni20::isfinite(*momentum) || *momentum < k.zero || *momentum > R{2} * k.pi)
              throw std::invalid_argument("require momentum in [0,2*pi]");
          }
          initial_.momentum = momentum;
          initial_.exchange = k.exchange;
          initial_.gap_cutoff = cutoff;
          initial_.kernel_tolerance = k.options.tolerance;
          initial_.max_kernel_evaluations = k.options.max_evaluations;
          initial_.max_kernel_order = k.options.max_order;
          initial_.order = order;
          initial_.max_nodes = max_nodes;
          initial_.unrestricted_static_divergence = momentum && *momentum == k.pi;
          rule_ = std::make_shared<bethe::detail::GaussLegendreRule<R>>(
              bethe::detail::gauss_legendre<R>(order, uni20::precision_of(k.one)));
          if (momentum && (*momentum == k.zero || *momentum == R{2} * k.pi))
          {
            initial_.status = TwoSpinonMomentStatus::physical_zero;
            return;
          }
          // Locate ceil(cutoff/2) using integer comparisons at retained
          // precision, without narrowing an MPFR coordinate to a native type.
          std::size_t left = 1, right = max_nodes / order;
          if (!right || cutoff > R{2} * uni20::scalar_like(k.one, right))
            throw std::length_error("two-spinon moment quadrature exceeds node budget");
          while (left < right)
          {
            std::size_t const mid = left + (right - left) / 2;
            if (cutoff <= R{2} * uni20::scalar_like(k.one, mid))
              right = mid;
            else
              left = mid + 1;
          }
          panels_ = left;
          R const coefficient = R{2} * k.one / R{3} * (uni20::math::log(R{2} * k.one) - k.one / R{4});
          R log_first = uni20::math::log(k.exchange);
          if (momentum)
          {
            R const q = std::min(*momentum, R{2} * k.pi - *momentum);
            log_sine_ = uni20::math::log(uni20::math::sin(q / R{2}));
            if (q != k.pi) log_cosine_ = uni20::math::log(uni20::math::sin((k.pi - q) / R{2}));
            log_first += uni20::math::log(R{2} * k.one) + R{2} * *log_sine_;
          }
          else
            initial_.sum_rule_denominators[0] = k.one / R{4};
          initial_.sum_rule_denominators[1] = detail::spinon_weighted_value(coefficient, log_first);
        }

        bethe::QuadratureObservation<R, 4> observe(std::size_t index, std::stop_token stop) const
        {
          bethe::QuadratureObservation<R, 4> observation;
          if (stop.stop_requested())
          {
            observation.cancelled = true;
            observation.reason = "two-spinon moment cancelled";
            return observation;
          }
          auto const& k = *kernel_;
          std::size_t const panel = index / rule_->x.size(), j = index % rule_->x.size();
          R const left = R{2} * uni20::scalar_like(k.one, panel);
          R const width = std::min(R{2} * k.one, R(initial_.gap_cutoff - left));
          R const gap = left + width * (k.one + rule_->x[j]) / R{2};
          auto const rate = k.rate(gap / (R{2} * k.pi));
          observation.work_units = rate.evaluations;
          if (!rate.log_rate || !rate.error)
          {
            observation.reason = "two-spinon transition kernel did not converge";
            return observation;
          }
          R const log_cosh = gap / R{2} + uni20::math::log1p(uni20::math::exp(-gap)) - uni20::math::log(R{2} * k.one);
          R const log_base = *rate.log_rate - R{2} * log_cosh + uni20::math::log(width * rule_->w[j] / R{2});
          std::array<R, 2> log_weights;
          if (initial_.momentum)
          {
            R log_d = R{2} * (*log_sine_ - log_cosh);
            if (log_cosine_)
            {
              R const other = R{2} * *log_cosine_, maximum = std::max(log_d, other);
              log_d = maximum + uni20::math::log1p(uni20::math::exp(std::min(log_d, other) - maximum));
            }
            log_weights = {*log_sine_ - uni20::math::log(R{8} * k.pi) - log_d / R{2},
                           uni20::math::log(k.exchange) + R{2} * *log_sine_ - uni20::math::log(R{8} * k.one)};
          }
          else
            log_weights = {uni20::math::log(gap) - uni20::math::log(uni20::math::tanh(gap / R{2})) -
                               uni20::math::log(R{8} * k.pi * k.pi),
                           uni20::math::log(k.exchange) - uni20::math::log(R{16} * k.one)};
          try
          {
            std::array<R, 4> values;
            for (std::size_t i = 0; i < 2; ++i)
            {
              R const log_weight = log_base + log_weights[i];
              values[i] = detail::spinon_weighted_value(k.one, log_weight);
              values[i + 2] = detail::spinon_weighted_value(*rate.error, log_weight);
            }
            observation.values = std::move(values);
          }
          catch (SpinonDensityError const& e)
          {
            observation.reason = e.what();
          }
          return observation;
        }
        static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<bethe::QuadratureReduction<R, 4>> input,
                                              Result result, uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          result.integral = co_await input;
          if (result.integral.status == bethe::QuadratureStatus::cancelled)
            result.status = TwoSpinonMomentStatus::cancelled;
          else if (result.integral.status != bethe::QuadratureStatus::complete)
            result.status = TwoSpinonMomentStatus::numerical_failure;
          else
          {
            if (result.status != TwoSpinonMomentStatus::physical_zero) result.status = TwoSpinonMomentStatus::estimated;
            for (std::size_t i = 0; i < 2; ++i)
              if (result.sum_rule_denominators[i])
                result.sum_rule_fractions[i] = detail::spinon_weighted_value(
                    (*result.integral.sums)[i], -uni20::math::log(*result.sum_rule_denominators[i]));
          }
          destination = std::move(result);
        }
    };
    Query fixed_moments(R momentum, std::size_t order, R gap_cutoff, std::size_t max_nodes = 100000) const
    {
      return Query(kernel_, std::move(momentum), order, std::move(gap_cutoff), max_nodes);
    }
    Query global_moments(std::size_t order, R gap_cutoff, std::size_t max_nodes = 100000) const
    {
      return Query(kernel_, {}, order, std::move(gap_cutoff), max_nodes);
    }
};
} // namespace bethe::heisenberg
