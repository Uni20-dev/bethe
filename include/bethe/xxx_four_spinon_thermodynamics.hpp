// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/spinon_weighted_value.hpp>
#include <bethe/detail/xxx_four_spinon_common_shift.hpp>
#include <bethe/detail/xxx_four_spinon_geometry.hpp>
#include <bethe/quadrature_accumulation.hpp>
#include <bethe/xxx_spinon_density_recovery.hpp>
#include <functional>

namespace bethe::heisenberg
{
enum class FourSpinonIntegralStatus
{
  estimated,
  physical_zero,
  threshold,
  numerical_failure,
  cancelled
};
enum class FourSpinonObservable
{
  spectrum,
  fixed_moments,
  global_moments
};

template <uni20::Real R> struct SpinonIntegralEstimate
{
    R value, contour_diagnostic, arithmetic_diagnostic, shift_quadrature_diagnostic, shift_tail_bound;
};

template <uni20::Real R, std::size_t N> struct FourSpinonIntegralResult
{
    FourSpinonIntegralStatus status = FourSpinonIntegralStatus::numerical_failure;
    FourSpinonObservable observable;
    // Spectrum: N=1. Moments: N=2, ordered as zeroth and first moments.
    // An incomplete prefix never supplies an estimate.
    std::optional<std::array<SpinonIntegralEstimate<R>, N>> estimates;
    std::array<std::optional<R>, N> sum_rule_denominators, sum_rule_fractions;
    std::optional<R> momentum, omega, log_radius_cutoff;
    std::optional<std::array<R, 3>> gap_cutoffs;
    R exchange;
    std::size_t order = 0, max_nodes = 0, geometry_bits = 0, total_nodes = 0, completed_nodes = 0, density_attempts = 0;
    std::optional<std::size_t> failed_node;
    std::optional<SpinonRecoveryStatus> recovery_failure;
    SpinonRecoveryPolicy recovery;
    bool unrestricted_static_divergence = false;
    double wall_seconds = 0;
    std::string reason;
};

template <uni20::Real R, std::size_t N> struct FourSpinonRefinementResult
{
    FourSpinonIntegralStatus status = FourSpinonIntegralStatus::numerical_failure;
    // Coarse/fine order on the base domain, then coarse/fine on the extended
    // domain. These are numerical comparisons, not certified total errors.
    std::array<std::optional<FourSpinonIntegralResult<R, N>>, 4> controls;
    std::optional<std::array<SpinonIntegralEstimate<R>, N>> estimates;
    std::array<std::optional<R>, N> sum_rule_denominators, sum_rule_fractions;
    std::optional<std::array<R, N>> quadrature_difference, base_quadrature_difference, cutoff_change;
    std::string reason;
};

/// Physical longitudinal W4 integrals. Geometry/accumulation retains at least
/// the largest density guard precision. Recovery checks the kernel at these
/// coordinates, not the geometry; repeat the complete geometry at higher
/// precision for that separate check. Every query is a FINITE-domain estimate.
template <uni20::Real Geometry, uni20::Real Low, uni20::Real High = Low> class ThermodynamicFourSpinon {
  public:
    using G = Geometry;
    using SpectrumResult = FourSpinonIntegralResult<G, 1>;
    using MomentsResult = FourSpinonIntegralResult<G, 2>;

  private:
    struct Kernel
    {
        XXXSpinonDensityRecovery<Low, High> recovery;
        G exchange, zero, one, pi;
        Kernel(SpinonRecoveryPolicy policy, G j)
            : recovery(std::move(policy)), exchange(std::move(j)), zero(bethe::detail::precision_zero(exchange)),
              one(zero + G{1}), pi(G{4} * uni20::math::atan(one))
        {
          if (!uni20::isfinite(exchange) || exchange <= zero) throw std::invalid_argument("require finite J>0");
          if (uni20::numeric_limits<G>::digits(one) < recovery.policy().stages.back().guard_bits)
            throw std::invalid_argument("geometry must retain the largest density guard precision");
        }
    };
    std::shared_ptr<Kernel const> kernel_;

    template <std::size_t N> struct Node
    {
        std::array<G, 4> beta;
        std::array<G, N> log_weights;
        // Extra absolute measure errors (arithmetic, shift quadrature, shift
        // tails), kept logarithmic until multiplication by the physical W4.
        std::array<std::array<std::optional<G>, 3>, N> log_extra_weights;
    };

  public:
    ThermodynamicFourSpinon(SpinonRecoveryPolicy recovery, G exchange)
        : kernel_(std::make_shared<Kernel>(std::move(recovery), std::move(exchange)))
    {}

    template <std::size_t N> class Query {
      public:
        using Result = FourSpinonIntegralResult<G, N>;
        std::size_t size() const { return initial_.total_nodes; }

        [[nodiscard]] uni20::async::Async<Result> evaluate_async(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          auto work = bethe::accumulate_quadrature_async<G, 5 * N, std::optional<SpinonRecoveryStatus>>(
              initial_.total_nodes,
              [kernel = kernel_, node = node_](std::size_t i, std::stop_token stop) {
                return observe(*kernel, node, i, stop);
              },
              kernel_->one, execution);
          uni20::async::Async<Result> result;
          auto task = finish(work.read(), initial_, result.write());
          bethe::detail::schedule_owned(scheduler, std::move(task));
          return result;
        }

        Result evaluate(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          auto result = evaluate_async(execution);
          return result.get_wait(scheduler);
        }

      private:
        friend class ThermodynamicFourSpinon;
        std::shared_ptr<Kernel const> kernel_;
        std::function<Node<N>(std::size_t)> node_;
        Result initial_;
        Query(std::shared_ptr<Kernel const> kernel, FourSpinonObservable observable, std::size_t order,
              std::size_t max_nodes)
            : kernel_(std::move(kernel))
        {
          initial_.observable = observable;
          initial_.order = order;
          initial_.max_nodes = max_nodes;
          initial_.geometry_bits = uni20::numeric_limits<G>::digits(kernel_->one);
          initial_.exchange = kernel_->exchange;
          initial_.recovery = kernel_->recovery.policy();
        }
        using Observation = bethe::QuadratureObservation<G, 5 * N, std::optional<SpinonRecoveryStatus>>;
        using Reduction = bethe::QuadratureReduction<G, 5 * N, std::optional<SpinonRecoveryStatus>>;
        static Observation observe(Kernel const& k, std::function<Node<N>(std::size_t)> const& source,
                                   std::size_t index, std::stop_token stop)
        {
          Observation observation;
          if (stop.stop_requested())
          {
            observation.cancelled = true;
            observation.reason = "integration cancelled";
            return observation;
          }
          try
          {
            auto const point = source(index);
            auto recovered = k.recovery.recover(
                [&](auto const& context) {
                  std::vector<decltype(context.parse("0"))> beta;
                  for (auto const& b : point.beta)
                    beta.push_back(context.import(b));
                  return beta;
                },
                stop);
            observation.work_units = recovered.attempts.size();
            if (!recovered.accepted())
            {
              observation.metadata = recovered.status;
              observation.cancelled = recovered.status == SpinonRecoveryStatus::cancelled;
              observation.reason = recovered.reason;
              if (!recovered.attempts.empty()) observation.reason += ": " + recovered.attempts.back().reason;
              return observation;
            }
            auto const& d = *recovered.estimate;
            G const density = uni20::scalar_like(k.one, d.value);
            G const contour = uni20::scalar_like(k.one, d.contour_tail_error);
            G const arithmetic = uni20::scalar_like(k.one, *recovered.arithmetic_screen);
            std::array<G, 5 * N> values;
            for (std::size_t i = 0; i < N; ++i)
            {
              auto scaled = [&](G const& x) { return detail::spinon_weighted_value(x, point.log_weights[i]); };
              std::array<G, 3> extra{k.zero, k.zero, k.zero};
              for (std::size_t j = 0; j < extra.size(); ++j)
                if (point.log_extra_weights[i][j])
                  extra[j] = detail::spinon_weighted_value(density, *point.log_extra_weights[i][j]);
              values[5 * i] = scaled(density);
              values[5 * i + 1] = scaled(contour);
              values[5 * i + 2] = scaled(arithmetic) + extra[0];
              values[5 * i + 3] = extra[1];
              values[5 * i + 4] = extra[2];
            }
            observation.values = std::move(values);
          }
          catch (std::runtime_error const& error)
          {
            observation.reason = error.what();
          }
          return observation;
        }
        static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<Reduction> work, Result result,
                                              uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          auto summary = co_await work;
          result.completed_nodes = summary.completed_nodes;
          result.density_attempts = summary.work_units;
          result.failed_node = summary.failed_node;
          if (summary.failure_metadata) result.recovery_failure = *summary.failure_metadata;
          result.wall_seconds = summary.wall_seconds;
          if (summary.status != bethe::QuadratureStatus::complete)
          {
            result.status = summary.status == bethe::QuadratureStatus::cancelled
                                ? FourSpinonIntegralStatus::cancelled
                                : FourSpinonIntegralStatus::numerical_failure;
            result.reason = summary.reason;
          }
          else if (result.total_nodes || result.status == FourSpinonIntegralStatus::physical_zero)
          {
            std::array<SpinonIntegralEstimate<G>, N> values;
            auto const& fields = *summary.sums;
            for (std::size_t i = 0; i < N; ++i)
              values[i] = {fields[5 * i], fields[5 * i + 1], fields[5 * i + 2], fields[5 * i + 3], fields[5 * i + 4]};
            result.estimates = std::move(values);
            if (result.total_nodes)
            {
              result.status = FourSpinonIntegralStatus::estimated;
              result.reason = "finite-domain quadrature estimate; order and cutoff refinement are separate checks";
            }
            try
            {
              for (std::size_t i = 0; i < N; ++i)
                if (result.sum_rule_denominators[i])
                  result.sum_rule_fractions[i] = detail::spinon_weighted_value(
                      (*result.estimates)[i].value, -uni20::math::log(*result.sum_rule_denominators[i]));
            }
            catch (SpinonDensityError const& error)
            {
              result.estimates.reset();
              result.sum_rule_fractions = {};
              result.status = FourSpinonIntegralStatus::numerical_failure;
              result.reason = error.what();
            }
          }
          destination = std::move(result);
        }
    };

    template <std::size_t N> class RefinedQuery {
      public:
        using Result = FourSpinonRefinementResult<G, N>;
        [[nodiscard]] uni20::async::Async<Result> evaluate_async(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          uni20::async::Async<Result> output;
          auto task = calculate(factory_, execution, output.write());
          bethe::detail::schedule_owned(scheduler, std::move(task));
          return output;
        }
        Result evaluate(bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          auto pending = evaluate_async(execution);
          return pending.get_wait(scheduler);
        }

      private:
        friend class ThermodynamicFourSpinon;
        std::function<Query<N>(std::size_t)> factory_;
        explicit RefinedQuery(std::function<Query<N>(std::size_t)> factory) : factory_(std::move(factory)) {}
        static uni20::async::AsyncTask calculate(std::function<Query<N>(std::size_t)> factory,
                                                 bethe::ExecutionOptions execution,
                                                 uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          Result result;
          // Retain one geometry at a time. Each inner calculation uses the
          // same bounded admission; no nested wait or extra thread pool.
          for (std::size_t i = 0; i < result.controls.size(); ++i)
          {
            auto query = factory(i);
            auto pending = query.evaluate_async(execution);
            {
              auto value = co_await pending.read();
              result.controls[i] = std::move(value);
            }
            auto const& control = *result.controls[i];
            if (!control.estimates)
            {
              result.status = control.status;
              result.reason = "refinement control " + std::to_string(i) + " failed: " + control.reason;
              destination = std::move(result);
              co_return;
            }
          }
          std::array<G, N> quadrature, base_quadrature, cutoff;
          for (std::size_t i = 0; i < N; ++i)
          {
            auto value = [&](std::size_t control) { return (*result.controls[control]->estimates)[i].value; };
            quadrature[i] = uni20::math::abs(value(3) - value(2));
            base_quadrature[i] = uni20::math::abs(value(1) - value(0));
            cutoff[i] = value(3) - value(1);
            if (!uni20::isfinite(quadrature[i]) || !uni20::isfinite(base_quadrature[i]) || !uni20::isfinite(cutoff[i]))
              throw std::overflow_error("refinement difference exceeds exponent range");
          }
          result.estimates = result.controls[3]->estimates;
          result.sum_rule_denominators = result.controls[3]->sum_rule_denominators;
          result.sum_rule_fractions = result.controls[3]->sum_rule_fractions;
          result.status = result.controls[3]->status;
          result.quadrature_difference = std::move(quadrature);
          result.base_quadrature_difference = std::move(base_quadrature);
          result.cutoff_change = std::move(cutoff);
          result.reason =
              "four complete controls; differences do not bound the omitted domain or establish convergence";
          destination = std::move(result);
        }
    };

    RefinedQuery<1> refined_spectrum(G momentum, G omega, std::size_t coarse_order, std::size_t fine_order, G cutoff,
                                     G extended_cutoff, std::size_t max_nodes = 100000) const
    {
      check_refinement(coarse_order, fine_order, cutoff, extended_cutoff, max_nodes);
      check_momentum(momentum);
      (void)normalized_frequency(omega);
      return RefinedQuery<1>([model = *this, momentum, omega, coarse_order, fine_order, cutoff, extended_cutoff,
                              max_nodes](std::size_t i) {
        return model.spectrum(momentum, omega, i % 2 ? fine_order : coarse_order, i >= 2 ? extended_cutoff : cutoff,
                              max_nodes);
      });
    }

    RefinedQuery<2> refined_fixed_moments(G momentum, std::size_t coarse_order, std::size_t fine_order, G cutoff,
                                          G extended_cutoff, std::size_t max_nodes = 100000) const
    {
      check_refinement(coarse_order, fine_order, cutoff, extended_cutoff, max_nodes);
      check_momentum(momentum);
      return RefinedQuery<2>(
          [model = *this, momentum, coarse_order, fine_order, cutoff, extended_cutoff, max_nodes](std::size_t i) {
            return model.fixed_moments(momentum, i % 2 ? fine_order : coarse_order, i >= 2 ? extended_cutoff : cutoff,
                                       max_nodes);
          });
    }

    RefinedQuery<2> refined_global_moments(std::size_t coarse_order, std::size_t fine_order, std::array<G, 3> cutoffs,
                                           std::array<G, 3> extended_cutoffs, std::size_t max_nodes = 100000) const
    {
      for (int axis = 0; axis < 3; ++axis)
        check_refinement(coarse_order, fine_order, cutoffs[axis], extended_cutoffs[axis], max_nodes);
      // The finest order and widest domain dominate all four gap-node budgets.
      detail::FourSpinonGapGeometry<G>::validate_parameters(fine_order, extended_cutoffs, max_nodes);
      return RefinedQuery<2>([model = *this, coarse_order, fine_order, cutoffs, extended_cutoffs,
                              max_nodes](std::size_t i) {
        return model.global_moments(i % 2 ? fine_order : coarse_order, i >= 2 ? extended_cutoffs : cutoffs, max_nodes);
      });
    }

    Query<1> spectrum(G momentum, G omega, std::size_t order, G log_radius_cutoff, std::size_t max_nodes = 100000) const
    {
      auto const& k = *kernel_;
      check_momentum(momentum);
      G const normalized = normalized_frequency(omega);
      validate(log_radius_cutoff);
      detail::FourSpinonGeometry<G> geometry(order, log_radius_cutoff);
      if (!max_nodes) throw std::invalid_argument("require a positive node budget");
      Query<1> query(kernel_, FourSpinonObservable::spectrum, order, max_nodes);
      auto& r = query.initial_;
      r.momentum = momentum;
      r.omega = omega;
      r.log_radius_cutoff = log_radius_cutoff;
      G const q = std::min(momentum, G{2} * k.pi - momentum);
      G const lower = q == k.pi ? k.zero : k.pi / G{2} * uni20::math::sin(q);
      G const upper = G{2} * k.pi * uni20::math::cos(q / G{4});
      if (q == k.zero || normalized < lower || normalized > upper)
      {
        r.status = FourSpinonIntegralStatus::physical_zero;
        r.reason = "outside the unrestricted continuum or at zero momentum";
        return query;
      }
      if (normalized == lower || normalized == upper)
      {
        r.status = FourSpinonIntegralStatus::threshold;
        r.reason = "continuum threshold requires a limiting calculation";
        return query;
      }
      auto nodes = std::make_shared<std::vector<typename detail::FourSpinonGeometry<G>::Node>>();
      try
      {
        geometry.spectrum(momentum, normalized, [&](auto const& node) {
          if (nodes->size() == max_nodes) throw std::runtime_error("spectral quadrature exceeds node budget");
          nodes->push_back(node);
        });
      }
      catch (std::runtime_error const& error)
      {
        r.reason = error.what();
        return query;
      }
      if (nodes->empty())
      {
        r.reason = "no geometric observations; not a physical zero";
        return query;
      }
      r.total_nodes = nodes->size();
      G const log_exchange = uni20::math::log(k.exchange);
      query.node_ = [nodes, log_exchange](std::size_t i) {
        auto const& point = nodes->at(i);
        return Node<1>{point.beta, {point.log_measure - log_exchange}, {}};
      };
      return query;
    }

    Query<2> fixed_moments(G momentum, std::size_t order, G log_radius_cutoff, std::size_t max_nodes = 100000) const
    {
      auto const& k = *kernel_;
      check_momentum(momentum);
      validate(log_radius_cutoff);
      detail::FourSpinonGeometry<G> geometry(order, log_radius_cutoff);
      if (!max_nodes) throw std::invalid_argument("require a positive node budget");
      Query<2> query(kernel_, FourSpinonObservable::fixed_moments, order, max_nodes);
      auto& r = query.initial_;
      r.momentum = momentum;
      r.log_radius_cutoff = log_radius_cutoff;
      r.unrestricted_static_divergence = momentum == k.pi;
      if (momentum == k.zero || momentum == G{2} * k.pi)
      {
        r.status = FourSpinonIntegralStatus::physical_zero;
        r.reason = "inelastic zero-momentum response";
        return query;
      }
      G const q = std::min(momentum, G{2} * k.pi - momentum);
      G const sine = uni20::math::sin(q / G{2});
      r.sum_rule_denominators[1] = detail::spinon_weighted_value(first_sum_rule(), uni20::math::log(G{2} * k.one) +
                                                                                       G{2} * uni20::math::log(sine));
      auto nodes = std::make_shared<std::vector<typename detail::FourSpinonGeometry<G>::Node>>();
      try
      {
        geometry.moments(momentum, [&](auto const& node) {
          if (nodes->size() == max_nodes) throw std::runtime_error("moment quadrature exceeds node budget");
          nodes->push_back(node);
        });
      }
      catch (std::runtime_error const& error)
      {
        r.reason = error.what();
        return query;
      }
      r.total_nodes = nodes->size();
      G const log_measure = -uni20::math::log(G{2} * k.pi), log_exchange = uni20::math::log(k.exchange);
      query.node_ = [nodes, log_measure, log_exchange](std::size_t i) {
        auto const& point = nodes->at(i);
        G const weight = point.log_measure + log_measure;
        return Node<2>{point.beta, {weight, weight + uni20::math::log(point.energy) + log_exchange}, {}};
      };
      return query;
    }

    Query<2> global_moments(std::size_t order, std::array<G, 3> gap_cutoffs, std::size_t max_nodes = 100000) const
    {
      auto const& k = *kernel_;
      for (auto const& cutoff : gap_cutoffs)
        validate(cutoff);
      auto geometry = std::make_shared<detail::FourSpinonGapGeometry<G>>(order, gap_cutoffs, max_nodes);
      Query<2> query(kernel_, FourSpinonObservable::global_moments, order, max_nodes);
      query.initial_.sum_rule_denominators = {k.one / G{4}, first_sum_rule()};
      query.initial_.gap_cutoffs = gap_cutoffs;
      query.initial_.total_nodes = geometry->size();
      G const log_measure = -G{2} * uni20::math::log(G{2} * k.pi), log_exchange = uni20::math::log(k.exchange);
      query.node_ = [geometry, log_measure, log_exchange](std::size_t i) {
        auto const point = geometry->node(i);
        G const weight = point.log_measure + log_measure;
        Node<2> node{point.beta,
                     {weight + uni20::math::log(point.zeroth_shift.value),
                      weight + uni20::math::log(point.first_shift) + log_exchange},
                     {}};
        std::array<G, 3> errors{point.zeroth_shift.arithmetic_indicator, point.zeroth_shift.quadrature_difference,
                                point.zeroth_shift.omitted_tail_bound};
        for (std::size_t j = 0; j < errors.size(); ++j)
          if (errors[j] > G{0}) node.log_extra_weights[0][j] = weight + uni20::math::log(errors[j]);
        return node;
      };
      return query;
    }

  private:
    G first_sum_rule() const
    {
      auto const& k = *kernel_;
      return detail::spinon_weighted_value(G{2} * k.one / G{3} * (uni20::math::log(G{2} * k.one) - k.one / G{4}),
                                           uni20::math::log(k.exchange));
    }
    void check_momentum(G const& momentum) const
    {
      validate(momentum);
      if (momentum < kernel_->zero || momentum > G{2} * kernel_->pi)
        throw std::invalid_argument("require momentum in [0,2*pi]");
    }
    G normalized_frequency(G const& omega) const
    {
      validate(omega);
      G const normalized = omega / kernel_->exchange;
      if (!uni20::isfinite(normalized) || (omega != kernel_->zero && normalized == kernel_->zero))
        throw std::overflow_error("frequency exceeds normalized exponent range");
      return normalized;
    }
    void check_refinement(std::size_t coarse, std::size_t fine, G const& cutoff, G const& extended,
                          std::size_t max_nodes) const
    {
      validate(cutoff);
      validate(extended);
      if (!max_nodes) throw std::invalid_argument("require a positive node budget");
      if (coarse < 2 || fine <= coarse || fine > 256 || cutoff <= kernel_->zero || extended <= cutoff)
        throw std::invalid_argument("require increasing quadrature orders 2..256 and increasing positive cutoffs");
    }
    void validate(G const& value) const
    {
      bethe::detail::require_same_precision(value, kernel_->one);
      if (!uni20::isfinite(value)) throw std::invalid_argument("require finite observable coordinates");
    }
};
} // namespace bethe::heisenberg
