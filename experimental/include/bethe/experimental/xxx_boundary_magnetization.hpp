// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/newton.hpp>
#include <bethe/experimental/xxx_boundary_dynamics.hpp>

namespace bethe::experimental::xxx::open
{
enum class BoundarySpinStatus
{
  converged,
  roots_unconverged,
  linear_solve_failed,
  precision_limit
};
template <uni20::Real R> class BoundarySpinFormFactors;

/// Immutable, calculator-bound root products and raising gauge. Copies share
/// owned data; a moved-from or rejected preparation is not ready.
template <uni20::Real R> class BoundarySpinState {
    friend class BoundarySpinFormFactors<R>;
    struct Snapshot
    {
        explicit Snapshot(R const& zero)
            : log_gaudin(zero), log_s(zero), log_p(zero), p_sum(zero), amplitude(zero), diagonal(zero), pivot(zero)
        {}
        std::shared_ptr<int const> owner;
        std::vector<R> roots, self_log_difference, self_log_shifted, kernel;
        std::vector<uni20::complex<R>> self_phase;
        R log_gaudin, log_s, log_p, p_sum, amplitude, diagonal, pivot;
        std::optional<R> solve_backward_error;
        int parity = 1, raising_sign = 0;
        BoundarySpinStatus status = BoundarySpinStatus::precision_limit;
    };
    std::shared_ptr<Snapshot const> data_;
    BoundarySpinStatus rejected_status_ = BoundarySpinStatus::precision_limit;
    bool attempted_ = false;

  public:
    bool attempted() const { return attempted_; }
    bool ready() const { return data_ && data_->status == BoundarySpinStatus::converged; }
    BoundarySpinStatus status() const { return data_ ? data_->status : rejected_status_; }
    std::span<R const> aba_rapidities() const
    {
      return ready() ? std::span<R const>(data_->roots) : std::span<R const>{};
    }
    std::optional<R> raising_amplitude() const
    {
      return ready() ? std::optional<R>(data_->amplitude) : std::optional<R>{};
    }
    std::optional<R> diagonal() const { return ready() ? std::optional<R>(data_->diagonal) : std::optional<R>{}; }
    std::optional<R> solve_backward_error() const { return data_ ? data_->solve_backward_error : std::optional<R>{}; }
    std::optional<R> pivot_margin() const { return data_ ? std::optional<R>(data_->pivot) : std::optional<R>{}; }
    std::optional<int> reflection_parity() const
    {
      return ready() ? std::optional<int>(data_->parity) : std::optional<int>{};
    }
};

template <uni20::Real R> struct BoundarySpinElement
{
    std::optional<R> value, pivot_margin;
    BoundarySpinStatus status = BoundarySpinStatus::precision_limit;
    bool converged() const { return status == BoundarySpinStatus::converged && value.has_value(); }
};

/// Same-sector S_1^z in the gauge with positive singlet-to-triplet S_1^+
/// amplitudes. Retains the directed boundary product: squared weights cannot
/// check it. Only the prepared calculator (or a copy) may combine its states.
template <uni20::Real R> class BoundarySpinFormFactors {
    std::size_t sites_;
    BoundaryFormFactors<R> raising_;
    std::shared_ptr<int const> owner_ = std::make_shared<int const>(0);

  public:
    BoundarySpinFormFactors(std::size_t sites, model::RealState<R> const& ground)
        : sites_(sites), raising_(sites, ground)
    {}
    bool ready() const { return owner_ && raising_.ready(); }

    [[nodiscard]] BoundarySpinState<R> prepare(model::RealState<R> const& state) const
    {
      using C = uni20::complex<R>;
      using uni20::math::abs;
      using uni20::math::log;
      BoundarySpinState<R> result;
      result.attempted_ = true;
      auto const raising = raising_(state); // Recheck equations; never trust cached convergence.
      if (!ready() || !raising.converged() || abs(raising.determinant_sign) != 1)
      {
        if (raising.status == BoundaryFormFactorStatus::roots_unconverged)
          result.rejected_status_ = BoundarySpinStatus::roots_unconverged;
        return result;
      }
      R const zero = *raising_.working_zero(), one = uni20::scalar_like(zero, 1);
      auto saved = std::make_shared<typename BoundarySpinState<R>::Snapshot>(zero);
      result.data_ = saved;
      auto& s = *saved;
      s.owner = owner_;
      s.amplitude = uni20::math::sqrt(*raising.weight);
      s.raising_sign = raising.determinant_sign;
      s.parity = real_root_reflection_parity(sites_, state.quantum_numbers);
      for (auto const& z : state.rapidities)
        s.roots.push_back(z / 2);
      auto const& u = s.roots;
      auto const m = u.size();
      auto const matrix = boundary_detail::gaudin_matrix<R>(sites_, u, zero);
      auto const norm = bethe::detail::log_determinant<R>(matrix, m, zero);
      s.pivot = std::min(*raising.pivot_margin, norm.pivot_margin);
      if (!norm.log_absolute) return result;
      s.log_gaudin = *norm.log_absolute;
      s.kernel.resize(m * m, zero);
      s.self_log_difference.resize(m, zero);
      s.self_log_shifted.resize(m, zero);
      s.self_phase.resize(m, C(zero, zero));
      bethe::detail::CompensatedSum<R> log_s(zero), log_p(zero), p_sum(zero);
      std::vector<R> rhs(m, zero);
      for (std::size_t i = 0; i < m; ++i)
      {
        R const denominator = u[i] * u[i] + one / 4;
        log_s.add(log(4 * u[i] * u[i]));
        log_p.add(log(denominator));
        p_sum.add(one / denominator);
        rhs[i] = 2 * u[i] / (denominator * denominator);
        bethe::detail::CompensatedSum<R> difference(zero), shifted(zero);
        C phase(one, zero);
        for (std::size_t j = 0; j < m; ++j)
        {
          R const d = u[i] - u[j], v = u[i] + u[j];
          s.kernel[i * m + j] = -2 * one / (one + d * d) + 2 * one / (one + v * v);
          C const q(d * v - one, -2 * u[i]);
          phase *= q / abs(q);
          shifted.add(boundary_detail::pair_log(u[i], u[j]) / 2);
          if (i != j) difference.add(log(abs(d * v)));
          if (j < i) log_s.add(2 * log(abs(d * v)) - boundary_detail::pair_log(u[i], u[j]));
        }
        s.self_log_difference[i] = difference.value();
        s.self_log_shifted[i] = shifted.value();
        s.self_phase[i] = phase;
      }
      s.log_s = log_s.value();
      s.log_p = log_p.value();
      s.p_sum = p_sum.value();
      auto response = rhs;
      // Shared square solve through normal Uni20 dispatch, at the root precision.
      if (m && !bethe::detail::newton_step(matrix, response, uni20::precision_of(zero)))
      {
        s.status = BoundarySpinStatus::linear_solve_failed;
        return result;
      }
      R matrix_norm = zero, rhs_norm = zero, response_norm = zero, residual_norm = zero;
      bethe::detail::CompensatedSum<R> contraction(zero);
      for (std::size_t i = 0; i < m; ++i)
      {
        if (!uni20::isfinite(response[i])) return result;
        contraction.add(u[i] * response[i]);
        rhs_norm = std::max(rhs_norm, abs(rhs[i]));
        response_norm = std::max(response_norm, abs(response[i]));
        bethe::detail::CompensatedSum<R> row_norm(zero), residual(zero);
        residual.add(-rhs[i]);
        for (std::size_t j = 0; j < m; ++j)
        {
          row_norm.add(abs(matrix[i * m + j]));
          residual.add(matrix[i * m + j] * response[j]);
        }
        matrix_norm = std::max(matrix_norm, row_norm.value());
        residual_norm = std::max(residual_norm, abs(residual.value()));
      }
      R const scale = matrix_norm * response_norm + rhs_norm;
      if (m && (!uni20::isfinite(scale) || scale <= zero || !uni20::isfinite(residual_norm))) return result;
      s.solve_backward_error = m ? R(residual_norm / scale) : zero;
      s.diagonal = one / 2 - 2 * contraction.value();
      R const gate = 4096 * uni20::scalar_like(zero, sites_) * uni20::numeric_limits<R>::epsilon(zero);
      if (!uni20::isfinite(s.diagonal) || abs(s.diagonal) > one / 2 + gate || *s.solve_backward_error > gate ||
          !uni20::isfinite(s.log_s) || !uni20::isfinite(s.log_p) || !uni20::isfinite(s.p_sum))
        return result;
      s.status = BoundarySpinStatus::converged;
      return result;
    }

    [[nodiscard]] BoundarySpinElement<R> operator()(BoundarySpinState<R> const& left,
                                                    BoundarySpinState<R> const& right) const
    {
      using C = uni20::complex<R>;
      using uni20::math::abs;
      using uni20::math::exp;
      using uni20::math::log;
      BoundarySpinElement<R> result;
      if (!ready() || !left.ready() || !right.ready()) return result;
      auto const& x = *left.data_;
      auto const& y = *right.data_;
      if (x.owner != owner_ || y.owner != owner_)
        throw std::invalid_argument("boundary spin states belong to another preparation");
      R const zero = *raising_.working_zero(), one = uni20::scalar_like(zero, 1);
      result.pivot_margin = std::min(x.pivot, y.pivot);
      auto const m = x.roots.size();
      if (x.roots == y.roots)
      {
        result.value = x.diagonal;
        result.status = BoundarySpinStatus::converged;
        return result;
      }
      std::vector<R> matrix(m * m, zero);
      R const pole_gate = 256 * uni20::scalar_like(zero, sites_) * uni20::numeric_limits<R>::epsilon(zero);
      R const dp = y.p_sum - x.p_sum, log_q = y.log_p - x.log_p, largest_y = y.roots.back();
      bethe::detail::CompensatedSum<R> row_scales(zero);
      for (std::size_t j = 0; j < m; ++j)
      {
        R const v = y.roots[j], p = one / (v * v + one / 4);
        int sign = 1;
        C phase(one, zero);
        bethe::detail::CompensatedSum<R> differences(zero), shifted(zero);
        for (auto const& u : x.roots)
        {
          R const d = v - u, s = v + u;
          if (abs(d) <= pole_gate * (one + abs(v) + abs(u))) return result;
          if (d < zero) sign = -sign;
          differences.add(log(abs(d * s)));
          C const q(d * s - one, -2 * v);
          phase *= q / abs(q);
          shifted.add(boundary_detail::pair_log(v, u) / 2);
        }
        if ((m - j - 1) % 2) sign = -sign;
        phase *= uni20::conj(y.self_phase[j]);
        R const log_diagonal =
            log(4 * v) + y.self_log_difference[j] - differences.value() + shifted.value() - y.self_log_shifted[j];
        R const rank = -8 * v * p * (dp - p);
        R const log_rank = rank == zero ? zero : log(abs(rank)) + log_q + log(largest_y);
        R const scale = std::max({zero, log_diagonal, log_rank});
        if (!uni20::isfinite(scale)) return result;
        row_scales.add(scale);
        R const diagonal = sign * phase.imag() * exp(log_diagonal - scale);
        R const rank_scaled = rank == zero ? zero : R((rank > zero ? 1 : -1) * exp(log_rank - scale));
        for (std::size_t k = 0; k < m; ++k)
          matrix[j * m + k] =
              y.kernel[j * m + k] * exp(-scale) + rank_scaled * (y.roots[k] / largest_y) + (j == k ? diagonal : zero);
      }
      auto const det = bethe::detail::log_determinant<R>(std::move(matrix), m, zero);
      result.pivot_margin = std::min(*result.pivot_margin, det.pivot_margin);
      if (!det.log_absolute || !det.phase) return result;
      bethe::detail::CompensatedSum<R> logarithm(zero);
      logarithm.add(*det.log_absolute + row_scales.value());
      logarithm.add(-log(2 * one) - (x.log_gaudin + y.log_gaudin) / 2);
      logarithm.add((x.log_s - y.log_s) / 2 + x.log_p - y.log_p);
      int const sign =
          (*det.phase > zero ? 1 : -1) * x.raising_sign * y.raising_sign * x.parity * y.parity * (m % 2 ? -1 : 1);
      R const value = sign * exp(logarithm.value());
      R const gate = 4096 * uni20::scalar_like(zero, sites_) * uni20::numeric_limits<R>::epsilon(zero);
      if (!uni20::isfinite(value) || value == zero || abs(value) > one / 2 + gate) return result;
      result.value = value;
      result.status = BoundarySpinStatus::converged;
      return result;
    }
};

template <uni20::Real R> struct BoundaryMagnetizationTerm
{
    R frequency, coefficient; // Real signed cosine coefficient, not a probability.
    int opposite_end_parity = 1;
};
template <uni20::Real R> struct BoundaryMagnetizationPair
{
    std::size_t left_state = 0, right_state = 0; // One-based canonical spectrum identifiers.
    std::optional<BoundarySpinElement<R>> element;
    std::optional<BoundaryMagnetizationTerm<R>> term;
    std::string diagnostic;
};
template <uni20::Real R> struct BoundaryMagnetization
{
    explicit BoundaryMagnetization(BoundarySpectrum<R> input) : spectrum(std::move(input)) {}
    BoundarySpectrum<R> spectrum;
    std::vector<BoundarySpinState<R>> states;
    std::vector<BoundaryMagnetizationPair<R>> pairs;
    std::optional<R> projected_norm;
    bool matrix_ready = false;
    bool converged() const { return matrix_ready && spectrum.converged(); }
    std::optional<R> omitted_norm() const
    {
      if (!matrix_ready || !projected_norm || !uni20::isfinite(*projected_norm) || *projected_norm < 0 ||
          *projected_norm > 1)
        return {};
      return uni20::scalar_like(*projected_norm, 1) - *projected_norm;
    }
    /// Uniform bound for the unnormalized projection; numerical error is separate.
    std::optional<R> truncation_bound() const
    {
      if (auto e = omitted_norm()) return uni20::math::sqrt(*e * (uni20::scalar_like(*e, 4) - 3 * *e)) / 2;
      return {};
    }
    [[nodiscard]] R at(R const& time, bool opposite_end = false) const
    {
      if (!converged() || !projected_norm)
        throw std::logic_error("boundary magnetization has unresolved states or pairs");
      bethe::detail::require_same_precision(time, *projected_norm);
      if (!uni20::isfinite(time)) throw std::invalid_argument("magnetization time must be finite");
      bethe::detail::CompensatedSum<R> sum(bethe::detail::precision_zero(time));
      for (auto const& pair : pairs)
      {
        if (!pair.term) throw std::logic_error("boundary magnetization has an unresolved pair");
        auto const& term = *pair.term;
        R const phase = term.frequency * time;
        if (!uni20::isfinite(phase)) throw std::overflow_error("magnetization phase overflow");
        sum.add(term.coefficient * (opposite_end ? term.opposite_end_parity : 1) * uni20::math::cos(phase));
      }
      if (!uni20::isfinite(sum.value())) throw std::overflow_error("magnetization sum overflow");
      return sum.value();
    }
};

namespace boundary_detail
{
inline std::size_t magnetization_pair_count(std::size_t count, std::size_t max_pairs)
{
  if (count == std::numeric_limits<std::size_t>::max()) throw std::length_error("boundary pair count overflow");
  auto const a = count % 2 ? count : count / 2, b = count % 2 ? (count + 1) / 2 : count + 1;
  if (!max_pairs || a > max_pairs / b) throw std::length_error("boundary magnetization pair budget exceeded");
  return a * b;
}

template <uni20::Real R, typename Evaluate>
uni20::async::AsyncTask magnetization_task(uni20::async::IAsyncScheduler& scheduler,
                                           uni20::async::ReadBuffer<BoundarySpectrum<R>> input, std::size_t pair_count,
                                           ExecutionOptions execution, Evaluate evaluate,
                                           uni20::async::WriteBuffer<BoundaryMagnetization<R>> output)
{
  auto const& spectrum = co_await input;
  auto destination = co_await output;
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary magnetization cancelled");
  BoundaryMagnetization<R> result(spectrum);
  auto const count = spectrum.records.size();
  result.states.resize(count);
  result.pairs.reserve(pair_count);
  for (std::size_t i = 0; i < count; ++i)
    for (std::size_t j = i; j < count; ++j)
      result.pairs.push_back({i + 1, j + 1, {}, {}, "state preparation unavailable"});
  if (spectrum.converged())
  {
    BoundarySpinFormFactors<R> calculator(spectrum.sites, spectrum.ground_state);
    auto completed = bethe::detail::continuous_work(
        scheduler, execution.max_in_flight, count,
        [&](std::size_t i, std::stop_token) { result.states[i] = calculator.prepare(*spectrum.records[i].state); },
        execution.stop);
    if (completed != count || execution.stop.stop_requested())
      throw ExecutionCancelled("boundary magnetization cancelled");
    if (std::all_of(result.states.begin(), result.states.end(), [](auto const& s) { return s.ready(); }))
    {
      R const zero = bethe::detail::precision_zero(spectrum.moments.weight);
      bethe::detail::CompensatedSum<R> norm(zero);
      for (auto const& s : result.states)
      {
        auto const a = *s.raising_amplitude();
        norm.add(2 * a * a);
      }
      result.projected_norm = norm.value();
      completed = bethe::detail::continuous_work(
          scheduler, execution.max_in_flight, pair_count,
          [&](std::size_t k, std::stop_token stop) {
            auto& pair = result.pairs[k];
            auto const i = pair.left_state - 1, j = pair.right_state - 1;
            pair.diagnostic = "signed boundary matrix element unresolved";
            pair.element = evaluate(calculator, result.states[i], result.states[j], stop);
            if (!pair.element->converged()) return;
            R const frequency = uni20::math::abs(spectrum.records[j].line->gap - spectrum.records[i].line->gap);
            R const coefficient = (i == j ? 2 : 4) * *result.states[i].raising_amplitude() *
                                  *result.states[j].raising_amplitude() * *pair.element->value;
            if (!uni20::isfinite(frequency) || !uni20::isfinite(coefficient))
            {
              pair.diagnostic = "nonfinite cosine term";
              return;
            }
            pair.term = BoundaryMagnetizationTerm<R>{
                frequency, coefficient, *result.states[i].reflection_parity() * *result.states[j].reflection_parity()};
            pair.diagnostic.clear();
          },
          execution.stop);
      if (completed != pair_count || execution.stop.stop_requested())
        throw ExecutionCancelled("boundary magnetization cancelled");
      R const gate = 4096 * uni20::scalar_like(zero, spectrum.sites) * uni20::numeric_limits<R>::epsilon(zero);
      result.matrix_ready =
          uni20::isfinite(*result.projected_norm) && *result.projected_norm >= zero &&
          *result.projected_norm <= uni20::scalar_like(zero, 1) + gate &&
          std::all_of(result.pairs.begin(), result.pairs.end(), [](auto const& p) { return p.term.has_value(); });
    }
  }
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary magnetization cancelled");
  destination = std::move(result);
}

template <uni20::Real R, typename Evaluate>
uni20::async::Async<BoundaryMagnetization<R>> submit_magnetization(std::size_t sites, SolverOptions<R> solver,
                                                                   std::size_t max_pairs, ExecutionOptions execution,
                                                                   Evaluate evaluate)
{
  (void)heisenberg::detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("boundary raised-singlet dynamics requires even N");
  auto const count = heisenberg::real_excitation_count(sites, uni20::half_int{1}, max_pairs);
  auto const pairs = magnetization_pair_count(count, max_pairs); // Before roots or quadratic allocation.
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto spectrum = boundary_real_triplet_spectrum_async<R>(sites, std::move(solver), count, execution);
  uni20::async::Async<BoundaryMagnetization<R>> output;
  auto task = magnetization_task(scheduler, spectrum.read(), pairs, execution, std::move(evaluate), output.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace boundary_detail

/// Owns the spectrum dependency and all admitted pair work. The scheduler is
/// borrowed until completion; dropping the result handle does not cancel.
template <uni20::Real R = double>
[[nodiscard]] uni20::async::Async<BoundaryMagnetization<R>>
boundary_real_triplet_magnetization_async(std::size_t sites, SolverOptions<R> solver = {},
                                          std::size_t max_pairs = 250000, ExecutionOptions execution = {})
{
  return boundary_detail::submit_magnetization(
      sites, std::move(solver), max_pairs, execution,
      [](auto const& calculator, auto const& x, auto const& y, std::stop_token) { return calculator(x, y); });
}
template <uni20::Real R = double>
[[nodiscard]] BoundaryMagnetization<R>
boundary_real_triplet_magnetization(std::size_t sites, SolverOptions<R> solver = {}, std::size_t max_pairs = 250000,
                                    ExecutionOptions execution = {})
{
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto pending = boundary_real_triplet_magnetization_async<R>(sites, std::move(solver), max_pairs, execution);
  return bethe::detail::take_result(pending, scheduler);
}
} // namespace bethe::experimental::xxx::open
