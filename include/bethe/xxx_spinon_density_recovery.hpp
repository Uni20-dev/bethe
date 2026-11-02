// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxx_spinon_density.hpp>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <optional>
#include <set>
#include <stop_token>
#include <utility>

namespace bethe::heisenberg
{
struct SpinonRecoveryStage
{
    int bits, guard_bits;
    std::array<SpinonDensityOptions, 2> contours;
};

struct SpinonRecoveryPolicy
{
    std::string relative_tolerance = "1e-8", absolute_tolerance = "0";
    std::vector<SpinonRecoveryStage> stages;
    double max_wall_seconds = 600;
};

enum class SpinonRecoveryStatus
{
  accepted,
  exhausted,
  cancelled,
  time_limit,
  invalid_input,
  unsupported_limit,
  inconsistent_controls
};

template <uni20::Real R> struct CheckedSpinonDensity
{
    SpinonDensityEstimate<R> estimate;
    R precision_difference, scalar_precision_difference, arithmetic_screen, requested_error;
    SpinonDensityFailure failure;
    // Low-precision evidence imported exactly at the guard precision. Keeping
    // these three values lets restore rebuild the arithmetic screens.
    R low_value, low_log_scale, low_arithmetic_indicator;
    bool accepted() const { return failure == SpinonDensityFailure::none; }
};

template <uni20::Real R> struct SpinonRecoveryAttempt
{
    std::size_t stage, control;
    SpinonDensityFailure failure = SpinonDensityFailure::none;
    std::optional<CheckedSpinonDensity<R>> checked;
    std::string reason;
    double completed_wall_seconds = 0;
};

template <uni20::Real R> struct SpinonRecoveryResult
{
    SpinonRecoveryStatus status = SpinonRecoveryStatus::exhausted;
    std::vector<SpinonRecoveryAttempt<R>> attempts;
    std::vector<std::size_t> conflicting_attempts;
    std::optional<SpinonDensityEstimate<R>> estimate;
    std::optional<R> arithmetic_screen;
    // The last genuinely failed control, only when the declared stages are
    // exhausted. Cancellation, time limits and inconsistent controls have
    // their own terminal reason rather than an earlier attempt's failure.
    std::optional<SpinonDensityFailure> density_failure;
    std::string reason;
    double wall_seconds = 0;
    // Bind the declared numerical policy. Generic input factories remain the
    // caller's identity contract; decimal-input calls also retain exact text.
    SpinonRecoveryPolicy policy;
    std::vector<std::string> decimal_inputs;
    bool accepted() const { return status == SpinonRecoveryStatus::accepted; }
};

namespace detail
{
inline bool same_spinon_recovery_policy(SpinonRecoveryPolicy const& a, SpinonRecoveryPolicy const& b)
{
  if (a.relative_tolerance != b.relative_tolerance || a.absolute_tolerance != b.absolute_tolerance ||
      a.max_wall_seconds != b.max_wall_seconds || a.stages.size() != b.stages.size())
    return false;
  for (std::size_t i = 0; i < a.stages.size(); ++i)
  {
    auto const& x = a.stages[i];
    auto const& y = b.stages[i];
    if (x.bits != y.bits || x.guard_bits != y.guard_bits) return false;
    for (int j = 0; j < 2; ++j)
      if (x.contours[j].residue_terms != y.contours[j].residue_terms ||
          x.contours[j].tail_order != y.contours[j].tail_order)
        return false;
  }
  return true;
}

// Compare decimal identities before rounding; exact coincidences need a limit,
// while distinct inputs that collide at a low precision may be recoverable.
inline std::string spinon_decimal_key(std::string_view text)
{
  if (text.empty()) throw std::invalid_argument("empty decimal rapidity");
  std::size_t i = 0;
  bool negative = false, dot = false;
  if (text[i] == '+' || text[i] == '-') negative = text[i++] == '-';
  std::string digits;
  std::int64_t fractional = 0;
  for (; i < text.size() && text[i] != 'e' && text[i] != 'E'; ++i)
  {
    char c = text[i];
    if (c == '.' && !dot)
    {
      dot = true;
      continue;
    }
    if (c < '0' || c > '9') throw std::invalid_argument("rapidities must be finite decimal strings");
    digits += c;
    if (dot) ++fractional;
  }
  if (digits.empty()) throw std::invalid_argument("decimal rapidity has no digits");
  std::int64_t exponent = 0;
  if (i < text.size())
  {
    std::size_t start = ++i;
    if (i < text.size() && (text[i] == '+' || text[i] == '-')) ++i;
    if (i == text.size()) throw std::invalid_argument("decimal exponent has no digits");
    for (; i < text.size(); ++i)
      if (text[i] < '0' || text[i] > '9') throw std::invalid_argument("invalid decimal exponent");
    if (text[start] == '+') ++start;
    auto [end, error] = std::from_chars(text.data() + start, text.data() + text.size(), exponent);
    if (error != std::errc{} || end != text.data() + text.size())
      throw std::invalid_argument("decimal exponent overflow");
  }
  auto first = digits.find_first_not_of('0');
  if (first == std::string::npos) return "0";
  digits.erase(0, first);
  if (exponent < std::numeric_limits<std::int64_t>::min() + fractional)
    throw std::invalid_argument("decimal exponent overflow");
  exponent -= fractional;
  while (digits.back() == '0')
  {
    digits.pop_back();
    if (exponent == std::numeric_limits<std::int64_t>::max()) throw std::invalid_argument("decimal exponent overflow");
    ++exponent;
  }
  return (negative ? "-" : "") + digits + "e" + std::to_string(exponent);
}

template <uni20::Real R> XXXSpinonDensity<R> spinon_density_context(int bits)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return XXXSpinonDensity<R>(uni20::Precision::bits(bits));
  else
  {
    if (bits != uni20::numeric_limits<R>::digits(R{1}))
      throw std::invalid_argument("declared recovery bits differ from native scalar precision");
    return XXXSpinonDensity<R>{};
  }
}
} // namespace detail

/// Prepared bounded recovery. Contexts are immutable; attempts and scratch are
/// local to each call. The input factory must reconstruct the SAME point at
/// each supplied precision, never draw a new random sample. Two contour budgets
/// and two distinct arithmetic precisions are mandatory for acceptance.
/// Diagnostics are numerical screens, not rigorous error enclosures.
template <uni20::Real Low, uni20::Real High = Low> class XXXSpinonDensityRecovery {
  public:
    using Result = SpinonRecoveryResult<High>;

    explicit XXXSpinonDensityRecovery(SpinonRecoveryPolicy policy) : policy_(std::move(policy))
    {
      if (policy_.stages.empty() || policy_.stages.size() > 64 || !std::isfinite(policy_.max_wall_seconds) ||
          policy_.max_wall_seconds <= 0)
        throw std::invalid_argument("recovery requires 1..64 stages and a finite positive wall budget");
      SpinonRecoveryStage const* previous = nullptr;
      for (auto const& s : policy_.stages)
      {
        if (s.bits < 24 || s.guard_bits <= s.bits)
          throw std::invalid_argument("recovery requires distinct low/guard precisions");
        for (auto const& c : s.contours)
          if (c.residue_terms < 2 || c.tail_order < 1 || c.residue_terms > std::numeric_limits<int>::max() - 6 ||
              c.tail_order > std::numeric_limits<int>::max() - 4)
            throw std::invalid_argument("invalid recovery contour budget");
        auto nondecreasing = [](auto const& a, auto const& b) {
          return b.residue_terms >= a.residue_terms && b.tail_order >= a.tail_order;
        };
        auto same = [](auto const& a, auto const& b) {
          return a.residue_terms == b.residue_terms && a.tail_order == b.tail_order;
        };
        if (!nondecreasing(s.contours[0], s.contours[1]) || same(s.contours[0], s.contours[1]))
          throw std::invalid_argument("each recovery stage requires two increasing contour budgets");
        if (previous && (s.bits < previous->bits || s.guard_bits < previous->guard_bits ||
                         !nondecreasing(previous->contours[1], s.contours[0]) ||
                         (s.bits == previous->bits && s.guard_bits == previous->guard_bits &&
                          same(previous->contours[1], s.contours[0]))))
          throw std::invalid_argument("successive recovery stages must strengthen controls");
        contexts_.push_back(
            {detail::spinon_density_context<Low>(s.bits), detail::spinon_density_context<High>(s.guard_bits)});
        auto const& context = contexts_.back().high;
        auto relative = context.parse(policy_.relative_tolerance), absolute = context.parse(policy_.absolute_tolerance);
        if (!uni20::isfinite(relative) || !uni20::isfinite(absolute) || relative < High{0} || absolute < High{0} ||
            (relative == High{0} && absolute == High{0}))
          throw std::invalid_argument("require finite nonnegative nonzero recovery tolerances");
        previous = &s;
      }
    }

    SpinonRecoveryPolicy const& policy() const { return policy_; }

    /// Original decimal values are parsed afresh in every context. Invalid or
    /// exactly coincident inputs are terminal without running numerical kernels.
    Result operator()(std::span<std::string const> beta, std::stop_token stop = {}) const
    {
      return resume(beta, Result{}, [](Result const&) {}, stop);
    }

    /// Resume the same decimal point, publishing a const progress view after
    /// each completed low/guard control. The observer must copy anything it
    /// retains; exceptions propagate. Original spellings bind input identity.
    template <class Observer>
    Result resume(std::span<std::string const> beta, Result const& previous, Observer&& observer,
                  std::stop_token stop = {}) const
    {
      Result result;
      result.policy = policy_;
      result.decimal_inputs.assign(beta.begin(), beta.end());
      bool const fresh = previous.policy.stages.empty();
      if (fresh && (!previous.attempts.empty() || previous.wall_seconds != 0 || !previous.decimal_inputs.empty()))
        throw std::invalid_argument("unbound recovery history");
      if (!fresh && previous.decimal_inputs != result.decimal_inputs)
        throw std::invalid_argument("recovery decimal input identity changed");
      if (!fresh) result = restore(previous);
      if (beta.size() != 2 && beta.size() != 4 && beta.size() != 6)
      {
        result.status = SpinonRecoveryStatus::invalid_input;
        result.reason = "require two, four or six rapidities";
        return result;
      }
      std::set<std::string> keys;
      for (auto const& text : beta)
      {
        try
        {
          if (!keys.insert(detail::spinon_decimal_key(text)).second)
          {
            result.status = SpinonRecoveryStatus::unsupported_limit;
            result.reason = "exact rapidity coincidence requires a limiting calculation";
            return result;
          }
        }
        catch (std::invalid_argument const& error)
        {
          result.status = SpinonRecoveryStatus::invalid_input;
          result.reason = error.what();
          return result;
        }
      }
      return run(
          [&](auto const& context) {
            using R = decltype(context.parse("0"));
            std::vector<R> coordinates;
            for (auto const& text : beta)
            {
              try
              {
                coordinates.push_back(context.parse(text));
              }
              catch (std::out_of_range const& error)
              {
                // Syntax was validated above; this is a representability failure,
                // including native decimal parsing underflow, not an invalid point.
                throw SpinonDensityError(SpinonDensityFailure::exponent_range, error.what());
              }
              catch (std::overflow_error const& error)
              {
                throw SpinonDensityError(SpinonDensityFailure::exponent_range, error.what());
              }
              catch (std::underflow_error const& error)
              {
                throw SpinonDensityError(SpinonDensityFailure::exponent_range, error.what());
              }
            }
            return coordinates;
          },
          std::move(result), std::forward<Observer>(observer), stop);
    }

    /// Validate and reassess a saved control prefix without numerical kernels.
    /// Top-level status, estimate and acceptance flags are not trusted.
    /// File/producer identity and crash-time accounting belong to the caller.
    Result restore(Result const& previous) const
    {
      if (!detail::same_spinon_recovery_policy(previous.policy, policy_) || !std::isfinite(previous.wall_seconds) ||
          previous.wall_seconds < 0 || previous.attempts.size() > 2 * contexts_.size())
        throw std::invalid_argument("recovery policy, history size or wall budget changed");
      Result result;
      result.policy = policy_;
      result.decimal_inputs = previous.decimal_inputs;
      std::size_t expected_amplitudes = 0;
      if (!result.decimal_inputs.empty())
      {
        auto const n = result.decimal_inputs.size();
        if (n != 2 && n != 4 && n != 6) throw std::invalid_argument("invalid saved rapidity count");
        expected_amplitudes = n == 2 ? 1 : n == 4 ? 4 : 15;
        std::set<std::string> keys;
        for (auto const& text : result.decimal_inputs)
          if (!keys.insert(detail::spinon_decimal_key(text)).second)
            throw std::invalid_argument("saved rapidities require an unsupported coincidence limit");
      }
      result.wall_seconds = previous.wall_seconds;
      result.attempts.reserve(previous.attempts.size());
      double completed = 0;
      bool terminal = false;
      for (std::size_t i = 0; i < previous.attempts.size(); ++i)
      {
        auto const& saved = previous.attempts[i];
        if (terminal || saved.stage != i / 2 || saved.control != i % 2 ||
            !std::isfinite(saved.completed_wall_seconds) || saved.completed_wall_seconds < completed ||
            saved.completed_wall_seconds > result.wall_seconds)
          throw std::invalid_argument("recovery control order or consumed wall budget changed");
        completed = saved.completed_wall_seconds;
        validate_attempt(saved);
        if (expected_amplitudes && saved.checked && saved.checked->estimate.amplitudes.size() != expected_amplitudes)
          throw std::invalid_argument("saved polarization count differs from rapidity count");
        result.attempts.push_back(saved);
        if (i % 2 == 1) terminal = assess(result, contexts_[i / 2].high);
      }
      // An over-budget prefix cannot gain acceptance just by restarting it.
      if (result.wall_seconds >= policy_.max_wall_seconds)
      {
        result.status = SpinonRecoveryStatus::time_limit;
        result.estimate.reset();
        result.arithmetic_screen.reset();
        result.conflicting_attempts.clear();
        result.reason = "recovery wall budget exhausted";
      }
      else if (!terminal && result.attempts.size() == 2 * contexts_.size())
        exhaust(result);
      return result;
    }

    /// Cancellation/time are checked between kernels; they cannot interrupt a
    /// running special function. Unexpected exceptions propagate, not retries.
    template <class Source> Result recover(Source const& source, std::stop_token stop = {}) const
    {
      Result result;
      result.policy = policy_;
      return run(source, std::move(result), [](Result const&) {}, stop);
    }

    /// Factory variant: caller guarantees the SAME point at every precision
    /// and across restarts. Use resume(decimal_inputs,...) for bound text inputs.
    template <class Source, class Observer>
    Result recover(Source const& source, Result const& previous, Observer&& observer, std::stop_token stop = {}) const
    {
      Result result;
      if (previous.policy.stages.empty())
      {
        if (!previous.attempts.empty() || previous.wall_seconds != 0 || !previous.decimal_inputs.empty())
          throw std::invalid_argument("unbound recovery history");
        result.policy = policy_;
      }
      else
        result = restore(previous);
      return run(source, std::move(result), std::forward<Observer>(observer), stop);
    }

  private:
    template <class Source, class Observer>
    Result run(Source const& source, Result result, Observer&& observer, std::stop_token stop) const
    {
      if (result.accepted() || result.status == SpinonRecoveryStatus::inconsistent_controls ||
          result.status == SpinonRecoveryStatus::time_limit || result.attempts.size() == 2 * contexts_.size())
        return result;
      double const consumed = result.wall_seconds;
      auto const start = std::chrono::steady_clock::now();
      struct Interrupted
      {
          SpinonRecoveryStatus status;
      };
      auto check_budget = [&] {
        result.wall_seconds =
            consumed + std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (stop.stop_requested()) throw Interrupted{SpinonRecoveryStatus::cancelled};
        if (result.wall_seconds >= policy_.max_wall_seconds) throw Interrupted{SpinonRecoveryStatus::time_limit};
      };
      try
      {
        for (std::size_t index = result.attempts.size(); index < 2 * contexts_.size(); ++index)
        {
          std::size_t const stage = index / 2, control = index % 2;
          auto const& contexts = contexts_[stage];
          check_budget();
          SpinonRecoveryAttempt<High> attempt{.stage = stage, .control = control};
          try
          {
            auto low_beta = source(contexts.low);
            auto high_beta = source(contexts.high);
            for (auto const& x : low_beta)
              if (!uni20::isfinite(x))
                throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                                         "input outside low precision exponent range");
            for (auto const& x : high_beta)
              if (!uni20::isfinite(x))
                throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                                         "input outside guard precision exponent range");
            check_budget();
            auto low = contexts.low.evaluate(low_beta, policy_.stages[stage].contours[control]);
            check_budget();
            auto high = contexts.high.evaluate(high_beta, policy_.stages[stage].contours[control]);
            check_budget();
            attempt.checked = check(contexts.high, low, std::move(high));
            attempt.failure = attempt.checked->failure;
            attempt.reason = attempt.checked->accepted() ? "precision and contour screens passed"
                             : attempt.failure == SpinonDensityFailure::contour_tail ? "contour-tail screen failed"
                                                                                     : "arithmetic screen failed";
          }
          catch (SpinonDensityError const& error)
          {
            attempt.failure = error.failure();
            attempt.reason = error.what();
          }
          result.wall_seconds =
              consumed + std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
          attempt.completed_wall_seconds = result.wall_seconds;
          result.attempts.push_back(std::move(attempt));
          observer(std::as_const(result));
          check_budget();
          if (control == 1 && assess(result, contexts.high)) return result;
        }
        exhaust(result);
      }
      catch (Interrupted const& interrupted)
      {
        result.status = interrupted.status;
        result.reason = interrupted.status == SpinonRecoveryStatus::cancelled ? "recovery cancelled"
                                                                              : "recovery wall budget exhausted";
      }
      return result;
    }

    struct Contexts
    {
        XXXSpinonDensity<Low> low;
        XXXSpinonDensity<High> high;
    };
    SpinonRecoveryPolicy policy_;
    std::vector<Contexts> contexts_;

    CheckedSpinonDensity<High> check(XXXSpinonDensity<High> const& context, SpinonDensityEstimate<Low> const& low,
                                     SpinonDensityEstimate<High> high) const
    {
      return check_imported(context, context.import(low.value), context.import(low.log_scale),
                            context.import(low.arithmetic_indicator), std::move(high));
    }

    CheckedSpinonDensity<High> check_imported(XXXSpinonDensity<High> const& context, High low_value, High low_scale,
                                              High low_arithmetic, SpinonDensityEstimate<High> high) const
    {
      using uni20::math::abs;
      High difference = abs(high.value - low_value);
      High scalar_change = abs(uni20::math::expm1(low_scale - high.log_scale));
      High scalar = scalar_change * abs(high.value);
      if (scalar_change != High{0} && scalar == High{0})
        throw SpinonDensityError(SpinonDensityFailure::exponent_range, "precision diagnostic underflow");
      High arithmetic = low_arithmetic + high.arithmetic_indicator + difference + scalar;
      High tolerance =
          context.parse(policy_.absolute_tolerance) + context.parse(policy_.relative_tolerance) * abs(high.value);
      if (!uni20::isfinite(scalar) || !uni20::isfinite(arithmetic) || !uni20::isfinite(tolerance) ||
          arithmetic <= High{0})
        throw SpinonDensityError(SpinonDensityFailure::exponent_range, "unrepresentable recovery diagnostic");
      auto failure = high.contour_tail_error + arithmetic <= tolerance ? SpinonDensityFailure::none
                     : high.contour_tail_error > tolerance             ? SpinonDensityFailure::contour_tail
                                                                       : SpinonDensityFailure::arithmetic_precision;
      return {std::move(high),
              difference,
              scalar,
              arithmetic,
              tolerance,
              failure,
              std::move(low_value),
              std::move(low_scale),
              std::move(low_arithmetic)};
    }

    static void exhaust(Result& result)
    {
      result.reason = "all declared recovery stages exhausted";
      for (auto attempt = result.attempts.rbegin(); attempt != result.attempts.rend(); ++attempt)
        if (attempt->failure != SpinonDensityFailure::none)
        {
          result.density_failure = attempt->failure;
          result.reason += ": " + attempt->reason;
          break;
        }
    }

    void validate_attempt(SpinonRecoveryAttempt<High> const& attempt) const
    {
      auto failure = static_cast<int>(attempt.failure);
      if (failure < 0 || failure > static_cast<int>(SpinonDensityFailure::rounded_collision))
        throw std::invalid_argument("invalid recovery failure code");
      if (!attempt.checked)
      {
        if (attempt.failure == SpinonDensityFailure::none || attempt.failure == SpinonDensityFailure::contour_tail)
          throw std::invalid_argument("recovery control lacks its numerical evidence");
        return;
      }
      auto const& x = *attempt.checked;
      auto const& e = x.estimate;
      auto const& contexts = contexts_[attempt.stage];
      auto const& options = policy_.stages[attempt.stage].contours[attempt.control];
      auto const zero = contexts.high.parse("0");
      auto valid = [&](High const& value, bool nonnegative) {
        bethe::detail::require_same_precision(value, zero);
        if (!uni20::isfinite(value) || (nonnegative && value < zero))
          throw std::invalid_argument("invalid recovery numerical evidence");
      };
      for (auto const& value :
           {e.value, e.contour_tail_error, e.arithmetic_indicator, e.contour_norm, e.contour_norm_error,
            e.contour_arithmetic_indicator, x.precision_difference, x.scalar_precision_difference, x.arithmetic_screen,
            x.requested_error, x.low_value, x.low_arithmetic_indicator})
        valid(value, true);
      valid(e.log_scale, false);
      valid(x.low_log_scale, false);
      if (e.value == zero || e.contour_norm == zero || x.low_value == zero || e.terms != options.residue_terms ||
          e.order != options.tail_order ||
          (e.amplitudes.size() != 1 && e.amplitudes.size() != 4 && e.amplitudes.size() != 15))
        throw std::invalid_argument("recovery estimate dimensions or controls changed");
      for (auto const& amplitude : e.amplitudes)
      {
        valid(amplitude.real(), false);
        valid(amplitude.imag(), false);
      }
      for (auto const& value : {x.low_value, x.low_log_scale, x.low_arithmetic_indicator})
        if (contexts.high.import(contexts.low.import(value)) != value)
          throw std::invalid_argument("low-precision recovery evidence is not exactly representable");
      auto rebuilt = check_imported(contexts.high, x.low_value, x.low_log_scale, x.low_arithmetic_indicator, e);
      if (x.failure != rebuilt.failure || attempt.failure != rebuilt.failure ||
          x.precision_difference != rebuilt.precision_difference ||
          x.scalar_precision_difference != rebuilt.scalar_precision_difference ||
          x.arithmetic_screen != rebuilt.arithmetic_screen || x.requested_error != rebuilt.requested_error)
        throw std::invalid_argument("recovery diagnostics or acceptance changed");
    }

    bool assess(Result& result, XXXSpinonDensity<High> const& context) const
    {
      auto passed = [](auto const& attempt) { return attempt.checked && attempt.checked->accepted(); };
      auto const n = result.attempts.size();
      if (!passed(result.attempts[n - 2]) || !passed(result.attempts[n - 1])) return false;
      auto const& final = *result.attempts.back().checked;
      for (std::size_t i = 0; i + 1 < n; ++i)
        if (passed(result.attempts[i]))
        {
          auto const& earlier = *result.attempts[i].checked;
          auto difference = uni20::math::abs(final.estimate.value - context.import(earlier.estimate.value));
          auto allowance = final.estimate.contour_tail_error + final.arithmetic_screen +
                           context.import(earlier.estimate.contour_tail_error) +
                           context.import(earlier.arithmetic_screen);
          if (!uni20::isfinite(allowance) || difference > allowance) result.conflicting_attempts.push_back(i);
        }
      if (!result.conflicting_attempts.empty())
      {
        result.status = SpinonRecoveryStatus::inconsistent_controls;
        result.reason = "individually passing controls disagree";
      }
      else
      {
        result.status = SpinonRecoveryStatus::accepted;
        result.estimate = final.estimate;
        result.arithmetic_screen = final.arithmetic_screen;
        result.reason = "both contour budgets and all earlier passing controls agree";
      }
      return true;
    }
};
} // namespace bethe::heisenberg
