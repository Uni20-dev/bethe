// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/controlled_hypergeometric_tail.hpp>
#include <bethe/detail/xxx_barnes_transition_rate.hpp>
#include <bethe/solver.hpp>
#include <span>
#include <string>
#include <string_view>
#include <uni20/core/complex_log_gamma.hpp>
#include <uni20/core/scalar_io.hpp>
#include <vector>

namespace bethe::heisenberg
{
enum class SpinonDensityFailure
{
  none,
  contour_tail,
  arithmetic_precision,
  tail_domain,
  exponent_range,
  special_function_budget,
  rounded_collision
};

class SpinonDensityError : public std::runtime_error {
  public:
    SpinonDensityError(SpinonDensityFailure failure, std::string const& reason)
        : std::runtime_error(reason), failure_(failure)
    {}
    SpinonDensityFailure failure() const noexcept { return failure_; }

  private:
    SpinonDensityFailure failure_;
};

template <uni20::Real R> struct SpinonDensityEstimate
{
    R value, contour_tail_error, arithmetic_indicator, log_scale;
    R contour_norm, contour_norm_error, contour_arithmetic_indicator;
    std::vector<uni20::complex<R>> amplitudes;
    int terms, order;
};

/// Per-evaluation contour work. Arithmetic precision belongs to the context;
/// increasing these budgets does not repair arithmetic cancellation.
struct SpinonDensityOptions
{
    int residue_terms = 32;
    int tail_order = 12;
};

/// Physical longitudinal W_n against labelled product(dp), with conservation
/// deltas unintegrated and exactly one permutation factorial already included.
/// beta inputs mean p=2 atan(exp(beta)); the rapidity Jacobian is NOT included.
/// Supports n=2,4,6 at distinct finite rapidities. The returned diagnostics are
/// estimates, not an interval certificate or a claim of integration convergence.
/// One immutable context per precision; concurrent evaluations own their scratch.
template <uni20::Real R> class XXXSpinonDensity {
  public:
    using C = uni20::complex<R>;
    explicit XXXSpinonDensity(uni20::log_gamma_options options = {})
      requires(!uni20::has_runtime_precision_v<R>)
        : XXXSpinonDensity(uni20::native_precision<R>(), options)
    {}

    explicit XXXSpinonDensity(uni20::Precision precision, uni20::log_gamma_options options = {})
        : functions_(precision, options), gamma_(precision, options), rate_(precision, options)
    {
      auto g = functions_.log_abs_barnes_g_right(C(R{1} / R{2}, R{0}));
      log_b_ = -R{2} * g.value;
      log_b_error_ = R{2} * (g.truncation_bound + g.roundoff_indicator);
      zero_ = g.value - g.value;
      one_ = zero_ + R{1};
      using std::atan;
      pi_ = R{4} * atan(one_);
      eps_ = uni20::numeric_limits<R>::epsilon(one_);
      gamma_seed_ = this->log_gamma_shifted(C(-one_ / R{2}, zero_));
    }

    uni20::Precision precision() const { return uni20::precision_of(one_); }

    R parse(std::string_view text) const
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        return R(text, one_.precision());
      else
        return uni20::parse_real<R>(text);
    }
    template <class T> R import(T const& value) const
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        return R(value, one_.precision());
      else
        return R(value);
    }
    auto bits() const
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        return one_.precision().bit_count();
      else
        return uni20::numeric_limits<R>::digits(one_);
    }

    SpinonDensityEstimate<R> evaluate(std::span<R const> beta, SpinonDensityOptions options = {}) const
    try
    {
      int const terms = options.residue_terms, order = options.tail_order;
      using bethe::detail::SpinonSum;
      using uni20::math::abs;
      using uni20::math::exp;
      using uni20::math::log;
      if ((beta.size() != 2 && beta.size() != 4 && beta.size() != 6) || terms < 2 || order < 1)
        throw std::invalid_argument("require two, four or six distinct real rapidities, terms>=2 and order>=1");
      int const n = static_cast<int>(beta.size());
      if (terms > std::numeric_limits<int>::max() - 6 || order > std::numeric_limits<int>::max() - 4)
        throw std::invalid_argument("residue terms/tail order exceed index representation");
      for (int j = 0; j < n; ++j)
      {
        if (!uni20::isfinite(beta[j])) throw std::invalid_argument("rapidities must be finite");
        // Reject unset, exact-only or mixed-precision numerical inputs.
        bethe::detail::require_same_precision(beta[j], one_);
        for (int k = 0; k < j; ++k)
          if (beta[k] == beta[j])
            throw SpinonDensityError(SpinonDensityFailure::rounded_collision,
                                     "coincident rapidities require a limiting calculation");
      }
      R norm = one_, norm_error = zero_, norm_round = zero_;
      std::vector<C> amplitudes{C(one_, zero_)};
      if (n != 2)
      {
        SpinonSum<R> center_sum(zero_);
        for (auto const& x : beta)
          center_sum.add(x);
        R center = center_sum.value() / R(n);
        C cz(zero_, zero_), imaginary(zero_, one_);
        using Moments = std::array<C, 4>;
        using Errors = std::array<R, 4>;
        // Point-local immutable plan and exclusive scratch: all 72 six-spinon
        // tails reuse these without any mutable state in the precision context.
        bethe::detail::ControlledTailPlan<R> tail_plan(one_, order);
        bethe::detail::ControlledTailWorkspace<R> tail_workspace;
        std::vector<Moments> values(n, Moments{cz, cz, cz, cz});
        std::vector<Errors> errors(n, Errors{zero_, zero_, zero_, zero_}), rounds = errors;
        R const two_pi = R{2} * pi_, half_pi = pi_ / R{2};
        std::vector<R> positive_exp, negative_exp, hyperbolic, ladder_imaginary;
        for (int j = 0; j < n; ++j)
        {
          positive_exp.push_back(exp(beta[j]));
          negative_exp.push_back(exp(-beta[j]));
          hyperbolic.push_back((positive_exp.back() + negative_exp.back()) / R{2});
        }
        for (int m = 0; m <= terms; ++m)
          ladder_imaginary.push_back(half_pi - two_pi * R(m));
        std::vector<C> differences(n * n, cz);
        for (int j = 0; j < n; ++j)
          for (int k = 0; k < n; ++k)
            differences[j * n + k] = C(zero_, (beta[j] - beta[k]) / two_pi);
        // Residues shared by all polarizations; only the initial Gamma seed
        // needs special functions. The rest of each ladder is rational.
        std::vector<std::vector<C>> residues(n, std::vector<C>(terms + 1, cz));
        std::vector<R> seed_relative_errors(n, zero_);
        for (int j = 0; j < n; ++j)
        {
          auto [seed, seed_error] = gamma_seed_;
          std::vector<C> x;
          for (int k = 0; k < n; ++k)
            if (k != j)
            {
              x.push_back(-differences[j * n + k]);
              for (auto const& z : {x.back(), C(-one_ / R{2} - x.back())})
              {
                auto [value, error] = this->log_gamma_shifted(z);
                seed += value;
                seed_error += error;
              }
            }
          // |exp(delta)-1| <= |delta| exp(|delta|). The special-function
          // roundoff component remains an indicator, not an interval bound.
          seed_relative_errors[j] = seed_error * exp(seed_error);
          residues[j][0] = two_pi * imaginary * exp(seed);
          if (!uni20::isfinite(residues[j][0]) || residues[j][0] == cz)
            throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                                     "unrepresentable Gamma residue seed; increase precision or exponent range");
          for (int m = 0; m < terms; ++m)
          {
            C ratio(-(one_ * R(m) - one_ / R{2}) / R(m + 1), zero_);
            for (auto const& a : x)
              ratio *= (R(m) - one_ / R{2} - a) / (a - R(m + 1));
            residues[j][m + 1] = residues[j][m] * ratio;
            if (!uni20::isfinite(residues[j][m + 1]) || residues[j][m + 1] == cz)
              throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                                       "unrepresentable Gamma residue ladder; increase precision or exponent range");
          }
        }
        std::vector<C> samples;
        samples.reserve(terms + 1);
        std::vector<C> a, b;
        a.reserve(n + 1);
        b.reserve(n + 1);
        for (int ell = 0; ell < n; ++ell)
        {
          std::array<std::vector<C>, 2> sums;
          std::array<std::vector<R>, 2> bounds, arithmetic;
          for (int j = 0; j < n; ++j)
          {
            samples.clear();
            for (int m = 0; m <= terms; ++m)
            {
              C alpha(beta[j], ladder_imaginary[m]);
              C polynomial(one_, zero_);
              for (int k = 0; k < ell; ++k)
                polynomial *= alpha - beta[k] + imaginary * half_pi;
              for (int k = ell + 1; k < n; ++k)
                polynomial *= beta[k] - alpha + imaginary * half_pi;
              samples.push_back(m == 0 && j > ell ? cz : residues[j][m] * imaginary * hyperbolic[j] * polynomial);
            }
            for (int moment = 0; moment < (n == 4 ? 1 : 2); ++moment)
            {
              SpinonSum<C> prefix(cz);
              SpinonSum<R> magnitudes(zero_);
              C last = cz;
              for (int m = 0; m <= terms; ++m)
              {
                C v = samples[m];
                if (moment) v *= C(beta[j] - center, ladder_imaginary[m]);
                if (m < terms)
                {
                  prefix.add(v);
                  magnitudes.add(abs(v));
                }
                else
                  last = v;
              }
              a.clear();
              b.clear();
              for (int k = 0; k < n; ++k)
              {
                C const& delta = differences[j * n + k];
                a.push_back(delta + (k < ell ? one_ / R{2} : -one_ / R{2}));
                b.push_back(delta + (k > ell ? zero_ : one_));
              }
              if (moment)
              {
                C offset(-one_ / R{4}, (beta[j] - center) / two_pi);
                a.push_back(offset + one_);
                b.push_back(offset);
              }
              auto tail = bethe::detail::controlled_hypergeometric_tail<R>(a, b, n / 2 + 1 - moment, terms, tail_plan,
                                                                           tail_workspace);
              sums[moment].push_back(prefix.value() + last * tail.multiplier);
              bounds[moment].push_back(abs(last) * tail.error_multiplier);
              arithmetic[moment].push_back(abs(last) * tail.arithmetic_indicator +
                                           (R{64} * eps_ * R(terms + n) + seed_relative_errors[j]) *
                                               (magnitudes.value() + abs(last * tail.multiplier)));
            }
          }
          for (int moment = 0; moment < (n == 4 ? 1 : 2); ++moment)
            for (int sign : {1, -1})
            {
              if (n == 4 && sign == -1) continue;
              SpinonSum<C> value(cz);
              SpinonSum<R> error(zero_), round(zero_);
              for (int j = 0; j < n; ++j)
              {
                C factor =
                    n == 4 ? C(one_, zero_) : R(sign) * imaginary * (sign == 1 ? positive_exp[j] : negative_exp[j]);
                value.add(factor * sums[moment][j]);
                error.add(abs(factor) * bounds[moment][j]);
                round.add(abs(factor) * arithmetic[moment][j]);
              }
              int index = moment + (sign == 1 ? 0 : 2);
              values[ell][index] = value.value();
              errors[ell][index] = error.value();
              rounds[ell][index] = round.value();
            }
        }
        amplitudes.clear();
        SpinonSum<R> ns(zero_), es(zero_), rs(zero_);
        auto add = [&](C const& v, R const& e, R const& r) {
          R m = abs(v);
          amplitudes.push_back(v);
          ns.add(m * m);
          es.add(R{2} * m * e + e * e);
          rs.add(R{2} * (m + e) * r + r * r);
        };
        if (n == 4)
          for (int ell = 0; ell < n; ++ell)
            add(values[ell][0], errors[ell][0], rounds[ell][0]);
        else
          for (int ell = 0; ell < n; ++ell)
            for (int other = ell + 1; other < n; ++other)
            {
              auto const& l = values[ell];
              auto const& r = values[other];
              C amplitude = ((l[1] + pi_ * imaginary * l[0]) * r[2] - l[0] * r[3] -
                             (l[3] + pi_ * imaginary * l[2]) * r[0] + l[2] * r[1]) /
                            R{2};
              SpinonSum<R> tail(zero_), round(zero_);
              for (auto ij : {std::pair{1, 2}, {0, 3}, {0, 2}, {3, 0}, {2, 1}, {2, 0}})
              {
                auto [i, j] = ij;
                R scale = ((i == 0 && j == 2) || (i == 2 && j == 0)) ? pi_ : one_;
                R x = abs(l[i]), y = abs(r[j]), ex = errors[ell][i], ey = errors[other][j];
                R rx = rounds[ell][i], ry = rounds[other][j];
                tail.add(scale * (x * ey + y * ex + ex * ey));
                round.add(scale * (x * ry + y * rx + rx * ry + ex * ry + rx * ey + R{192} * eps_ * x * y));
              }
              add(amplitude, tail.value() / R{2}, round.value() / R{2});
            }
        norm = ns.value();
        norm_error = es.value();
        norm_round = rs.value();
      }
      R coefficient = -log(R{4} * one_);
      if (n == 4) coefficient = -R{15} * log(R{2} * one_) - log(R{24} * one_) - R{14} * log(pi_) - R{8} * log_b_;
      if (n == 6) coefficient = -R{41} * log(R{2} * one_) - log(R{720} * one_) - R{42} * log(pi_) - R{24} * log_b_;
      SpinonSum<R> pairs(zero_);
      SpinonSum<R> pair_errors(zero_);
      for (int j = 0; j < n; ++j)
        for (int k = j + 1; k < n; ++k)
        {
          auto rate = rate_((beta[j] - beta[k]) / (R{2} * pi_));
          pairs.add(rate.value);
          pair_errors.add(rate.truncation_bound + rate.roundoff_indicator);
        }
      R log_scale = coefficient + pairs.value(), scale = exp(log_scale), value = scale * norm;
      R log_scale_error = pair_errors.value() + R(n == 2 ? 0 : n == 4 ? 8 : 24) * log_b_error_;
      R arithmetic = scale * norm_round + (R{64} * eps_ * (one_ + R(n * (n - 1) / 2) + abs(log_scale)) +
                                           log_scale_error * exp(log_scale_error)) *
                                              abs(value);
      R tail = scale * norm_error;
      for (auto const& x : {value, tail, arithmetic, norm, norm_error, norm_round})
        if (!uni20::isfinite(x) || x < zero_)
          throw SpinonDensityError(SpinonDensityFailure::arithmetic_precision,
                                   "nonfinite or negative density/diagnostic");
      // Distinct finite rapidities have a positive pair scale. Never accept
      // an underflowed zero as a successful physical zero.
      if (scale == zero_ || norm == zero_ || value == zero_)
        throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                                 "density underflow or unresolved zero contour norm; increase precision");
      return {value, tail, arithmetic, log_scale, norm, norm_error, norm_round, amplitudes, terms, order};
    }

    catch (bethe::detail::TailDomainError const& e)
    {
      throw SpinonDensityError(SpinonDensityFailure::tail_domain, e.what());
    }
    catch (bethe::detail::SpinonFunctionError const& e)
    {
      throw SpinonDensityError(e.failure() == bethe::detail::SpinonFunctionFailure::budget
                                   ? SpinonDensityFailure::special_function_budget
                                   : SpinonDensityFailure::arithmetic_precision,
                               e.what());
    }
    catch (SpinonDensityError const&)
    {
      throw;
    }
    catch (std::overflow_error const& e)
    {
      throw SpinonDensityError(SpinonDensityFailure::exponent_range, e.what());
    }
    catch (std::underflow_error const& e)
    {
      throw SpinonDensityError(SpinonDensityFailure::exponent_range, e.what());
    }
    catch (std::runtime_error const& e)
    {
      // Uni20's log-Gamma context reports exhausted recurrence budgets as
      // runtime_error. Invalid arguments and allocation failures still propagate.
      throw SpinonDensityError(SpinonDensityFailure::special_function_budget, e.what());
    }

  private:
    bethe::detail::SpinonBarnesContext<R> functions_;
    uni20::log_gamma_context<R> gamma_;
    heisenberg::detail::XXXBarnesTransitionRate<R> rate_;
    R log_b_, log_b_error_, zero_, one_, pi_, eps_;
    std::pair<C, R> gamma_seed_;
    std::pair<C, R> log_gamma_shifted(C const& z) const
    {
      // All seed arguments have Re(z)=0 or -1/2; one recurrence shift
      // enters the supported right half plane. Only exp of this log is used.
      auto shifted = gamma_.log_gamma_right(z + one_);
      C logarithm = uni20::math::log(z);
      R round = R{64} * eps_ * (one_ + uni20::math::abs(shifted.value) + uni20::math::abs(logarithm));
      return {shifted.value - logarithm, shifted.truncation_bound + shifted.roundoff_indicator + round};
    }
};

} // namespace bethe::heisenberg
