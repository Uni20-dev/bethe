// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/experimental/detail/scaled_complex_product.hpp>
#include <bethe/solver.hpp>
#include <bethe/spinon.hpp>
#include <span>
#include <vector>

namespace bethe::experimental::xxz
{
struct VertexOptions
{
    // Zero selects scalar-dependent counts. Explicit counts are useful for
    // independent product/residue refinements; neither is a certified bound.
    std::size_t products = 0, residues = 0;
};

template <uni20::Real R> struct VertexAmplitudes
{
    std::array<std::array<uni20::complex<R>, 4>, 2> value{};
    std::array<std::array<R, 4>, 2> arithmetic_indicator{};
    unsigned polarizations = 0;
};

/// Original finite-q two/four-spinon vertex amplitudes, including both vacua.
/// xi=i exp(-i phi), u=-xi^2, q=-exp(-acosh Delta). One contour for n=4.
/// The m=0/1 starting residues are literal q-products; later levels use the
/// independently derived cancelled recurrence (limiting ratio q^6).
template <uni20::Real R> class VertexEvaluator {
    using C = uni20::complex<R>;
    using S = bethe::experimental::detail::ScaledComplexProduct<R>;

  public:
    explicit VertexEvaluator(R delta, VertexOptions options = {})
        : zero_(bethe::detail::precision_zero(delta)), pi_(bethe::detail::spinon_pi(delta))
    {
      if (!uni20::isfinite(delta) || delta <= R{1}) throw std::invalid_argument("vertex requires finite Delta>1");
      R const eta = uni20::math::acosh(delta);
      q_ = -uni20::math::exp(-eta);
      q2_ = q_ * q_;
      r_ = q2_ * q2_;
      if (!(r_ > R{0} && r_ < R{1})) throw std::overflow_error("vertex nome unresolved in this scalar");
      R const digits = -uni20::math::log(uni20::numeric_limits<R>::epsilon(zero_));
      auto count = [&](std::size_t requested, unsigned divisor, unsigned extra) {
        if (requested) return requested;
        R const estimate = uni20::math::ceil(digits / (R(divisor) * eta)) + R(extra);
        // Check in Real before converting: sufficiently small eta can put
        // the automatic count outside even size_t's representable range.
        if (!uni20::isfinite(estimate) || estimate > R{20000})
          throw std::invalid_argument("automatic vertex count exceeds 20000; use an asymptotic evaluator");
        // This is an integer in [0,20000], exactly representable in binary64.
        // Only the bounded work count crosses formats, never amplitude data.
        return std::size_t(static_cast<double>(estimate));
      };
      options_.products = count(options.products, 2, 8);
      options_.residues = count(options.residues, 6, 4);
      if (options_.products < 4 || options_.products > 20000 || options_.residues < 1 || options_.residues > 20000)
        throw std::invalid_argument("vertex requires products 4..20000 and residues 1..20000");
      for (R base : {r_, q2_, r_ * r_})
      {
        std::vector<R> powers(options_.products, R{1});
        for (std::size_t h = 1; h < powers.size(); ++h)
          powers[h] = powers[h - 1] * base;
        powers_.push_back(std::move(powers));
      }
      qp_r_ = qp(C{r_});
      rho_ = dp(C{r_});
      rho_ /= dp(C{r_ * q2_});
      theta_constant_ = qp(C{r_ * r_}, 2);
    }
    VertexOptions options() const { return options_; }
    R nome() const { return q_; }
    uni20::Precision precision() const { return uni20::precision_of(zero_); }

    VertexAmplitudes<R> operator()(std::span<R const> angles) const
    {
      std::size_t const n = angles.size();
      if (n != 2 && n != 4) throw std::invalid_argument("vertex requires two or four angles");
      std::array<C, 4> xi{}, u{};
      C product{uni20::scalar_like(zero_, 1), zero_};
      R separation = R{1} - q2_;
      for (std::size_t j = 0; j < n; ++j)
      {
        bethe::detail::require_same_precision(angles[j], zero_);
        if (!uni20::isfinite(angles[j]) || uni20::math::abs(angles[j]) > pi_ / R{2})
          throw std::invalid_argument("spectral angles must be finite in [-pi/2,pi/2]");
        xi[j] = C{uni20::math::sin(angles[j]), uni20::math::cos(angles[j])};
        u[j] = -xi[j] * xi[j];
        product *= u[j];
        for (std::size_t k = 0; k < j; ++k)
          separation = std::min(separation, uni20::math::abs(u[j] / u[k] - R{1}));
      }
      R const eps = uni20::numeric_limits<R>::epsilon(zero_);
      if (separation < R{64} * eps)
        throw std::domain_error("coincident/unresolved spectral nodes require a confluent formula or higher precision");
      S common(C{-q_});
      if (n == 4) common *= C{R{1} - R{1} / q2_};
      common *= qp(C{q2_});
      for (std::size_t j = 1; j < n; ++j)
        common *= qp_r_;
      for (std::size_t j = 0; j < n; ++j)
      {
        common *= rho_;
        common /= qp(-q_ * u[j]);
        common /= qp(-q_ * q2_ / u[j]);
        for (std::size_t k = 0; k < j; ++k)
        {
          C const a = u[j] / u[k];
          common *= dp(r_ * a);
          common *= dp(R{1} / a);
          common /= dp(r_ * q2_ * a);
          common /= dp(q2_ / a);
        }
      }
      VertexAmplitudes<R> result;
      for (auto& vacuum : result.value)
        vacuum.fill(C{zero_, zero_});
      for (auto& vacuum : result.arithmetic_indicator)
        vacuum.fill(zero_);
      result.polarizations = n == 2 ? 1 : 4;
      R const screen =
          R{128} * eps * (R(options_.products) * R(options_.products) + R(options_.residues) + R{1}) / separation;
      // Starting q-products do not depend on the polarization ell. Cache
      // them for the two possible starting pole levels, sharing both vacua.
      std::array<std::array<S, 2>, 4> starting;
      if (n == 4)
        for (unsigned j = 0; j < 4; ++j)
          for (unsigned start = 0; start < 2; ++start)
          {
            if (start >= options_.residues) continue;
            C const v = u[j] * (start ? r_ : R{1});
            S residue = common;
            residue *= C{start ? -r_ / (R{1} - r_) : R{1}};
            residue /= qp_r_;
            for (unsigned k = 0; k < 4; ++k)
            {
              if (k != j) residue /= qp(u[k] / v);
              residue /= qp(v / (q2_ * u[k]));
            }
            residue *= qp(-v / q_, 1);
            residue *= qp(-q_ * q2_ / v, 1);
            starting[j][start] = residue;
          }
      for (unsigned vacuum = 0; vacuum < 2; ++vacuum)
      {
        S base = common;
        if (vacuum) base /= C{n == 2 ? -q_ : q2_};
        std::array<std::array<S, 2>, 4> vacuum_start = starting;
        if (n == 4)
          for (unsigned j = 0; j < 4; ++j)
            for (unsigned start = 0; start < 2; ++start)
            {
              if (start >= options_.residues) continue;
              C const v = u[j] * (start ? r_ : R{1});
              if (vacuum)
              {
                vacuum_start[j][start] /= C{q2_};
                vacuum_start[j][start] *= -v;
              }
              vacuum_start[j][start] *= theta(-(vacuum ? q2_ : R{1} / q2_) * v * v / product);
            }
        for (unsigned ell = 0; ell < result.polarizations; ++ell)
        {
          S monomial = n == 2 ? base : S{};
          for (std::size_t j = 0; j < n; ++j)
          {
            int const power = int(j) + 1 - int(vacuum) - (n == 4 && j == ell ? 1 : 0);
            for (int h = 0; h < std::abs(power); ++h)
              monomial *= power < 0 ? R{1} / xi[j] : xi[j];
          }
          if (n == 2)
          {
            monomial *= theta(-(vacuum ? r_ : R{1}) / product);
            result.value[vacuum][ell] = monomial.value();
            result.arithmetic_indicator[vacuum][ell] = screen * uni20::math::abs(result.value[vacuum][ell]);
            if (!uni20::isfinite(result.arithmetic_indicator[vacuum][ell]))
              throw std::overflow_error("nonfinite vertex arithmetic indicator");
            continue;
          }
          bethe::detail::CompensatedSum<C> sum(C{zero_, zero_});
          bethe::detail::CompensatedSum<R> absolute(zero_);
          for (unsigned j = 0; j < 4; ++j)
          {
            unsigned const start = j > ell ? 1 : 0;
            if (start >= options_.residues) continue;
            C const v = u[j] * (start ? r_ : R{1});
            S residue = vacuum_start[j][start];
            residue *= monomial;
            for (unsigned k = 0; k < ell; ++k)
              residue *= q_ / v - R{1} / (q_ * u[k]);
            for (unsigned k = ell + 1; k < 4; ++k)
              residue *= R{1} / v - R{1} / u[k];
            C term = residue.value();
            sum.add(term);
            absolute.add(uni20::math::abs(term));
            C level_v = v;
            R level_r = start ? r_ : R{1};
            for (std::size_t m = start; m + 1 < options_.residues; ++m)
            {
              C ratio = -q2_ * r_ * product / (u[j] * (R{1} - level_r * r_));
              for (unsigned k = 0; k < 4; ++k)
              {
                ratio *= R{1} - level_v / (q2_ * u[k]);
                if (k != j) ratio /= r_ * level_v - u[k];
                if (k != ell)
                {
                  C const divisor = (k < ell ? q2_ : R{1}) * u[k];
                  ratio *= (R{1} - r_ * level_v / divisor) / (R{1} - level_v / divisor);
                }
              }
              term *= ratio;
              sum.add(term);
              absolute.add(uni20::math::abs(term));
              level_v *= r_;
              level_r *= r_;
            }
          }
          result.value[vacuum][ell] = sum.value();
          result.arithmetic_indicator[vacuum][ell] = screen * absolute.value();
          if (!uni20::isfinite(uni20::math::abs(sum.value())) || !uni20::isfinite(screen * absolute.value()))
            throw std::overflow_error("nonfinite vertex amplitude or arithmetic indicator");
        }
      }
      return result;
    }

  private:
    S qp(C z, unsigned base = 0) const
    {
      S result;
      for (R p : powers_[base])
        result *= R{1} - z * p;
      return result;
    }
    S dp(C z) const
    {
      // Reverse cumulative products implement powers (h+1) in O(N) work.
      S result, cumulative;
      for (auto h = powers_[0].size(); h > 0; --h)
      {
        cumulative *= R{1} - z * powers_[0][h - 1];
        result *= cumulative;
      }
      return result;
    }
    S theta(C z) const
    {
      S result = theta_constant_;
      result *= qp(z, 2);
      result *= qp(r_ * r_ / z, 2);
      return result;
    }
    VertexOptions options_;
    R zero_, pi_, q_, q2_, r_;
    std::vector<std::vector<R>> powers_;
    S qp_r_, rho_, theta_constant_;
};
} // namespace bethe::experimental::xxz
