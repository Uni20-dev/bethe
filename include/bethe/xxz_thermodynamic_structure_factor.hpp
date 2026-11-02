// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/carlson_rf.hpp>
#include <bethe/xxz_dispersion.hpp>
#include <uni20/core/elliptic_theta_jet.hpp>

namespace bethe::xxz
{
enum class SpectralDensityStatus
{
  converged,
  outside_continuum,
  threshold,
  numerical_failure
};
template <uni20::Real Real> struct SpectralDensity
{
    std::optional<Real> value, error;
    std::size_t evaluations = 0;
    SpectralDensityStatus status = SpectralDensityStatus::numerical_failure;
    bool converged() const { return status != SpectralDensityStatus::numerical_failure; }
};

template <uni20::Real Real> struct TwoSpinonStructureFactorOptions
{
    Real tolerance = Real{4096} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_terms = 200000;
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct TwoSpinonStructureFactorOptions<uni20::mpreal>
{
    explicit TwoSpinonStructureFactorOptions(uni20::Precision p)
        : tolerance(4096 * uni20::numeric_limits<uni20::mpreal>::epsilon(p)), precision_(p)
    {
      if (p.is_exact()) throw std::invalid_argument("XXZ spectral kernel requires finite precision");
    }
    uni20::mpreal tolerance;
    std::size_t max_terms = 200000;
    uni20::Precision precision() const { return precision_; }

  private:
    uni20::Precision precision_;
};
#endif

template <uni20::Real Real> struct SpectralSheet
{
    int sign;
    Real lower, upper;
    bool present;
};

/// Infinite massive XXZ, zero field/temperature, H=J sum(SxSx+SySy+Delta SzSz).
/// Returns the transverse S^{-+}_2, with no broadening or sum-rule rescaling.
/// Caux--Mossel--Perez Castillo (2008), PDF Eqs. (44), (51), (56).
/// Physical momentum k=pi-Q; Neville theta_d is theta_3(z)/theta_3(0).
template <uni20::Real Real> class ThermodynamicTwoSpinonStructureFactor {
  public:
    explicit ThermodynamicTwoSpinonStructureFactor(Real delta, Real exchange = Real{1})
        : ThermodynamicTwoSpinonStructureFactor(delta, exchange, default_options(delta))
    {}
    ThermodynamicTwoSpinonStructureFactor(Real delta, Real exchange, TwoSpinonStructureFactorOptions<Real> options)
        : zero_(bethe::detail::precision_zero(delta)), pi_(bethe::detail::spinon_pi(delta)),
          band_(checked_delta(delta), exchange), options_(options), eta_(uni20::math::acosh(delta)),
          scale_(band_.maximum_energy()), kp_(band_.gap() / scale_), kappa_(zero_), critical_(zero_), K_(zero_),
          theta0_(uni20::elliptic_theta_constant(3, eta_ / pi_))
    {
      bethe::detail::require_same_precision(options.tolerance, zero_);
      if (options.precision() != uni20::precision_of(zero_))
        throw std::invalid_argument("XXZ spectral options must match the anisotropy precision");
      if (!uni20::isfinite(options.tolerance) || options.tolerance < Real{64} * eps())
        throw std::invalid_argument("XXZ log-kernel tolerance must be finite and >=64 epsilon");
      if (!(kp_ > Real{0} && kp_ < Real{1}))
        throw std::underflow_error("XXZ bandwidth or mass unresolved; use higher precision");
      if (!uni20::isfinite(Real{2} * scale_)) throw std::overflow_error("two-spinon energy exceeds scalar range");
      kappa_ = (Real{1} - kp_) / (Real{1} + kp_);
      critical_ = Real{2} * uni20::math::atan(uni20::math::sqrt(kp_));
      K_ = pi() / Real{2} * uni20::math::exp(Real{2} * theta0_.log_scale) * theta0_.scaled_value * theta0_.scaled_value;
    }

    Real gap() const { return band_.gap(); }
    Real maximum_energy() const { return scale_; }

    /// Separate sheets are essential: the energy support alone does not
    /// specify their unequal intensities. Signs follow the paper's C_sigma.
    std::array<SpectralSheet<Real>, 2> sheets(Real momentum) const
    {
      Real const k = folded(momentum);
      auto plus = [&](Real Q, int sign) -> SpectralSheet<Real> {
        if (Q < critical_) return {sign, zero_, zero_, false};
        Real const s = Q == Real{0} ? Real{0} : uni20::math::sin(Q / Real{2});
        Real const c = Q == pi() ? Real{0} : uni20::math::cos(Q / Real{2});
        Real const upper = Real{2} * scale_ * uni20::math::hypot(s, kp_ * c);
        Real const distance_from_pi = sign == 1 ? k : pi() - k;
        Real const lower = distance_from_pi >= critical_ ? scale_ * (Real{1} + kp_) * uni20::math::sin(k)
                                                         : Real{2} * scale_ * uni20::math::hypot(c, kp_ * s);
        return {sign, lower, upper, true};
      };
      return {plus(pi() - k, 1), plus(k, -1)};
    }

    SpectralDensity<Real> operator()(Real momentum, Real omega) const
    {
      bethe::detail::require_same_precision(omega, zero_);
      if (!uni20::isfinite(omega)) throw std::invalid_argument("frequency must be finite");
      Real const k = folded(momentum), w = omega / scale_;
      SpectralDensity<Real> out{zero_, zero_, 0, SpectralDensityStatus::outside_continuum};
      for (auto const& sheet : sheets(momentum))
      {
        if (!sheet.present || omega < sheet.lower || omega > sheet.upper) continue;
        if (omega == sheet.lower || omega == sheet.upper)
          return {{}, {}, out.evaluations, SpectralDensityStatus::threshold};
        // Rationalize both cancellations in Eq. (44)/(54). In particular
        // W_sigma and 1+sigma*cos(Q) vanish together at k=0,pi; compute their
        // quotient directly. This also includes the exact endpoint momenta.
        Real const sh = k == Real{0} ? Real{0} : uni20::math::sin(k / Real{2});
        Real const ch = k == pi() ? Real{0} : uni20::math::cos(k / Real{2});
        Real const a = Real{2} * (sheet.sign == 1 ? sh * sh : ch * ch);
        Real const other = Real{2} * (sheet.sign == 1 ? ch * ch : sh * sh);
        Real const s2 = a * other, w0 = (Real{1} + kp_) * uni20::math::sqrt(s2);
        Real const above = (w - w0) * (w + w0);
        Real const t2 = (w - kappa_ * w0) * (w + kappa_ * w0) * above;
        if (!(t2 > Real{0})) return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        Real const T = uni20::math::sqrt(t2), b = kappa_ * w0 * w0;
        Real const z2 = (above + Real{2} * kp_ * s2 + T) / (Real{2} * other * other);
        Real const signed_cos = (a - other) / Real{2};
        Real const A = b - signed_cos * w * w, B = b + signed_cos * w * w;
        Real const upper = sheet.upper / scale_;
        Real const upper_distance = (upper - w) * (upper + w);
        Real ratio2;
        if (B < Real{0})
        {
          Real const opposite = Real{2} * uni20::math::sqrt((a + kp_ * kp_ * other) / Real{2});
          ratio2 = other * other * upper_distance * (w - opposite) * (w + opposite) / ((A + T) * (T - B));
        }
        else
          ratio2 = (other / a) * upper_distance * (B + T) / (w * w * (A + T));
        if (!(ratio2 > Real{0})) return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        Real const numerator =
            sheet.sign == 1 ? kp_ * (w * w + b + T) : above + Real{2} * kp_ * w0 * w0 / (Real{1} + kp_) + T;
        Real const prefactor = numerator / (Real{2} * w * T * uni20::math::sqrt(ratio2));
        Real conditioning = Real{1} + w * w / upper_distance + w * w * w * w / t2;
        Real const complement = Real{1} - z2, distance = z2 - kp_ * kp_;
        if (!(complement > Real{0} && distance > Real{0} && prefactor > Real{0}))
          return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        Real const modulus2 = (Real{1} - kp_) * (Real{1} + kp_);
        Real const sin_amplitude = uni20::math::sqrt(complement / modulus2);
        Real const beta = sin_amplitude * bethe::detail::carlson_rf_unit(distance / modulus2, z2);
        auto const rate = log_rate(beta);
        out.evaluations += rate.evaluations;
        if (!rate.value) return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        Real const value = uni20::math::exp(*rate.value + uni20::math::log(prefactor) - uni20::math::log(scale_));
        conditioning += Real{1} / complement + Real{1} / distance;
        Real const error = value * (*rate.error + Real{128} * eps() * conditioning);
        if (!uni20::isfinite(value) || !uni20::isfinite(error) || value < positive_minimum())
          return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        *out.value += value;
        *out.error += error;
        if (!uni20::isfinite(*out.value) || !uni20::isfinite(*out.error))
          return {{}, {}, out.evaluations, SpectralDensityStatus::numerical_failure};
        out.status = SpectralDensityStatus::converged;
      }
      return out;
    }

  private:
    static Real checked_delta(Real delta)
    {
      if (!uni20::isfinite(delta) || delta <= Real{1})
        throw std::invalid_argument("massive XXZ structure factor requires finite Delta>1");
      return delta;
    }
    static TwoSpinonStructureFactorOptions<Real> default_options(Real const& delta)
    {
      if constexpr (uni20::has_runtime_precision_v<Real>)
        return TwoSpinonStructureFactorOptions<Real>(uni20::precision_of(delta));
      else
        return {};
    }
    Real pi() const { return pi_; }
    Real eps() const { return uni20::numeric_limits<Real>::epsilon(zero_); }
    Real positive_minimum() const
    {
      if constexpr (uni20::has_runtime_precision_v<Real>)
        return uni20::math::next_up(zero_);
      else
        return uni20::numeric_limits<Real>::min();
    }
    Real folded(Real k) const
    {
      bethe::detail::require_same_precision(k, zero_);
      if (!uni20::isfinite(k) || k < Real{0} || k > Real{2} * pi())
        throw std::invalid_argument("momentum must be finite and in [0,2*pi]");
      return std::min(k, Real{2} * pi() - k);
    }
    SpectralDensity<Real> log_rate(Real beta) const
    {
      // Extract the nondecaying harmonic series exactly:
      // -2 sum cos(n*x)/n = 2 log(2 |sin(x/2)|). The remainder is geometric.
      Real const x = Real{2} * pi() * beta / K_;
      Real const leading = Real{2} * uni20::math::log(Real{2} * uni20::math::abs(uni20::math::sin(x / Real{2})));
      bethe::detail::CompensatedSum<Real> sum(zero_);
      Real absolute = uni20::math::abs(leading), tail = zero_;
      std::size_t count = 0;
      bool complete = false;
      Real const decay = uni20::math::exp(-Real{2} * eta_);
      for (std::size_t index = 0; index < options_.max_terms; ++index)
      {
        std::size_t const n = index + 1;
        Real const t = uni20::math::exp(-Real{2} * Real(n) * eta_);
        Real const term = Real{2} * t * ((-Real{1} + Real{2} * t + t * t) * uni20::math::cos(Real(n) * x) - Real{2}) /
                          (Real(n) * (Real{1} + t) * (Real{1} - t * t));
        sum.add(term);
        absolute += uni20::math::abs(term);
        count = n;
        Real const next = t * decay;
        // |(-1+2t+t^2)cos(nx)-2| <= 6 for 0<=t<1.
        tail =
            Real{12} * next / ((Real(n) + Real{1}) * (Real{1} - next * next) * (-uni20::math::expm1(-Real{2} * eta_)));
        if (tail < options_.tolerance / Real{4})
        {
          complete = true;
          break;
        }
      }
      if (!complete) return {{}, {}, count, SpectralDensityStatus::numerical_failure};
      auto const theta = uni20::elliptic_theta_jet(3, uni20::complex<Real>(beta / (Real{2} * K_), zero_), eta_ / pi());
      Real const log_d =
          theta.log_scale - theta0_.log_scale + uni20::math::log(theta.derivative[0].real() / theta0_.scaled_value);
      Real const result = leading - sum.value() - Real{2} * log_d;
      Real const error = tail + Real{64} * eps() * (Real{1} + absolute + Real{2} * uni20::math::abs(log_d));
      return {result, error, count, SpectralDensityStatus::converged};
    }
    Real zero_, pi_;
    SpinonDispersion<Real> band_;
    TwoSpinonStructureFactorOptions<Real> options_;
    Real eta_, scale_, kp_, kappa_, critical_, K_;
    uni20::theta_constant_result<Real> theta0_;
};
} // namespace bethe::xxz
