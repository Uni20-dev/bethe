// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/experimental/detail/xxz_spectral_kinematics.hpp>
#include <bethe/experimental/xxz_vertex.hpp>
#include <limits>

namespace bethe::experimental::xxz
{
template <uni20::Real R> struct VertexMoments
{
    R weight{}, first_moment{}, arithmetic_indicator{}, first_moment_arithmetic_indicator{};
    std::size_t nodes = 0;
};

template <uni20::Real R> struct VertexCut
{
    std::vector<R> density, arithmetic_indicator;
    VertexMoments<R> moments;
    std::array<R, 2> sector_weights{};
};

/// Fixed deterministic quadrature of original finite-q amplitudes. No fitted
/// normalization and no claim of quadrature convergence from a single order.
/// R is used throughout geometry, products, amplitudes and accumulation.
template <uni20::Real R> class ThermodynamicVertexSpectrum {
  public:
    explicit ThermodynamicVertexSpectrum(R delta, R exchange = R{1}, VertexOptions options = {})
        : kinematics_(delta, exchange), vertex_(delta, options)
    {}
    SpectralKinematics<R> const& kinematics() const { return kinematics_; }
    VertexOptions options() const { return vertex_.options(); }

    VertexMoments<R> integrated(unsigned n, std::size_t order, R offset = R{0}) const
    {
      offset = checked_offset(offset);
      validate(n, order, offset, false);
      validate_energy_range(n);
      std::size_t count = 1;
      for (unsigned j = 0; j < n; ++j)
      {
        if (count > std::numeric_limits<std::size_t>::max() / order)
          throw std::length_error("vertex grid node count overflows size_t");
        count *= order;
      }
      std::array<std::vector<R>, 4> angle, energy;
      for (unsigned j = 0; j < n; ++j)
        for (std::size_t i = 0; i < order; ++i)
        {
          R const a = node_angle(i, j, order, offset);
          angle[j].push_back(a);
          energy[j].push_back(kinematics_.energy(kinematics_.momentum(a)));
        }
      std::array<bethe::detail::CompensatedSum<R>, 4> sums;
      sums.fill(bethe::detail::CompensatedSum<R>(kinematics_.zero()));
      for (std::size_t index = 0; index < count; ++index)
      {
        std::array<R, 4> angles{};
        R e = kinematics_.zero();
        auto rest = index;
        for (unsigned j = 0; j < n; ++j)
        {
          auto const i = rest % order;
          rest /= order;
          angles[j] = angle[j][i];
          e += energy[j][i];
        }
        auto const a = vertex_(std::span<R const>{angles.data(), n});
        R weight = kinematics_.zero(), error = kinematics_.zero();
        for (unsigned v = 0; v < 2; ++v)
          for (unsigned l = 0; l < a.polarizations; ++l)
          {
            R const absolute = uni20::math::abs(a.value[v][l]), screen = a.arithmetic_indicator[v][l];
            weight += absolute * absolute;
            error += R{2} * absolute * screen + screen * screen;
          }
        weight /= R(2 * factorial(n));
        error /= R(2 * factorial(n));
        sums[0].add(weight);
        sums[1].add(weight * e);
        sums[2].add(error);
        sums[3].add(error * e);
      }
      for (auto const& sum : sums)
        require_finite(sum.value());
      return {sums[0].value() / R(count), sums[1].value() / R(count), sums[2].value() / R(count),
              sums[3].value() / R(count), count};
    }

    /// Executor(count, body) must invoke body(row) once for each row and join
    /// before returning. Frontends can pass a Uni20 scheduler; the library
    /// acquires no threading dependency. Reduction is in fixed row order.
    template <typename Execute>
    VertexCut<R> broadened_cut(unsigned n, R momentum, std::span<R const> omega, R width, std::size_t order,
                               Execute&& execute, R offset = R{0}) const
    {
      offset = checked_offset(offset);
      bethe::detail::require_same_precision(momentum, kinematics_.zero());
      bethe::detail::require_same_precision(width, kinematics_.zero());
      for (auto const& w : omega)
        bethe::detail::require_same_precision(w, kinematics_.zero());
      validate(n, order, offset, true);
      validate_energy_range(n);
      if (!uni20::isfinite(momentum) || momentum < R{0} || momentum > R{2} * pi() || !uni20::isfinite(width) ||
          width <= R{0})
        throw std::invalid_argument("invalid momentum or Gaussian width");
      // The pair quadrature explicitly resolves Q=pi-k. A nonzero k lost
      // in that subtraction would silently become an endpoint calculation.
      R const Q = pi() - momentum;
      if ((momentum != R{0} && Q == pi()) || (momentum != R{2} * pi() && Q == -pi()))
        throw std::underflow_error("physical momentum shift unresolved; increase precision");
      R const gaussian_normalization = uni20::math::sqrt(R{2} * pi()) / width;
      if (!uni20::isfinite(gaussian_normalization))
        throw std::overflow_error("Gaussian normalization outside scalar range; increase width or precision");
      for (R w : omega)
        if (!uni20::isfinite(w)) throw std::invalid_argument("nonfinite frequency");
      auto const rows = n == 4 ? order : std::size_t{1};
      if (omega.size() > 1000000 / rows) throw std::length_error("spectral scratch grid exceeds one million values");
      auto const rule = bethe::detail::gauss_legendre<R>(order, kinematics_.precision());
      std::array<std::vector<R>, 2> free_phi, free_p, free_e;
      if (n == 4)
        for (unsigned j = 0; j < 2; ++j)
          for (std::size_t i = 0; i < order; ++i)
          {
            R const phi = node_angle(i, j, order, offset), p = kinematics_.momentum(phi);
            free_phi[j].push_back(phi);
            free_p[j].push_back(p);
            free_e[j].push_back(kinematics_.energy(p));
          }
      std::vector<VertexCut<R>> partial(rows);
      execute(rows, [&](std::size_t row) {
        auto& out = partial.at(row);
        std::vector<bethe::detail::CompensatedSum<R>> density(omega.size(),
                                                              bethe::detail::CompensatedSum<R>(kinematics_.zero())),
            error(omega.size(), bethe::detail::CompensatedSum<R>(kinematics_.zero()));
        std::array<bethe::detail::CompensatedSum<R>, 6> moments;
        moments.fill(bethe::detail::CompensatedSum<R>(kinematics_.zero()));
        for (std::size_t col = 0; col < rows; ++col)
        {
          R const pfree = n == 4 ? free_p[0][row] + free_p[1][col] : kinematics_.zero();
          R const efree = n == 4 ? free_e[0][row] + free_e[1][col] : kinematics_.zero();
          for (unsigned sector = 0; sector < 2; ++sector)
          {
            int const sign = sector ? 1 : -1;
            R const raw = Q + R(sector) * pi() - pfree;
            R const pair = raw - R{2} * pi() * uni20::math::floor(raw / (R{2} * pi()));
            R const low = std::max(kinematics_.zero(), pair - pi()), high = std::min(pi(), pair);
            if (low == high) continue; // genuine zero-length integration interval
            R const phi_low = kinematics_.angle(low), phi_high = kinematics_.angle(high);
            R const length = (phi_low - phi_high) / pi();
            for (std::size_t i = 0; i < order; ++i)
            {
              R const phi = phi_high + (phi_low - phi_high) * (rule.x[i] + R{1}) / R{2};
              R const p = kinematics_.momentum(phi), other = pair - p;
              if (!(p > low && p < high && other > R{0} && other < pi()))
                throw std::underflow_error("resolved pair lost endpoint distance; increase precision");
              R const ja = kinematics_.jacobian(p), jb = kinematics_.jacobian(other);
              R const factor =
                  pi() / R(2 * factorial(n)) * rule.w[i] / (R{2} * R(rows) * R(rows)) * length * ja * jb / (ja + jb);
              R const energy = efree + kinematics_.energy(p) + kinematics_.energy(other);
              R const phi_other = kinematics_.angle(other);
              for (unsigned swap = 0; swap < 2; ++swap)
              {
                std::array<R, 4> angles{};
                if (n == 4)
                {
                  angles[0] = free_phi[0][row];
                  angles[1] = free_phi[1][col];
                }
                angles[n - 2] = swap ? phi_other : phi;
                angles[n - 1] = swap ? phi : phi_other;
                auto const amplitudes = vertex_(std::span<R const>{angles.data(), n});
                R weight = kinematics_.zero(), screen = kinematics_.zero();
                for (unsigned l = 0; l < amplitudes.polarizations; ++l)
                {
                  R const a = uni20::math::abs(amplitudes.value[0][l] - R(sign) * amplitudes.value[1][l]);
                  R const e = amplitudes.arithmetic_indicator[0][l] + amplitudes.arithmetic_indicator[1][l];
                  weight += a * a;
                  screen += R{2} * a * e + e * e;
                }
                weight *= factor;
                screen *= factor;
                if (!uni20::isfinite(weight) || !uni20::isfinite(screen))
                  throw std::overflow_error("invalid spectral contribution; no zero substitution");
                moments[0].add(weight);
                moments[1].add(energy * weight);
                moments[2].add(screen);
                moments[3].add(energy * screen);
                moments[4 + sector].add(weight);
                ++out.moments.nodes;
                for (std::size_t j = 0; j < omega.size(); ++j)
                {
                  // Energy and width are positive. Scale opposite signs
                  // separately so a finite Gaussian argument is not lost to
                  // overflow in omega-energy; retain the difference otherwise.
                  R const x = omega[j] < R{0} ? omega[j] / width - energy / width : (omega[j] - energy) / width;
                  R const gaussian = gaussian_normalization * uni20::math::exp(-x * x / R{2}); // includes 2pi
                  density[j].add(gaussian * weight);
                  error[j].add(gaussian * screen);
                }
              }
            }
          }
        }
        out.moments.weight = moments[0].value();
        out.moments.first_moment = moments[1].value();
        out.moments.arithmetic_indicator = moments[2].value();
        out.moments.first_moment_arithmetic_indicator = moments[3].value();
        out.sector_weights = {moments[4].value(), moments[5].value()};
        for (std::size_t j = 0; j < omega.size(); ++j)
        {
          out.density.push_back(density[j].value());
          out.arithmetic_indicator.push_back(error[j].value());
        }
      });
      VertexCut<R> result;
      std::vector<bethe::detail::CompensatedSum<R>> density(omega.size(),
                                                            bethe::detail::CompensatedSum<R>(kinematics_.zero())),
          error(omega.size(), bethe::detail::CompensatedSum<R>(kinematics_.zero()));
      std::array<bethe::detail::CompensatedSum<R>, 6> sums;
      sums.fill(bethe::detail::CompensatedSum<R>(kinematics_.zero()));
      for (auto const& out : partial)
      {
        for (std::size_t j = 0; j < omega.size(); ++j)
        {
          density[j].add(out.density.at(j));
          error[j].add(out.arithmetic_indicator.at(j));
        }
        sums[0].add(out.moments.weight);
        sums[1].add(out.moments.first_moment);
        sums[2].add(out.moments.arithmetic_indicator);
        sums[3].add(out.moments.first_moment_arithmetic_indicator);
        sums[4].add(out.sector_weights[0]);
        sums[5].add(out.sector_weights[1]);
        result.moments.nodes += out.moments.nodes;
      }
      result.moments.weight = sums[0].value();
      result.moments.first_moment = sums[1].value();
      result.moments.arithmetic_indicator = sums[2].value();
      result.moments.first_moment_arithmetic_indicator = sums[3].value();
      result.sector_weights = {sums[4].value(), sums[5].value()};
      for (auto const& sum : sums)
        require_finite(sum.value());
      for (std::size_t j = 0; j < omega.size(); ++j)
      {
        require_finite(density[j].value());
        require_finite(error[j].value());
        result.density.push_back(density[j].value());
        result.arithmetic_indicator.push_back(error[j].value());
      }
      return result;
    }
    VertexCut<R> broadened_cut(unsigned n, R momentum, std::span<R const> omega, R width, std::size_t order,
                               R offset = R{0}) const
    {
      return broadened_cut(
          n, momentum, omega, width, order,
          [](std::size_t rows, auto body) {
            for (std::size_t j = 0; j < rows; ++j)
              body(j);
          },
          offset);
    }

  private:
    R checked_offset(R value) const
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (value.initialized() && value.is_exact() && value == R{0}) value = kinematics_.zero();
      bethe::detail::require_same_precision(value, kinematics_.zero());
      return value;
    }
    static void require_finite(R value)
    {
      if (!uni20::isfinite(value)) throw std::overflow_error("nonfinite spectral sum; no zero substitution");
    }
    void validate_energy_range(unsigned n) const
    {
      if (!uni20::isfinite(R(n) * kinematics_.maximum_energy()))
        throw std::overflow_error("spinon energy sum outside scalar range");
    }
    R pi() const { return bethe::detail::spinon_pi(kinematics_.zero()); }
    static unsigned factorial(unsigned n) { return n == 2 ? 2 : 24; }
    static void validate(unsigned n, std::size_t order, R offset, bool paired)
    {
      if ((n != 2 && n != 4) || order < 2 || order > 512 || !uni20::isfinite(offset) || (paired && order % 2))
        throw std::invalid_argument("require spinons 2/4, order 2..512 (even for pair cuts), finite offset");
    }
    R node_angle(std::size_t i, unsigned dimension, std::size_t order, R offset) const
    {
      static constexpr unsigned shift[4] = {173, 389, 617, 823};
      R const raw = R(shift[dimension]) / R{1000} + offset, frac = raw - uni20::math::floor(raw);
      return pi() * ((R(i) + frac) / R(order) - (R{1} / R{2}));
    }
    SpectralKinematics<R> kinematics_;
    VertexEvaluator<R> vertex_;
};
} // namespace bethe::experimental::xxz
