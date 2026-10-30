// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/detail/spinon_numerics.hpp>
#include <bethe/solver.hpp>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace bethe
{
/// Exact midpoint uniforms on a declared binary grid. Retain this identity
/// when retrying geometry at a higher precision: a retry must not redraw.
template <std::size_t Dimension> class UniformSample {
  public:
    using Words = std::array<std::vector<std::uint64_t>, Dimension>;

    UniformSample(unsigned bits, Words words) : bits_(bits), words_(std::move(words))
    {
      if (!bits_ || bits_ >= unsigned(std::numeric_limits<int>::max()))
        throw std::invalid_argument("invalid uniform bit count");
      std::size_t const count = (std::size_t(bits_) + 63) / 64;
      unsigned const first = (bits_ - 1) % 64 + 1;
      for (auto const& coordinate : words_)
        if (coordinate.size() != count || (first < 64 && coordinate.front() >> first))
          throw std::invalid_argument("uniform words disagree with bit count");
    }

    unsigned bits() const { return bits_; }
    Words const& words() const { return words_; }

    template <uni20::Real R> std::array<R, Dimension> at_precision(R const& exemplar) const
    {
      R const one = uni20::scalar_like(exemplar, 1), zero = uni20::scalar_like(exemplar, 0);
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (one.precision().is_exact()) throw std::invalid_argument("uniform conversion requires finite precision");
      if (uni20::numeric_limits<R>::digits(one) < bits_ + 1)
        throw std::invalid_argument("precision cannot represent the declared midpoint grid");
      std::array<R, Dimension> result;
      for (std::size_t d = 0; d < Dimension; ++d)
      {
        R value = zero;
        for (auto const word : words_[d])
          value = uni20::math::ldexp(value, 64) + uni20::scalar_like(one, word);
        result[d] = uni20::math::ldexp(value + one / R{2}, -int(bits_));
        if (!(zero < result[d] && result[d] < one)) throw std::runtime_error("uniform midpoint lost open support");
      }
      return result;
    }

  private:
    unsigned bits_;
    Words words_;
};

/// One task-local stream per replicate. Its identity is (master seed,
/// replicate id, uniform bits); completion order never changes the draws.
/// Uses the standardized mt19937_64/seed_seq engines, not implementation-defined
/// uniform_real_distribution. This is ordinary Monte Carlo, not scrambled QMC.
template <std::size_t Dimension> class MonteCarloStream {
  public:
    MonteCarloStream(std::uint64_t seed, std::uint64_t replicate, unsigned bits) : bits_(bits)
    {
      if (!bits_ || bits_ >= unsigned(std::numeric_limits<int>::max()))
        throw std::invalid_argument("invalid uniform bit count");
      std::seed_seq sequence{std::uint32_t(0x5350494e), std::uint32_t(1),         std::uint32_t(seed),
                             std::uint32_t(seed >> 32), std::uint32_t(replicate), std::uint32_t(replicate >> 32)};
      engine_.seed(sequence);
    }

    UniformSample<Dimension> next()
    {
      typename UniformSample<Dimension>::Words words;
      std::size_t const count = (std::size_t(bits_) + 63) / 64;
      unsigned const first = (bits_ - 1) % 64 + 1;
      for (auto& coordinate : words)
      {
        coordinate.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
          coordinate.push_back(engine_());
        if (first < 64) coordinate.front() &= (std::uint64_t{1} << first) - 1;
      }
      return {bits_, std::move(words)};
    }

  private:
    unsigned bits_;
    std::mt19937_64 engine_;
};

enum class ReplicateStatus
{
  incomplete,
  failed,
  complete
};

template <uni20::Real R> struct ReplicateSummary
{
    ReplicateStatus status;
    std::size_t received, failed;
    std::optional<std::vector<R>> mean;
    // Row-major covariance OF THE MEAN across complete, equally weighted
    // independent replicates. Absent for fewer than two replicates.
    std::optional<std::vector<R>> covariance;
    std::vector<std::string> failures;
};

/// A collector for complete equal-size replicate means. Missing/failed
/// replicates never become zeros or a successful mean over the survivors.
/// Writes belong to the collector thread; workers return their own results.
template <uni20::Real R> class ReplicateStatistics {
  public:
    ReplicateStatistics(std::size_t replicates, std::size_t samples, std::size_t dimension, R const& exemplar)
        : samples_(samples), dimension_(dimension), zero_(uni20::scalar_like(exemplar, 0)), rows_(replicates)
    {
      if (!replicates || !samples || !dimension || dimension > std::numeric_limits<std::size_t>::max() / dimension)
        throw std::invalid_argument("invalid replicate dimensions or work count");
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (zero_.precision().is_exact()) throw std::invalid_argument("statistics require explicit finite precision");
    }

    void accept(std::size_t replicate, std::size_t samples, std::vector<R> values)
    {
      auto& row = vacant(replicate);
      if (samples != samples_ || values.size() != dimension_)
        throw std::invalid_argument("replicate is incomplete or has inconsistent dimensions");
      for (auto const& value : values)
      {
        if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite replicate mean");
        bethe::detail::require_same_precision(value, zero_);
      }
      row.values = std::move(values);
    }

    void reject(std::size_t replicate, std::string reason)
    {
      auto& row = vacant(replicate);
      if (reason.empty()) throw std::invalid_argument("failed replicate requires a reason");
      row.failure = std::move(reason);
    }

    ReplicateSummary<R> summary() const
    {
      ReplicateSummary<R> result{ReplicateStatus::incomplete, 0, 0, {}, {}, {}};
      for (std::size_t i = 0; i < rows_.size(); ++i)
      {
        auto const& row = rows_[i];
        if (row.values || !row.failure.empty()) ++result.received;
        if (!row.failure.empty())
        {
          ++result.failed;
          result.failures.push_back("replicate " + std::to_string(i) + ": " + row.failure);
        }
      }
      if (result.failed) result.status = ReplicateStatus::failed;
      if (result.received != rows_.size() || result.failed) return result;

      // Reduce in logical id order even when work finished out of order.
      auto const& anchor = *rows_.front().values;
      R const count = uni20::scalar_like(zero_, rows_.size());
      std::vector<R> shifts, means;
      for (std::size_t j = 0; j < dimension_; ++j)
      {
        bethe::detail::SpinonSum<R> sum(zero_);
        for (auto const& row : rows_)
          sum.add((*row.values)[j] - anchor[j]);
        shifts.push_back(sum.value() / count);
        if (sum.value() != zero_ && shifts.back() == zero_)
          throw std::runtime_error("replicate mean shift underflowed");
        means.push_back(anchor[j] + shifts.back());
        if (!uni20::isfinite(means.back())) throw std::runtime_error("replicate mean exceeded exponent range");
      }
      if (rows_.size() > 1)
      {
        R const denominator = count * (count - R{1});
        std::vector<R> covariance(dimension_ * dimension_, zero_);
        for (std::size_t j = 0; j < dimension_; ++j)
          for (std::size_t k = 0; k <= j; ++k)
          {
            bethe::detail::SpinonSum<R> sum(zero_);
            for (auto const& row : rows_)
            {
              R const a = ((*row.values)[j] - anchor[j]) - shifts[j];
              R const b = ((*row.values)[k] - anchor[k]) - shifts[k];
              R const product = a * b;
              if (a != zero_ && b != zero_ && product == zero_)
                throw std::runtime_error("replicate covariance product underflowed");
              sum.add(product);
            }
            R const value = sum.value() / denominator;
            if (!uni20::isfinite(value) || (j == k && value < zero_) || (sum.value() != zero_ && value == zero_))
              throw std::runtime_error("replicate covariance unresolved");
            covariance[j * dimension_ + k] = covariance[k * dimension_ + j] = value;
          }
        result.covariance = std::move(covariance);
      }
      result.mean = std::move(means);
      result.status = ReplicateStatus::complete;
      return result;
    }

  private:
    struct Row
    {
        std::optional<std::vector<R>> values;
        std::string failure;
    };
    Row& vacant(std::size_t id)
    {
      if (id >= rows_.size()) throw std::out_of_range("replicate id outside requested work");
      if (rows_[id].values || !rows_[id].failure.empty()) throw std::invalid_argument("duplicate replicate id");
      return rows_[id];
    }
    std::size_t samples_, dimension_;
    R zero_;
    std::vector<Row> rows_;
};
} // namespace bethe
