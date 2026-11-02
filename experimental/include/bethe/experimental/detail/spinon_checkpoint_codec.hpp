// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxx_spinon_density_recovery.hpp>
#include <bit>
#include <limits>

namespace bethe::experimental::xxx
{
struct RecoveryPoint
{
    std::string id;
    std::vector<std::string> beta;
};

struct RecoveryChunkBinding
{
    /// Caller-supplied executable/kernel revision and complete request provenance.
    std::string producer, provenance;
    heisenberg::SpinonRecoveryPolicy policy;
    std::vector<RecoveryPoint> points;
};

namespace checkpoint_detail
{
inline constexpr std::size_t max_points = 64, max_stages = 16, max_bytes = 16 * 1024 * 1024;
inline constexpr std::size_t max_string = 65536;
inline constexpr int max_bits = 65536;

inline void require(bool condition, char const* message)
{
  if (!condition) throw std::invalid_argument(message);
}

// Detect accidental truncation/corruption, not hostile edits or rollback.
inline std::uint64_t checksum(std::string_view text)
{
  std::uint64_t hash = 14695981039346656037ULL;
  for (unsigned char c : text)
    hash = (hash ^ c) * 1099511628211ULL;
  return hash;
}

struct Encoder
{
    std::string data;
    void bytes(std::string_view s)
    {
      if (s.size() > max_bytes - data.size()) throw std::length_error("recovery chunk exceeds 16 MiB");
      data.append(s);
    }
    void integer(std::uint64_t x)
    {
      char buffer[8];
      for (int i = 0; i < 8; ++i)
        buffer[i] = static_cast<char>((x >> (8 * i)) & 255);
      bytes({buffer, 8});
    }
    void text(std::string_view s)
    {
      if (s.size() > max_string) throw std::length_error("recovery string exceeds 64 KiB");
      integer(s.size());
      bytes(s);
    }
    void clock(double x)
    {
      static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
      require(std::isfinite(x) && x >= 0, "invalid checkpoint time");
      integer(std::bit_cast<std::uint64_t>(x));
    }
    template <uni20::Real R> void real(R const& x)
    {
      auto const p = uni20::precision_of(x);
      require(!p.is_exact() && uni20::isfinite(x), "checkpoint requires finite approximate scalars");
      integer(p.bit_count());
      auto const parts = uni20::math::frexp(x);
      integer(std::bit_cast<std::uint64_t>(parts.exponent));
      text(uni20::format_real(parts.fraction, {.normalize_negative_zero = false}));
    }
    std::string finish()
    {
      integer(checksum(data));
      return std::move(data);
    }
};

struct Decoder
{
    std::string_view data;
    std::size_t at = 0;
    std::string_view bytes(std::size_t n)
    {
      require(n <= data.size() - at, "truncated recovery chunk");
      auto result = data.substr(at, n);
      at += n;
      return result;
    }
    std::uint64_t integer()
    {
      auto s = bytes(8);
      std::uint64_t x = 0;
      for (int i = 0; i < 8; ++i)
        x |= std::uint64_t(static_cast<unsigned char>(s[i])) << (8 * i);
      return x;
    }
    std::size_t count(std::size_t limit)
    {
      auto n = integer();
      require(n <= limit, "recovery count exceeds bound");
      return static_cast<std::size_t>(n);
    }
    bool flag() { return count(1) != 0; }
    std::string text() { return std::string(bytes(count(max_string))); }
    double clock()
    {
      double x = std::bit_cast<double>(integer());
      require(std::isfinite(x) && x >= 0, "invalid checkpoint time");
      return x;
    }
    template <uni20::Real R> R real(int bits)
    {
      require(integer() == static_cast<std::uint64_t>(bits), "changed scalar precision");
      auto exponent = std::bit_cast<std::int64_t>(integer());
      auto mantissa = text();
      auto fraction = [&] {
        if constexpr (uni20::has_runtime_precision_v<R>)
          return uni20::parse_real<R>(mantissa, uni20::Precision::bits(bits));
        else
        {
          require(uni20::numeric_limits<R>::digits() == bits, "native arithmetic changed");
          return uni20::parse_real<R>(mantissa);
        }
      }();
      require(uni20::isfinite(fraction), "nonfinite checkpoint mantissa");
      require(uni20::format_real(fraction, {.normalize_negative_zero = false}) == mantissa,
              "noncanonical checkpoint mantissa");
      auto value = uni20::math::ldexp(fraction, exponent);
      require(uni20::isfinite(value), "checkpoint exponent overflow");
      auto rebuilt = uni20::math::frexp(value);
      require(rebuilt.exponent == exponent &&
                  uni20::format_real(rebuilt.fraction, {.normalize_negative_zero = false}) == mantissa,
              "checkpoint scalar did not round trip");
      return value;
    }
};

template <uni20::Real R> void arithmetic_tag(Encoder& e)
{
  e.integer(uni20::has_runtime_precision_v<R>);
  if constexpr (!uni20::has_runtime_precision_v<R>)
  {
    e.integer(uni20::numeric_limits<R>::digits());
    e.integer(std::bit_cast<std::uint64_t>(std::int64_t(uni20::numeric_limits<R>::min_exponent)));
    e.integer(uni20::numeric_limits<R>::max_exponent);
  }
}

template <uni20::Real Low, uni20::Real High> std::string encode_binding(RecoveryChunkBinding const& b)
{
  require(!b.producer.empty() && !b.provenance.empty(), "producer and request provenance are required");
  require(!b.points.empty() && b.points.size() <= max_points, "recovery chunk requires 1 to 64 points");
  require(!b.policy.stages.empty() && b.policy.stages.size() <= max_stages, "recovery chunk requires 1 to 16 stages");
  Encoder e;
  arithmetic_tag<Low>(e);
  arithmetic_tag<High>(e);
  e.text(b.producer);
  e.text(b.provenance);
  e.text(b.policy.relative_tolerance);
  e.text(b.policy.absolute_tolerance);
  e.clock(b.policy.max_wall_seconds);
  e.integer(b.policy.stages.size());
  for (auto const& stage : b.policy.stages)
  {
    require(stage.bits > 0 && stage.bits <= max_bits && stage.guard_bits > 0 && stage.guard_bits <= max_bits,
            "checkpoint precision must be in [1,65536] bits");
    e.integer(stage.bits);
    e.integer(stage.guard_bits);
    for (auto const& control : stage.contours)
    {
      e.integer(control.residue_terms);
      e.integer(control.tail_order);
    }
  }
  e.integer(b.points.size());
  std::set<std::string> identities;
  for (auto const& point : b.points)
  {
    require(!point.id.empty() && identities.insert(point.id).second, "empty or duplicate recovery point identity");
    require(point.beta.size() == 2 || point.beta.size() == 4 || point.beta.size() == 6,
            "recovery points require two, four or six rapidities");
    e.text(point.id);
    e.integer(point.beta.size());
    std::set<std::string> rapidities;
    for (auto const& value : point.beta)
    {
      require(rapidities.insert(heisenberg::detail::spinon_decimal_key(value)).second,
              "coincident rapidities need an unsupported limit");
      e.text(value);
    }
  }
  return std::move(e.data);
}

template <uni20::Real R> void encode_attempt(Encoder& e, heisenberg::SpinonRecoveryAttempt<R> const& a)
{
  e.integer(a.stage);
  e.integer(a.control);
  e.integer(static_cast<int>(a.failure));
  e.clock(a.completed_wall_seconds);
  e.text(a.reason);
  e.integer(a.checked.has_value());
  if (!a.checked) return;
  auto const& c = *a.checked;
  auto const& x = c.estimate;
  for (auto const* v :
       {&x.value, &x.contour_tail_error, &x.arithmetic_indicator, &x.log_scale, &x.contour_norm, &x.contour_norm_error,
        &x.contour_arithmetic_indicator, &c.precision_difference, &c.scalar_precision_difference, &c.arithmetic_screen,
        &c.requested_error, &c.low_value, &c.low_log_scale, &c.low_arithmetic_indicator})
    e.real(*v);
  e.integer(x.terms);
  e.integer(x.order);
  e.integer(x.amplitudes.size());
  for (auto const& amplitude : x.amplitudes)
  {
    e.real(amplitude.real());
    e.real(amplitude.imag());
  }
  e.integer(static_cast<int>(c.failure));
}

template <uni20::Real R>
heisenberg::SpinonRecoveryAttempt<R> decode_attempt(Decoder& d, heisenberg::SpinonRecoveryPolicy const& policy)
{
  using Failure = heisenberg::SpinonDensityFailure;
  heisenberg::SpinonRecoveryAttempt<R> a;
  a.stage = d.count(policy.stages.size() - 1);
  a.control = d.count(1);
  a.failure = static_cast<Failure>(d.count(static_cast<int>(Failure::rounded_collision)));
  a.completed_wall_seconds = d.clock();
  a.reason = d.text();
  if (!d.flag()) return a;
  heisenberg::CheckedSpinonDensity<R> c;
  auto& x = c.estimate;
  for (auto* v :
       {&x.value, &x.contour_tail_error, &x.arithmetic_indicator, &x.log_scale, &x.contour_norm, &x.contour_norm_error,
        &x.contour_arithmetic_indicator, &c.precision_difference, &c.scalar_precision_difference, &c.arithmetic_screen,
        &c.requested_error, &c.low_value, &c.low_log_scale, &c.low_arithmetic_indicator})
    *v = d.real<R>(policy.stages[a.stage].guard_bits);
  x.terms = static_cast<int>(d.count(std::numeric_limits<int>::max()));
  x.order = static_cast<int>(d.count(std::numeric_limits<int>::max()));
  auto n = d.count(15);
  x.amplitudes.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    auto re = d.real<R>(policy.stages[a.stage].guard_bits);
    auto im = d.real<R>(policy.stages[a.stage].guard_bits);
    x.amplitudes.emplace_back(std::move(re), std::move(im));
  }
  c.failure = static_cast<Failure>(d.count(static_cast<int>(Failure::rounded_collision)));
  a.checked = std::move(c);
  return a;
}
} // namespace checkpoint_detail
} // namespace bethe::experimental::xxx
