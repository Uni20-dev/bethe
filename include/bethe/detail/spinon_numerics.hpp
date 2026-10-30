// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <cstddef>
#include <stdexcept>
#include <string>
#include <uni20/core/math.hpp>
#include <utility>

namespace bethe::detail
{
enum class SpinonFunctionFailure
{
  arithmetic,
  budget
};
class SpinonFunctionError : public std::runtime_error {
  public:
    SpinonFunctionError(SpinonFunctionFailure failure, std::string const& message)
        : std::runtime_error(message), failure_(failure)
    {}
    SpinonFunctionFailure failure() const noexcept { return failure_; }

  private:
    SpinonFunctionFailure failure_;
};

template <uni20::RealOrComplex V> struct SpinonFunctionEstimate
{
    V value;
    uni20::make_real_t<V> truncation_bound, roundoff_indicator;
    std::size_t terms, shifts;
};

// Neumaier's correction, including when a new contour term exceeds the
// accumulated magnitude. The Kahan accumulator used by other Bethe solvers
// has different cancellation behaviour; retain this algorithm for residue sums.
template <class T> class SpinonSum {
  public:
    explicit SpinonSum(T const& initial) : value_(initial), correction_(initial - initial) {}
    void add(T const& x)
    {
      T next = value_ + x;
      if (uni20::math::abs(value_) >= uni20::math::abs(x))
        correction_ += (value_ - next) + x;
      else
        correction_ += (x - next) + value_;
      value_ = std::move(next);
    }
    T value() const { return value_ + correction_; }

  private:
    T value_, correction_;
};
} // namespace bethe::detail
