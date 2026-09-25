// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <bit>
#include <complex>
#include <stdexcept>
#include <uni20/linalg/ops/nonsymmetric_eigen.hpp>
#include <uni20/tensor/tensor.hpp>
#include <vector>

namespace test_support
{
// Independent physical configuration-space generator; deliberately test-only.
inline std::vector<std::complex<double>> exclusion_spectrum(unsigned l, unsigned n, double right = 1, double left = 0)
{
  if (l < 2 || l > 12 || n > l) throw std::invalid_argument("small exclusion reference requires 2<=L<=12, N<=L");
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1u << l);
  for (unsigned a = 0; a < (1u << l); ++a)
    if (std::popcount(a) == int(n))
    {
      index[a] = basis.size();
      basis.push_back(a);
    }
  uni20::DenseMatrix<double> m(basis.size(), basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      m[i, j] = 0;
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (unsigned j = 0; j < l; ++j)
      for (bool forward : {true, false})
      {
        auto const k = forward ? (j + 1) % l : (j + l - 1) % l, a = basis[i];
        double const rate = forward ? right : left;
        if ((a & (1u << j)) && !(a & (1u << k)))
        {
          m[index[a ^ (1u << j) ^ (1u << k)], i] += rate;
          m[i, i] -= rate;
        }
      }
  std::vector<std::complex<double>> values(basis.size());
  uni20::DenseMatrix<std::complex<double>> vectors(basis.size(), basis.size());
  uni20::linalg::nonsymmetric_eigen(m, std::span<std::complex<double>>(values), vectors, false);
  std::sort(values.begin(), values.end(), [](auto a, auto b) { return a.real() > b.real(); });
  return values;
}
} // namespace test_support
