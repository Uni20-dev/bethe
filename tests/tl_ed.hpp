// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <stdexcept>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <utility>
#include <vector>

namespace bethe::test
{
inline std::vector<double> tl_ed_eigenvalues(uni20::DenseMatrix<double> h)
{
  auto const result = uni20::linalg::eigh(std::move(h));
  std::vector<double> values(result.eigenvalues.extent(0));
  for (std::size_t i = 0; i < values.size(); ++i)
    values[i] = result.eigenvalues[i];
  std::sort(values.begin(), values.end());
  return values;
}

// Independent full spin-1 basis: form -(S.S)^2 by multiplying the physical
// local spin matrix, NOT from the singlet projector or an XXZ energy shift.
inline std::vector<double> biquadratic_ed(unsigned n)
{
  if (n < 2 || n > 6) throw std::invalid_argument("small biquadratic ED oracle domain");
  std::vector<unsigned> powers(n + 1, 1);
  for (unsigned i = 1; i <= n; ++i)
    powers[i] = 3 * powers[i - 1];
  std::array<std::array<double, 9>, 9> dot{}, local{};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
    {
      dot[3 * a + b][3 * a + b] = (a - 1) * (b - 1);
      if (a < 2 && b > 0) dot[3 * (a + 1) + b - 1][3 * a + b] = 1;
      if (a > 0 && b < 2) dot[3 * (a - 1) + b + 1][3 * a + b] = 1;
    }
  for (unsigned a = 0; a < 9; ++a)
    for (unsigned b = 0; b < 9; ++b)
      for (unsigned k = 0; k < 9; ++k)
        local[a][b] -= dot[a][k] * dot[k][b];
  auto const dim = powers[n];
  uni20::DenseMatrix<double> h(dim, dim);
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      h[i, j] = 0;
  for (unsigned col = 0; col < dim; ++col)
    for (unsigned j = 0; j < n - 1; ++j)
    {
      auto const a = (col / powers[j]) % 3, b = (col / powers[j + 1]) % 3;
      for (unsigned aa = 0; aa < 3; ++aa)
        for (unsigned bb = 0; bb < 3; ++bb)
          h[col - a * powers[j] - b * powers[j + 1] + aa * powers[j] + bb * powers[j + 1], col] +=
              local[3 * aa + bb][3 * a + b];
    }
  return tl_ed_eigenvalues(std::move(h));
}

// Independent spin-half Sz basis, INCLUDING the opposite end fields.
inline std::vector<double> quantum_group_xxz_ed(unsigned n, unsigned down, double delta, bool fields = true)
{
  if (n < 2 || n > 10 || down > n) throw std::invalid_argument("small quantum-group XXZ ED oracle domain");
  std::vector<unsigned> basis, index(1U << n);
  for (unsigned s = 0; s < index.size(); ++s)
    if (std::popcount(s) == int(down))
    {
      index[s] = basis.size();
      basis.push_back(s);
    }
  auto const dim = basis.size();
  uni20::DenseMatrix<double> h(dim, dim);
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      h[i, j] = 0;
  for (std::size_t col = 0; col < dim; ++col)
  {
    auto const s = basis[col];
    // A set bit denotes a down spin.
    if (fields) h[col, col] = std::sqrt(delta * delta - 1) / 2 * (double((s >> (n - 1)) & 1U) - double(s & 1U));
    for (unsigned j = 0; j < n - 1; ++j)
    {
      bool const opposite = ((s >> j) & 1U) != ((s >> (j + 1)) & 1U);
      h[col, col] += delta * (opposite ? -0.25 : 0.25);
      if (opposite) h[index[s ^ (1U << j) ^ (1U << (j + 1))], col] += 0.5;
    }
  }
  return tl_ed_eigenvalues(std::move(h));
}
} // namespace bethe::test
