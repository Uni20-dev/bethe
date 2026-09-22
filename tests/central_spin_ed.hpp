// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bit>
#include <span>
#include <stdexcept>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <vector>

namespace bethe::test
{
inline double central_spin_exact_ground(std::span<double const> a, double b, unsigned m)
{
  if (a.size() > 9 || m > a.size() + 1) throw std::invalid_argument("small central-spin ED oracle domain");
  std::vector<unsigned> basis, index(1U << (a.size() + 1));
  for (unsigned s = 0; s < index.size(); ++s)
    if (std::popcount(s) == int(m))
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
    auto const s = basis[col], z = s & 1U;
    h[col, col] = b * (double(z) - 0.5);
    for (unsigned j = 0; j < a.size(); ++j)
    {
      auto const w = (s >> (j + 1)) & 1U;
      h[col, col] += a[j] * (double(z) - 0.5) * (double(w) - 0.5);
      if (z != w) h[index[s ^ 1U ^ (1U << (j + 1))], col] += a[j] / 2;
    }
  }
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}
} // namespace bethe::test
