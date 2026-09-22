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
inline double richardson_exact_ground(std::span<double const> levels, unsigned pairs, double g,
                                      std::span<std::size_t const> blocked = {})
{
  if (levels.size() > 10) throw std::invalid_argument("small Richardson ED oracle domain");
  std::vector<double> active;
  double offset = 0;
  for (std::size_t j = 0; j < levels.size(); ++j)
    if (std::find(blocked.begin(), blocked.end(), j) != blocked.end())
      offset += levels[j];
    else
      active.push_back(levels[j]);
  if (pairs > active.size()) throw std::invalid_argument("invalid ED pair count");
  std::vector<unsigned> basis, index(1U << active.size());
  for (unsigned s = 0; s < index.size(); ++s)
    if (std::popcount(s) == int(pairs))
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
    h[col, col] = offset - g * pairs; // Includes diagonal i=j pair scattering.
    for (unsigned i = 0; i < active.size(); ++i)
      if (s & (1U << i))
      {
        h[col, col] += 2 * active[i];
        for (unsigned j = 0; j < active.size(); ++j)
          if (!(s & (1U << j))) h[index[s ^ (1U << i) ^ (1U << j)], col] -= g;
      }
  }
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}
} // namespace bethe::test
