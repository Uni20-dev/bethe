// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace test_support
{
// Independent small-chain oracle: construct H directly in the Sz bit basis.
// Optionally add a*(T+T^-1)/2, which shifts an eigenvalue by a*cos(P), to
// check momentum as well as energy without a momentum-space Bethe formula.
// Double is intentional only in this ED oracle; separate analytic regressions
// in the test programs check the solver at the selected precision.
inline std::vector<double> exact_spectrum(unsigned n, unsigned down, double translation_weight = 0,
                                          bool periodic = true)
{
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1U << n);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == static_cast<int>(down))
    {
      index[bits] = basis.size();
      basis.push_back(bits);
    }
  auto const dim = basis.size();
  std::vector<double> matrix(dim * dim, 0);
  auto at = [&](std::size_t i, std::size_t j) -> double& { return matrix[i * dim + j]; };
  for (std::size_t i = 0; i < dim; ++i)
  {
    auto const bits = basis[i];
    for (unsigned site = 0; site < (periodic ? n : n - 1); ++site)
    {
      unsigned const next = (site + 1) % n;
      bool const opposite = ((bits >> site) & 1U) != ((bits >> next) & 1U);
      at(i, i) += opposite ? -0.25 : 0.25;
      if (opposite) at(index[bits ^ (1U << site) ^ (1U << next)], i) += 0.5;
    }
    unsigned const translated = ((bits << 1) & ((1U << n) - 1)) | (bits >> (n - 1));
    at(index[translated], i) += translation_weight / 2;
    at(i, index[translated]) += translation_weight / 2;
  }
  // Cyclic Jacobi diagonalization of a real symmetric matrix.
  bool diagonal = false;
  for (int sweep = 0; sweep < 80; ++sweep)
  {
    double largest = 0;
    for (std::size_t p = 0; p < dim; ++p)
      for (std::size_t q = p + 1; q < dim; ++q)
      {
        double const off = at(p, q);
        largest = std::max(largest, std::abs(off));
        if (std::abs(off) < 1e-14) continue;
        double const tau = (at(q, q) - at(p, p)) / (2 * off);
        double const t = std::copysign(1.0, tau) / (std::abs(tau) + std::hypot(1.0, tau));
        double const c = 1 / std::sqrt(1 + t * t);
        double const s = t * c;
        at(p, p) -= t * off;
        at(q, q) += t * off;
        at(p, q) = at(q, p) = 0;
        for (std::size_t r = 0; r < dim; ++r)
          if (r != p && r != q)
          {
            double const rp = at(r, p), rq = at(r, q);
            at(r, p) = at(p, r) = c * rp - s * rq;
            at(r, q) = at(q, r) = s * rp + c * rq;
          }
      }
    if (largest < 1e-13)
    {
      diagonal = true;
      break;
    }
  }
  if (!diagonal) throw std::runtime_error("ED oracle diagonalization failed");
  std::vector<double> energies(dim);
  for (std::size_t i = 0; i < dim; ++i)
    energies[i] = at(i, i);
  std::sort(energies.begin(), energies.end());
  return energies;
}
} // namespace test_support
