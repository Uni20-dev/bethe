// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bit>
#include <cmath>
#include <stdexcept>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <vector>

namespace bethe::test
{
struct TJExactGround
{
    double energy;
    std::vector<double> momentum_weights;
};

// Independent projected Fock-space oracle. Site-major orbital order:
// (0 up,0 down,1 up,1 down,...). Boundary hops retain fermionic signs.
inline TJExactGround tj_exact_ground(unsigned n, unsigned up, unsigned down)
{
  if (n < 3 || n > 9 || up > n || down > n - up) throw std::invalid_argument("small t-J ED oracle domain");
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1U << (2 * n));
  for (unsigned a = 0; a < (1U << n); ++a)
    if (std::popcount(a) == int(up))
      for (unsigned b = 0; b < (1U << n); ++b)
        if (!(a & b) && std::popcount(b) == int(down))
        {
          unsigned bits = 0;
          for (unsigned j = 0; j < n; ++j)
            bits |= (((a >> j) & 1) << (2 * j)) | (((b >> j) & 1) << (2 * j + 1));
          index[bits] = basis.size();
          basis.push_back(bits);
        }
  auto const dim = basis.size();
  uni20::DenseMatrix<double> h(dim, dim);
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      h[i, j] = 0;
  auto bilinear = [](unsigned bits, unsigned from, unsigned to) {
    int parity = std::popcount(bits & ((1U << from) - 1));
    bits ^= 1U << from;
    parity += std::popcount(bits & ((1U << to) - 1));
    return std::pair{bits | (1U << to), parity % 2 ? -1.0 : 1.0};
  };
  for (std::size_t c = 0; c < dim; ++c)
    for (unsigned j = 0; j < n; ++j)
    {
      unsigned const k = (j + 1) % n, bits = basis[c], nj = (bits >> (2 * j)) & 3, nk = (bits >> (2 * k)) & 3;
      if (nj && nk && nj != nk)
      {
        h[c, c] -= 1; // J*(SzSz-nn/4), J=2.
        auto const a = bilinear(bits, 2 * k + (nk == 2), 2 * k + (nk == 1));
        auto const b = bilinear(a.first, 2 * j + (nj == 2), 2 * j + (nj == 1));
        h[index[b.first], c] += a.second * b.second; // J/2*(S+S-+S-S+).
      }
      for (unsigned spin = 0; spin < 2; ++spin)
        for (bool reverse : {false, true})
        {
          unsigned const from = reverse ? k : j, to = reverse ? j : k;
          if (((bits >> (2 * to)) & 3) || !(bits & (1U << (2 * from + spin)))) continue;
          auto const hop = bilinear(bits, 2 * from + spin, 2 * to + spin);
          h[index[hop.first], c] -= hop.second;
        }
    }
  auto const eig = uni20::linalg::eigh(std::move(h));
  std::vector<std::size_t> translation(dim);
  std::vector<double> sign(dim), weights(n, 0);
  for (std::size_t j = 0; j < dim; ++j)
  {
    unsigned const bits = basis[j], wrapped = bits >> (2 * (n - 1));
    translation[j] = index[((bits << 2) & ((1U << (2 * n)) - 1)) | wrapped];
    sign[j] = wrapped && (up + down - 1) % 2 ? -1 : 1;
  }
  double const pi = 4 * std::atan(1.0);
  for (std::size_t v = 0; v < dim && eig.eigenvalues[v] - eig.eigenvalues[0] < 1e-9; ++v)
    for (std::size_t col = 0; col < dim; ++col)
    {
      auto row = col;
      double parity = 1;
      for (unsigned power = 0; power < n; ++power)
      {
        double const contribution = parity * eig.eigenvectors[row, v] * eig.eigenvectors[col, v] / n;
        for (unsigned q = 0; q < n; ++q)
          weights[q] += contribution * std::cos(2 * pi * q * power / n);
        parity *= sign[row];
        row = translation[row];
      }
    }
  return {eig.eigenvalues[0], std::move(weights)};
}
} // namespace bethe::test
