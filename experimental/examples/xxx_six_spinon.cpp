// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/xxx_six_spinon.hpp>
#include <iostream>
#include <uni20/core/scalar_io.hpp>

int main()
{
  namespace xxx = bethe::heisenberg;
  namespace experimental = bethe::experimental::xxx;
  namespace six = experimental::six_spinon;
  auto const singlet = xxx::sector_ground_state<double>(12, uni20::half_int{0});
  experimental::PreparedSinglet<double> const ground(12, singlet);
  if (!ground.ready()) return 1;
  six::Configuration const request{
      12, six::Family::two_pairs, {-3, -2, -1, 0, 1, 3}, {uni20::from_twice(-1), uni20::from_twice(3)}};
  auto const state = six::solve<double>(request);
  if (!state.converged())
  {
    std::cerr << state.diagnostic << '\n';
    return 1;
  }
  auto const factor = six::raising_weight(ground, state);
  if (!factor.converged()) return 1;
  std::cout << "momentum_index = " << state.state->momentum_index << " of 12\n"
            << "energy = " << uni20::format_scalar(state.state->energy) << '\n'
            << "gap = " << uni20::format_scalar(state.state->energy - *ground.energy()) << '\n'
            << "Fourier raising weight = " << uni20::format_scalar(*factor.weight) << '\n';
}
