// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/xxz_vertex_spectrum.hpp>
#include <bethe/observable_grid.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

int main()
{
  namespace model = bethe::experimental::xxz;
  using R = double;
  model::ThermodynamicVertexSpectrum<R> spectrum(1.5);
  std::array<R, 3> omega{1.8, 2.2, 2.6};
  uni20::async::TbbScheduler scheduler(2);
  bethe::ExecutionOptions execution{.scheduler = &scheduler, .max_in_flight = 2};
  auto execute = [&](std::size_t rows, auto const& body) {
    // The grid waits for its final result, including drainage on failure.
    // Independent rows are accumulated later in their fixed original order.
    bethe::evaluate_observable_grid(
        rows,
        [&](std::size_t row) {
          body(row);
          return row;
        },
        [](std::size_t, std::size_t) {}, execution);
  };
  std::cout << "Fixed-order finite-q demonstration: Delta=1.5, J=1, k=1.5, Gaussian width=.125, order=4\n"
               "Transverse raising response; product/residue/grid errors require independent refinement.\n"
               "spinons,omega,density,arithmetic_indicator\n";
  for (unsigned spinons : {2u, 4u})
  {
    auto const cut = spectrum.broadened_cut(spinons, 1.5, omega, .125, 4, execute);
    auto const serial = spectrum.broadened_cut(spinons, 1.5, omega, .125, 4);
    if (cut.density != serial.density || cut.moments.weight != serial.moments.weight)
      throw std::runtime_error("serial/Async reduction mismatch");
    for (std::size_t i = 0; i < omega.size(); ++i)
      std::cout << spinons << ',' << omega[i] << ',' << uni20::format_scalar(cut.density[i]) << ','
                << uni20::format_scalar(cut.arithmetic_indicator[i]) << '\n';
    auto const full = spectrum.integrated(spinons, 4);
    std::cout << "# spinons=" << spinons << ", fixed-k weight=" << uni20::format_scalar(cut.moments.weight)
              << ", full-zone weight=" << uni20::format_scalar(full.weight)
              << ", full-zone first moment=" << uni20::format_scalar(full.first_moment) << '\n';
  }
}
