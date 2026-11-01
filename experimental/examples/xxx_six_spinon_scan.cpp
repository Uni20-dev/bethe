// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/xxx_six_spinon_scan.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>
#include <uni20/core/scalar_io.hpp>

int main()
{
  namespace six = bethe::experimental::xxx::six_spinon;
  namespace xxx = bethe::experimental::xxx;
  // Preflight labels and storage before preparing the ground or solving states.
  six::MomentumCatalogue const catalogue(12, 3); // q=pi/2, absolute Pn=3*pi/2.
  auto requests = catalogue.requests(100);
  auto const state = bethe::heisenberg::sector_ground_state<double>(12, uni20::half_int{0});
  xxx::PreparedSinglet<double> ground(12, state);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = six::scan_async(std::move(ground), std::move(requests), {&scheduler, 1, 2});
  auto const& result = pending.get_wait(scheduler);
  std::cout << "N=12, q=pi/2; J=1; raising Fourier weights, no rescaling\n"
            << "nominal label candidates: " << result.records.size() << "\naccepted lines: " << result.lines.size()
            << "\nThis selected family is not the complete six-spinon spectrum.\n";
  for (auto const& row : result.records)
  {
    std::cout << row.candidate_id;
    if (row.line)
      std::cout << " gap=" << uni20::format_scalar(row.line->gap)
                << " raising_weight=" << uni20::format_scalar(row.line->weight);
    else
      std::cout << " unavailable: " << row.diagnostic;
    std::cout << '\n';
  }
  std::cout << "accepted m0=" << uni20::format_scalar(result.moments[3].weight)
            << " m1=" << uni20::format_scalar(result.moments[3].first_moment) << '\n';
  return result.converged() ? 0 : 2;
}
