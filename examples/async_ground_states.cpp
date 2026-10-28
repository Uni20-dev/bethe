// SPDX-License-Identifier: GPL-3.0-or-later
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/observable_grid.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>
#include <uni20/core/scalar_io.hpp>

namespace model = bethe::lieb_liniger::thermo;

template <uni20::Real Real> int sweep(model::Options<Real> options)
{
  uni20::async::TbbScheduler scheduler(4);
  std::vector<model::Background<Real>> results;
  auto complete = bethe::evaluate_observable_grid_async(
      4,
      [options](std::size_t i) {
        auto const c = uni20::scalar_at_precision<Real>(4 * (i + 1), options.precision());
        auto const density = uni20::scalar_at_precision<Real>(1, options.precision());
        return model::ground_state(c, density, options);
      },
      [&results](std::size_t, auto result) { results.push_back(std::move(result)); },
      {.scheduler = &scheduler, .max_in_flight = 2});

  // One final result covers every admitted solve and the ordered collection,
  // including on failure. Numerical nonconvergence remains in each result.
  (void)complete.get_wait(scheduler);
  std::cout << "# Lieb-Liniger ground states; density=1; working bits=" << options.precision().bit_count()
            << "\n# c energy_per_length chemical_potential\n";
  for (auto const& result : results)
  {
    if (!result.converged) throw std::runtime_error("ground state did not converge within its numerical budget");
    std::cout << uni20::format_real(result.interaction) << ' ' << uni20::format_real(*result.energy_per_length) << ' '
              << uni20::format_real(*result.chemical_potential) << '\n';
  }
  return 0;
}

int main()
{
  try
  {
#if UNI20_HAS_MPLAPACK_MPREAL
    return sweep(model::Options<uni20::mpreal>(uni20::Precision::bits(160)));
#else
    return sweep(model::Options<double>{});
#endif
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
