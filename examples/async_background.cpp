// SPDX-License-Identifier: GPL-3.0-or-later
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

int main()
{
  namespace model = bethe::lieb_liniger::thermo;
  try
  {
    uni20::async::TbbScheduler scheduler(4);
    bethe::ExecutionOptions execution{.scheduler = &scheduler, .max_in_flight = 4};
    auto background = bethe::prepare_background_async([] { return model::Solver<double>(4.0, 1.0); }, execution);
    constexpr std::size_t points = 17;
    std::array<std::vector<model::Point<double>>, 2> curves;
    std::array<std::optional<uni20::async::Async<std::size_t>>, 2> completed;
    double const endpoint = 8.0 * std::atan(1.0); // 2*pi at density one
    std::exception_ptr error;
    try
    {
      for (std::size_t b = 0; b < curves.size(); ++b)
      {
        auto branch = b == 0 ? model::Branch::type_i : model::Branch::type_ii;
        completed[b].emplace(bethe::evaluate_observable_grid_async(
            background.read(), points,
            [branch, endpoint](model::Solver<double> const& solver, std::size_t i) {
              return solver.at_momentum(branch, endpoint * double(i) / double(points - 1));
            },
            [&curve = curves[b]](std::size_t, auto point) { curve.push_back(std::move(point)); }, execution));
      }
    }
    catch (...)
    {
      error = std::current_exception();
    }
    // Both curves are submitted before observing the background or either grid.
    // Observe both completions even if one fails: collectors borrow curves.
    for (auto& result : completed)
    {
      try
      {
        if (result) (void)result->get_wait(scheduler);
      }
      catch (...)
      {
        if (!error) error = std::current_exception();
      }
    }
    // Preparation may have been accepted even if the first grid was rejected.
    try
    {
      (void)background.get_wait(scheduler);
    }
    catch (...)
    {
      if (!error) error = std::current_exception();
    }
    if (error) std::rethrow_exception(error);
    std::cout << "# c=4 n=1; one prepared background; energy is the fixed-N gap\n"
              << "# branch momentum energy\n"
              << std::setprecision(17);
    for (std::size_t b = 0; b < curves.size(); ++b)
      for (auto const& point : curves[b])
      {
        if (!point.converged || !point.energy) throw std::runtime_error("unconverged dispersion point");
        std::cout << (b == 0 ? "I" : "II") << ' ' << point.momentum << ' ' << *point.energy << '\n';
      }
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
