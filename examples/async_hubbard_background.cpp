// SPDX-License-Identifier: GPL-3.0-or-later
#include <bethe/hubbard_doped.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

#include <array>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <vector>

int main()
{
  namespace model = bethe::hubbard::thermo;
  try
  {
    uni20::async::TbbScheduler scheduler(4);
    bethe::ExecutionOptions execution{.scheduler = &scheduler, .max_in_flight = 4};
    auto background = bethe::prepare_background_async([] { return model::DopedSolver<double>(4.0, 0.5); }, execution);
    constexpr std::size_t points = 17;
    std::array const branches{model::DopedBranch::spinon, model::DopedBranch::holon,
                              model::DopedBranch::charge_particle};
    std::array const names{"spinon", "holon", "charge-particle"};
    std::array<std::vector<model::DopedPoint<double>>, 3> curves;
    std::array<std::optional<uni20::async::Async<std::size_t>>, 3> completed;
    std::exception_ptr error;
    try
    {
      for (std::size_t b = 0; b < curves.size(); ++b)
        completed[b].emplace(bethe::evaluate_observable_grid_async(
            background.read(), points,
            [branch = branches[b]](model::DopedSolver<double> const& solver, std::size_t i) {
              auto const [low, high] = solver.momentum_range(branch);
              auto const p = i == points - 1 ? high : low + (high - low) * double(i) / double(points - 1);
              return solver.at_momentum(branch, p, model::Convention::symmetric, model::EnergyReference::fermi);
            },
            [&curve = curves[b]](std::size_t, auto point) { curve.push_back(std::move(point)); }, execution));
    }
    catch (...)
    {
      error = std::current_exception();
    }
    // All three grids depend on the same prepared meshes. Wait on each final
    // result, even after a failure: their collectors borrow the curve vectors.
    for (auto& result : completed)
      try
      {
        if (result) (void)result->get_wait(scheduler);
      }
      catch (...)
      {
        if (!error) error = std::current_exception();
      }
    // Preparation can be accepted even if the first grid submission fails.
    try
    {
      if (!background.get_wait(scheduler).background().converged)
        throw std::runtime_error("unconverged Hubbard background");
    }
    catch (...)
    {
      if (!error) error = std::current_exception();
    }
    if (error) std::rethrow_exception(error);
    std::cout << "# U=4 n=0.5 t=1; one prepared background; energies relative to the Fermi level\n"
              << "# branch momentum energy\n"
              << std::setprecision(17);
    for (std::size_t b = 0; b < curves.size(); ++b)
      for (auto const& point : curves[b])
      {
        if (!point.converged || !point.energy) throw std::runtime_error("unconverged dispersion point");
        std::cout << names[b] << ' ' << point.momentum << ' ' << *point.energy << '\n';
      }
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
