// SPDX-License-Identifier: GPL-3.0-or-later
#include <bethe/central_spin.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/richardson.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>
#include <uni20/core/scalar_io.hpp>

// Independent parameter points share one scheduler. Each point's continuation
// remains sequential, with its own numerical workspace and explicit precision.
template <uni20::Real Real>
int sweep(bethe::richardson::SolverOptions<Real> pairing, bethe::central_spin::SolverOptions<Real> central)
{
  uni20::async::TbbScheduler scheduler(4);
  auto v = [p = pairing.precision()](int n) { return uni20::scalar_at_precision<Real>(n, p); };
  using Result = std::pair<bethe::richardson::State<Real>, bethe::central_spin::State<Real>>;
  std::vector<Result> results;
  auto complete = [&] {
    std::vector<Real> levels{v(0), v(1), v(2), v(3)};
    std::vector<Real> couplings{v(3) / v(10), v(7) / v(10), v(1)};
    return bethe::evaluate_observable_grid_async(
        6,
        [levels = std::move(levels), couplings = std::move(couplings), pairing, central](std::size_t i) {
          Real const parameter = uni20::scalar_at_precision<Real>(i + 1, pairing.precision()) / Real{4};
          // The span-based APIs borrow these owned vectors only during each solve.
          return Result{
              bethe::richardson::ground_state<Real>(levels, 2, parameter, {}, pairing),
              bethe::central_spin::sector_ground_state<Real>(couplings, parameter, uni20::half_int{0}, central)};
        },
        [&results](std::size_t, Result state) { results.push_back(std::move(state)); },
        {.scheduler = &scheduler, .max_in_flight = 3});
  }();

  // Completion includes all admitted work, also when a solve throws. Numerical
  // nonconvergence is carried by each state. Print only after ordered collection.
  (void)complete.get_wait(scheduler);
  std::cout << "# Richardson/central-spin parameter sweep; working bits=" << pairing.precision().bit_count()
            << "\n# coupling_or_field pairing_energy central_spin_energy\n";
  for (auto const& [a, b] : results)
  {
    if (!a.converged || !b.converged)
      throw std::runtime_error("ground state did not converge within its numerical budget");
    std::cout << uni20::format_real(a.coupling) << ' ' << uni20::format_real(a.energy) << ' '
              << uni20::format_real(*b.energy) << '\n';
  }
  return 0;
}

int main()
{
  try
  {
#if UNI20_HAS_MPREAL
    auto const p = uni20::Precision::bits(160);
    return sweep(bethe::richardson::SolverOptions<uni20::mpreal>(p),
                 bethe::central_spin::SolverOptions<uni20::mpreal>(p));
#else
    return sweep(bethe::richardson::SolverOptions<double>{}, bethe::central_spin::SolverOptions<double>{});
#endif
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
