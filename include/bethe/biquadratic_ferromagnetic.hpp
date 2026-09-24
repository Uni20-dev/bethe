// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/biquadratic_qsystem.hpp>
#include <bethe/xxz_open_three_string.hpp>
#include <bethe/xxz_open_two_string.hpp>

namespace bethe::biquadratic::ferromagnetic
{
namespace detail
{
inline void validate_sites(std::size_t sites)
{
  if (sites < 2) throw std::invalid_argument("ferromagnetic free-end chain requires N>=2");
}

// Apply only to states freshly produced by the AF API. Auxiliary XXZ energies,
// roots and convergence diagnostics retain their original convention.
template <typename State> State reverse(State state)
{
  if constexpr (requires { *state.energy; })
  {
    if (state.energy) state.energy = -*state.energy;
    if (state.tl_energy) state.tl_energy = -*state.tl_energy;
  }
  else
  {
    state.energy = -state.energy;
    state.tl_energy = -state.tl_energy;
  }
  state.exchange = Exchange::ferromagnetic;
  return state;
}
} // namespace detail

/// H=+sum (S.S)^2=(N-1)+sum e_i. All e_i annihilate the ground space.
/// Works for odd and even N; multiplicity is F_(2N+2), not one SU(2) multiplet.
template <uni20::Real Real> struct GroundSpace
{
    std::size_t sites{}, through_lines{};
    Real energy{};
    std::optional<std::uint64_t> multiplicity;
};

template <uni20::Real Real = double> [[nodiscard]] GroundSpace<Real> ground_space(std::size_t sites)
{
  detail::validate_sites(sites);
  return {.sites = sites,
          .through_lines = sites,
          .energy = Real(sites - 1),
          .multiplicity = temperley_lieb::spin_chain_multiplicity(3, sites)};
}

/// Gap above the ENTIRE degenerate ground space, not between its zero modes.
/// Koma-Nachtergaele Prop. 2, rescaled by 2*Delta=3 via TL equivalence.
template <uni20::Real Real = double> [[nodiscard]] Real spectral_gap(std::size_t sites)
{
  detail::validate_sites(sites);
  Real const s = std::sin((Real{2} * std::atan(Real{1})) / Real(sites));
  return Real{1} + Real{4} * s * s; // 3-2*cos(pi/N)
}

template <uni20::Real Real> struct OneDefectLevel
{
    std::size_t sites{}, through_lines{}, mode{};
    /// OBC standing-wave coordinate pi*mode/N, NOT a translation eigenvalue.
    Real wave_number{}, energy{}, gap{};
    /// Copies per TL eigenvector; accidental coincidences add other modules.
    std::optional<std::uint64_t> multiplicity;
};

/// Complete ell=N-2 module: j=1,...,N-1, gap=3+2*cos(pi*j/N).
/// Descending j gives ascending energy. Exact formula, rounded at native precision.
/// One TL singlet defect is NOT one physical spin lowering (the latter costs zero).
template <uni20::Real Real = double>
[[nodiscard]] OneDefectLevel<Real> one_defect_level(std::size_t sites, std::size_t mode)
{
  detail::validate_sites(sites);
  if (mode == 0 || mode >= sites) throw std::invalid_argument("one-defect mode must satisfy 1<=j<N");
  Real const pi = Real{4} * std::atan(Real{1});
  // Compute the gap independently of the extensive energy, including near j=N.
  Real const s = std::sin((pi / Real{2}) * (Real(sites - mode) / Real(sites)));
  Real const gap = Real{1} + Real{4} * s * s;
  return {.sites = sites,
          .through_lines = sites - 2,
          .mode = mode,
          .wave_number = pi * (Real(mode) / Real(sites)),
          .energy = Real(sites - 1) + gap,
          .gap = gap,
          .multiplicity = temperley_lieb::spin_chain_multiplicity(3, sites - 2)};
}

/// Selected odd/even-chain Bethe eigenstate; no lowest-state claim.
/// tl_energy=+sum e_i=E-E0; publish a gap only if reference.converged.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve_real(std::size_t sites, std::span<uni20::half_int const> numbers,
                                     SolverOptions<Real> const& options = {})
{
  return detail::reverse(biquadratic::solve_real<Real>(sites, numbers, options));
}

namespace detail
{
// Shared scan for the complete positive-real family or a selected high-label
// window. Bethe equations and the bounded combination heap remain model-independent.
template <uni20::Real Real = double>
RealExcitationScan<Real> real_scan(std::size_t sites, std::size_t through_lines, std::optional<std::size_t> window,
                                   RealExcitationOptions const& scan, SolverOptions<Real> const& solver)
{
  namespace qg = xxz::quantum_group;
  auto const m = qg::detail::real_sector_roots(sites, through_lines);
  auto slots = sites - m;
  std::int64_t first = 2;
  if (window)
  {
    if (m == 0 || *window < m || *window > slots)
      throw std::invalid_argument("real-root label window requires 1<=M<=width<=N-M");
    first = 2 * std::int64_t(slots - *window + 1);
    slots = *window;
  }
  // A window of exactly M labels has only one candidate, even for enormous
  // M. Validate its Newton dimensions before the scanner allocates labels.
  qg::detail::check_matrix_size<Real>(m);
  Real const delta = Real{3} / Real{2};
  auto input = bethe::detail::scan_real_combinations<bethe::RealExcitationScan<qg::RealState<Real>>>(
      slots, m, first, scan, [&](auto const& numbers) { return qg::solve_real<Real>(sites, delta, numbers, solver); },
      [&] { return qg::solve_real<Real>(sites, delta, {}, solver); }, bethe::detail::EnergyOrder::descending,
      [](auto const& state) { return state.energy_shift; });
  return temperley_lieb::detail::map_real_scan(
      std::move(input),
      [](auto state) {
        return detail::reverse(
            biquadratic::detail::from_tl(temperley_lieb::detail::from_reference(std::move(state), Real{3})));
      },
      -Real{2});
}
} // namespace detail

/// Restricted finite-real family, ordered in FERROMAGNETIC energy before truncation.
/// Odd/even N; excludes complex-root module minima. The vacuum is the exact
/// global ground reference, even when the numerical iteration budget is zero.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<Real> real_excitations(std::size_t sites, std::size_t through_lines,
                                                        RealExcitationOptions const& scan = {},
                                                        SolverOptions<Real> const& solver = {})
{
  return detail::real_scan<Real>(sites, through_lines, std::nullopt, scan, solver);
}

/// Scan M-tuples in the HIGHEST width integer labels, I=N-M-width+1,...,N-M.
/// A low-energy scattering window, NOT a complete sector or globally ranked list.
/// Work scales with choose(width,M), not choose(N-M,M). Gaps and ordering use
/// direct reference energy shifts, independently of the extensive total energy.
template <uni20::Real Real = double>
[[nodiscard]] RealExcitationScan<Real>
real_excitations_window(std::size_t sites, std::size_t through_lines, std::size_t width,
                        RealExcitationOptions const& scan = {}, SolverOptions<Real> const& solver = {})
{
  return detail::real_scan<Real>(sites, through_lines, width, scan, solver);
}

template <uni20::Real Real, typename Reference> struct BoundClusterState
{
    std::size_t sites{}, through_lines{}, mode{};
    Real energy{}, tl_energy{}; // tl_energy=E-E0; an estimate unless reference.converged.
    std::optional<std::uint64_t> multiplicity;
    Reference reference;
};
template <uni20::Real Real> using BoundPairState = BoundClusterState<Real, xxz::quantum_group::two_string::State<Real>>;
template <uni20::Real Real>
using BoundTripleState = BoundClusterState<Real, xxz::quantum_group::three_string::State<Real>>;

namespace detail
{
template <uni20::Real Real, typename Reference>
BoundClusterState<Real, Reference> from_cluster(Reference ref, std::size_t defects, std::size_t mode)
{
  auto const sites = ref.sites;
  Real const gap = -Real{2} * ref.energy_shift;
  return {.sites = sites,
          .through_lines = sites - 2 * defects,
          .mode = mode,
          .energy = Real(sites - 1) + gap,
          .tl_energy = gap,
          .multiplicity = temperley_lieb::spin_chain_multiplicity(3, sites - 2 * defects),
          .reference = std::move(ref)};
}
} // namespace detail

/// One targeted two-singlet bound-pair mode, ell=N-4; odd/even N>=4.
/// mode=1,...,N-3 follows the two-string family from its low-energy edge.
/// Does not enumerate scattering levels, other modules, or physical SU(2) spins.
template <uni20::Real Real = double>
[[nodiscard]] BoundPairState<Real> bound_pair(std::size_t sites, std::size_t mode = 1,
                                              SolverOptions<Real> const& options = {})
{
  return detail::from_cluster<Real>(
      xxz::quantum_group::two_string::bound_pair<Real>(sites, Real{3} / Real{2}, mode, options), 2, mode);
}

/// Selected three-singlet droplet family; ell=N-6, mode=1,...,N-5, odd/even.
/// No scattering spectrum or full-sector ranking is implied.
template <uni20::Real Real = double>
[[nodiscard]] BoundTripleState<Real> bound_triple(std::size_t sites, std::size_t mode = 1,
                                                  SolverOptions<Real> const& options = {})
{
  return detail::from_cluster<Real>(
      xxz::quantum_group::three_string::bound_triple<Real>(sites, Real{3} / Real{2}, mode, options), 3, mode);
}

namespace qsystem
{
template <uni20::Real Real = double>
[[nodiscard]] biquadratic::qsystem::State<Real> solve(std::size_t sites, std::span<Real const> seed,
                                                      SolverOptions<Real> const& options = {})
{
  return detail::reverse(biquadratic::qsystem::solve<Real>(sites, seed, options));
}

/// Sign reversal preserves discovery and completeness diagnostics. Sorting is by
/// physical ferro energy; incomplete searches are NOT guaranteed lowest levels.
template <uni20::Real Real = double>
[[nodiscard]] biquadratic::qsystem::Spectrum<Real>
spectrum(std::size_t sites, std::size_t through_lines, xxz::quantum_group::qsystem::SearchOptions const& search = {},
         SolverOptions<Real> const& options = {})
{
  auto result = biquadratic::qsystem::spectrum<Real>(sites, through_lines, search, options);
  for (auto& state : result.states)
    state = detail::reverse(std::move(state));
  std::reverse(result.states.begin(), result.states.end());
  return result;
}
} // namespace qsystem
} // namespace bethe::biquadratic::ferromagnetic
