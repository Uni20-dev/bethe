// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/extrema.hpp>
#include <bethe/hubbard_thermo.hpp>

namespace bethe::hubbard::thermo
{
enum class ChargeChannel
{
  spinon_holon,
  spinon_antiholon,
  holon_antiholon
};
template <uni20::Real Real> struct ChargeContinuumOptions
{
    // Search acts on E/max(1,U); coordinate tolerance is in momentum units.
    bethe::detail::ExtremaOptions<Real> search;
    Options<Real> constituent;
    std::size_t max_quadrature_evaluations = 200000000;
};
#if UNI20_HAS_MPREAL
template <> struct ChargeContinuumOptions<uni20::mpreal>
{
    explicit ChargeContinuumOptions(uni20::Precision p) : search(p), constituent(p) {}
    bethe::detail::ExtremaOptions<uni20::mpreal> search;
    Options<uni20::mpreal> constituent;
    std::size_t max_quadrature_evaluations = 200000000;
};
#endif
template <uni20::Real Real> struct ChargeContinuum
{
    ChargeChannel channel = ChargeChannel::spinon_holon;
    Convention convention = Convention::symmetric;
    Real interaction{}, momentum{}, energy_scale{}, energy_offset{};
    int delta_particles = 0;
    uni20::half_int spin{0};
    std::optional<Real> lower_energy, upper_energy, lower_symmetric_energy, upper_symmetric_energy;
    std::array<Real, 2> lower_momenta{}, upper_momenta{};
    // Vertical estimates exclude horizontal momentum-inversion error.
    Real lower_error{};
    Real upper_error{};
    Real lower_momentum_error{};
    Real upper_momentum_error{};
    std::size_t quadrature_evaluations = 0, constituent_iterations = 0;
    std::size_t objective_evaluations = 0, search_iterations = 0, intervals = 0, meshes = 0;
    bool converged = false;
    bethe::detail::ExtremaStatus search_status = bethe::detail::ExtremaStatus::mesh_limit;
    Status constituent_status = Status::converged;
};

/// Numerically resolved family edges, not certified global extrema or the
/// minimum over arbitrary multiparticle states. Half filling, U>0, t=1, h=0.
template <uni20::Real Real = double>
[[nodiscard]] ChargeContinuum<Real> charge_continuum(ChargeChannel channel, Real interaction, Real momentum,
                                                     Convention convention = Convention::symmetric,
                                                     ChargeContinuumOptions<Real> const& options = {})
{
  detail::validate(interaction, options.constituent, convention);
  bethe::detail::require_same_precision(momentum, interaction);
  Real const zero = uni20::scalar_like(interaction, 0);
  Real const infinity = uni20::scalar_like(interaction, std::numeric_limits<double>::infinity());
  Real const pi = detail::pi(interaction);
  if (!uni20::isfinite(momentum) || uni20::math::abs(momentum) > pi)
    throw std::invalid_argument("charge-continuum total momentum must be in [-pi,pi]");
  if (channel != ChargeChannel::spinon_holon && channel != ChargeChannel::spinon_antiholon &&
      channel != ChargeChannel::holon_antiholon)
    throw std::invalid_argument("invalid Hubbard charge-continuum channel");
  ChargeContinuum<Real> out;
  out.lower_error = out.upper_error = out.lower_momentum_error = out.upper_momentum_error = infinity;
  out.lower_momenta.fill(zero);
  out.upper_momenta.fill(zero);
  if (options.search.precision() != options.constituent.precision())
    throw std::invalid_argument("Hubbard continuum search and constituents require matching precision");
  out.channel = channel;
  out.convention = convention;
  out.interaction = interaction;
  out.momentum = momentum;
  out.energy_scale = std::max(uni20::scalar_like(interaction, 1), interaction);
  bool const mixed = channel != ChargeChannel::holon_antiholon;
  out.delta_particles = mixed ? (channel == ChargeChannel::spinon_holon ? -1 : 1) : 0;
  out.spin = mixed ? uni20::from_twice(1) : uni20::half_int{0};
  out.energy_offset = convention == Convention::unshifted ? interaction / Real{2} * Real(out.delta_particles) : zero;
  Branch const first = mixed ? Branch::spinon : Branch::holon;
  Branch const second = channel == ChargeChannel::spinon_holon ? Branch::holon : Branch::antiholon;
  std::map<std::pair<Branch, Real>, Point<Real>> cache;
  auto point = [&](Branch branch, Real p) -> Point<Real> const* {
    // Reuse exact charge-doublet and spin-reflection symmetries in the cache.
    if (branch == Branch::antiholon)
    {
      branch = Branch::holon;
      p = detail::wrap(p + pi);
    }
    if (branch == Branch::spinon) p = std::min(p, pi - p);
    auto const key = std::pair{branch, p};
    if (auto it = cache.find(key); it != cache.end()) return &it->second;
    auto controls = options.constituent;
    controls.max_evaluations =
        std::min(controls.max_evaluations, options.max_quadrature_evaluations - out.quadrature_evaluations);
    auto state = dispersion(branch, interaction, p, Convention::symmetric, controls);
    out.quadrature_evaluations += state.evaluations;
    out.constituent_iterations += state.iterations;
    if (!state.converged)
    {
      out.constituent_status = state.status;
      return nullptr;
    }
    return &cache.emplace(key, std::move(state)).first->second;
  };
  auto objective = [&](Real p) -> std::optional<bethe::detail::ObjectiveSample<Real>> {
    auto const* a = point(first, p);
    if (!a) return {};
    auto const* b = point(second, detail::wrap(momentum - p));
    if (!b) return {};
    // Normalize before addition to avoid an overflowing intermediate energy.
    Real const value = *a->energy / out.energy_scale + *b->energy / out.energy_scale;
    Real const error = a->energy_error / out.energy_scale + b->energy_error / out.energy_scale +
                       Real{8} * uni20::numeric_limits<Real>::epsilon(interaction) * uni20::math::abs(value);
    return bethe::detail::ObjectiveSample<Real>{value, error};
  };
  auto const search = bethe::detail::bounded_extrema(objective, mixed ? zero : -pi, pi, options.search);
  out.search_status = search.status;
  out.objective_evaluations = search.evaluations;
  out.search_iterations = search.iterations;
  out.intervals = search.intervals;
  out.meshes = search.meshes;
  if (!search.converged) return out;
  Real const lower = search.minimum->value * out.energy_scale, upper = search.maximum->value * out.energy_scale;
  if (!uni20::isfinite(lower) || !uni20::isfinite(upper) || !uni20::isfinite(lower + out.energy_offset) ||
      !uni20::isfinite(upper + out.energy_offset))
  {
    out.search_status = bethe::detail::ExtremaStatus::precision_limit;
    return out;
  }
  out.lower_momenta = {search.minimum->position, detail::wrap(momentum - search.minimum->position)};
  out.upper_momenta = {search.maximum->position, detail::wrap(momentum - search.maximum->position)};
  auto horizontal = [&](std::array<Real, 2> const& ps) {
    // Selected witnesses are already cached; no unbudgeted recomputation.
    return point(first, ps[0])->momentum_error + point(second, ps[1])->momentum_error;
  };
  out.lower_momentum_error = horizontal(out.lower_momenta);
  out.upper_momentum_error = horizontal(out.upper_momenta);
  out.lower_symmetric_energy = lower;
  out.upper_symmetric_energy = upper;
  out.lower_energy = lower + out.energy_offset;
  out.upper_energy = upper + out.energy_offset;
  Real const eps = uni20::numeric_limits<Real>::epsilon(interaction);
  out.lower_error = search.minimum->error * out.energy_scale +
                    Real{4} * eps * (uni20::math::abs(lower) + uni20::math::abs(out.energy_offset));
  out.upper_error = search.maximum->error * out.energy_scale +
                    Real{4} * eps * (uni20::math::abs(upper) + uni20::math::abs(out.energy_offset));
  out.converged = true;
  return out;
}
} // namespace bethe::hubbard::thermo
