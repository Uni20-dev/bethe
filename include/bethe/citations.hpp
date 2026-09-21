// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Generated from data/citations.json by scripts/generate_citations.py. Do not edit.
// clang-format off
#pragma once

#include <array>
#include <span>
#include <stdexcept>
#include <string_view>

namespace bethe::citations
{
struct Link { std::string_view label, url; };
struct Reference
{
  std::string_view id, authors, title, publication;
  int year; // 0 means unspecified, e.g. undated web notes.
  std::span<Link const> links;
};
struct Use { Reference const* reference; std::string_view context; };
inline constexpr std::array<Link, 2> links_0{{
  {"DOI", "https://doi.org/10.1063/1.168740"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9809163"},
}};
inline constexpr std::array<Link, 1> links_1{{
  {"arXiv", "https://arxiv.org/abs/1702.06550"},
}};
inline constexpr std::array<Link, 1> links_2{{
  {"XXX spinons", "https://integrability.org/g_h_e.html"},
}};
inline constexpr std::array<Link, 1> links_3{{
  {"XXZ spinons", "https://integrability.org/g_sc_p_e.html"},
}};
inline constexpr std::array<Link, 1> links_4{{
  {"arXiv", "https://arxiv.org/abs/1609.08045"},
}};
inline constexpr std::array<Link, 3> links_5{{
  {"Real rapidities", "https://integrability.org/c_h_e_rr.html"},
  {"SU(2) descendants", "https://integrability.org/c_h_e_rr_10.html"},
  {"Two-string states", "https://integrability.org/c_h_e_s_2.html"},
}};
inline constexpr std::array<Link, 2> links_6{{
  {"DOI", "https://doi.org/10.21468/SciPostPhys.1.1.008"},
  {"arXiv", "https://arxiv.org/abs/1606.09516"},
}};
inline constexpr std::array<Link, 1> links_7{{
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/0207529v2"},
}};

inline constexpr std::array<Reference, 8> references{{
  {"karbach-1998", "Michael Karbach, Kun Hu, and Gerhard Müller", "Introduction to the Bethe ansatz II", "Computers in Physics 12, 565", 1998, links_0},
  {"groha-2017", "Stefan Groha and Fabian H. L. Essler", "Spinon decay in the spin-1/2 Heisenberg chain with weak next nearest neighbour exchange", "J. Phys. A 50, 334002", 2017, links_1},
  {"caux-xxx-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXX spinons", "Online notes", 0, links_2},
  {"caux-xxz-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXZ spinons", "Online notes", 0, links_3},
  {"mei-2017", "Zhongtao Mei and C. J. Bolech", "Derivation of matrix product states for the Heisenberg spin chain with open boundary conditions", "Phys. Rev. E 95, 032127", 2017, links_4},
  {"caux-xxx-states", "Jean-Sébastien Caux", "The Bethe Ansatz: real rapidities, SU(2) descendants, and two-string states", "Online notes", 0, links_5},
  {"vlijm-2016", "R. Vlijm, I. S. Eliëns, and J.-S. Caux", "Correlations of zero-entropy critical states in the XXZ model: integrability and Luttinger theory far from the ground state", "SciPost Phys. 1, 008", 2016, links_6},
  {"lieb-wu-2003", "Elliott H. Lieb and F. Y. Wu", "The one-dimensional Hubbard model: A reminiscence", "Physica A 321, 1-27", 2003, links_7},
}};

inline constexpr std::array<Use, 4> uses_xxx_pbc{{
  {&references[0], "Periodic XXX equations, energy normalization, and sector quantum numbers; Eqs. (6)-(9), (16), Table I."},
  {&references[1], "Odd-chain one-spinon states and hole labels, Sec. 3; only the unperturbed integrable chain is used."},
  {&references[2], "Thermodynamic XXX spinon dispersion and the bulk reference for --spinons."},
  {&references[5], "Highest-weight and string-classification background; strings and infinite-root descendants are not implemented."},
}};
inline constexpr std::array<Use, 2> uses_xxx_obc{{
  {&references[4], "Free-end XXX equations: rational limit of Eqs. (11)-(12) and footnote 2; our spin-1/2 normalization differs."},
  {&references[5], "SU(2) and string-classification background, not the boundary equations; scans include only the supported finite-real family."},
}};
inline constexpr std::array<Use, 3> uses_xxz_pbc{{
  {&references[6], "Periodic XXZ equations, energy, and momentum, Eqs. (1), (3)-(5); our energy includes the N*Delta/4 shift."},
  {&references[0], "XXX limit at Delta=1 and the conventional real-root quantum-number window."},
  {&references[5], "Background for the XXX window used to restrict the XXZ scan; not a complete XXZ state classification."},
}};
inline constexpr std::array<Use, 1> uses_xxz_obc{{
  {&references[4], "Free-end XXZ equations and boundary reflection phase, Eqs. (11)-(12); our Hamiltonian is divided by four and shifted."},
}};
inline constexpr std::array<Use, 1> uses_hubbard_pbc{{
  {&references[7], "Lieb-Wu ground-state equations and quantum-number parity, Eqs. (1), (11), (14)-(18); t=1 and unshifted U*n_up*n_down."},
}};

enum class Tool { xxx_pbc, xxx_obc, xxz_pbc, xxz_obc, hubbard_pbc };

[[nodiscard]] constexpr std::span<Use const> for_tool(Tool tool)
{
  switch (tool)
  {
    case Tool::xxx_pbc: return uses_xxx_pbc;
    case Tool::xxx_obc: return uses_xxx_obc;
    case Tool::xxz_pbc: return uses_xxz_pbc;
    case Tool::xxz_obc: return uses_xxz_obc;
    case Tool::hubbard_pbc: return uses_hubbard_pbc;
  }
  throw std::invalid_argument("unknown citation tool");
}

[[nodiscard]] constexpr Reference const* find(std::string_view id)
{
  for (auto const& reference : references)
    if (reference.id == id) return &reference;
  return nullptr;
}
} // namespace bethe::citations
// clang-format on
