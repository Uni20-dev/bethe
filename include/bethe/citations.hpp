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
inline constexpr std::array<Link, 2> links_7{{
  {"arXiv v2", "https://arxiv.org/abs/1508.05741v2"},
  {"HTML", "https://arxiv.org/html/1508.05741v2"},
}};
inline constexpr std::array<Link, 1> links_8{{
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/0207529v2"},
}};
inline constexpr std::array<Link, 2> links_9{{
  {"DOI", "https://doi.org/10.1088/1742-5468/ac98be"},
  {"arXiv", "https://arxiv.org/abs/2206.07985"},
}};
inline constexpr std::array<Link, 1> links_10{{
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/9704138v2"},
}};
inline constexpr std::array<Link, 1> links_11{{
  {"DOI", "https://doi.org/10.1103/PhysRev.125.164"},
}};
inline constexpr std::array<Link, 1> links_12{{
  {"DOI", "https://doi.org/10.1103/PhysRev.130.1605"},
}};
inline constexpr std::array<Link, 1> links_13{{
  {"DOI", "https://doi.org/10.1103/PhysRev.130.1616"},
}};
inline constexpr std::array<Link, 1> links_14{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.4.386"},
}};
inline constexpr std::array<Link, 1> links_15{{
  {"DOI", "https://doi.org/10.1016/0375-9601(67)90193-4"},
}};
inline constexpr std::array<Link, 1> links_16{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.19.1312"},
}};
inline constexpr std::array<Link, 1> links_17{{
  {"DOI", "https://doi.org/10.1103/PhysRev.168.1920"},
}};
inline constexpr std::array<Link, 1> links_18{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.12.3795"},
}};
inline constexpr std::array<Link, 2> links_19{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.46.9147"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9207007"},
}};
inline constexpr std::array<Link, 1> links_20{{
  {"DOI", "https://doi.org/10.1016/0375-9601(82)90403-0"},
}};
inline constexpr std::array<Link, 2> links_21{{
  {"DOI", "https://doi.org/10.1088/1742-5468/2014/05/P05009"},
  {"arXiv", "https://arxiv.org/abs/1401.4450"},
}};
inline constexpr std::array<Link, 2> links_22{{
  {"DOI", "https://doi.org/10.1103/RevModPhys.76.643"},
  {"arXiv", "https://arxiv.org/abs/nucl-th/0405011"},
}};
inline constexpr std::array<Link, 2> links_23{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.88.085323"},
  {"arXiv", "https://arxiv.org/abs/1306.2541"},
}};
inline constexpr std::array<Link, 1> links_24{{
  {"arXiv", "https://arxiv.org/abs/1011.0128"},
}};
inline constexpr std::array<Link, 2> links_25{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.73.021602"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0505632"},
}};
inline constexpr std::array<Link, 2> links_26{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.60.9236"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9901168"},
}};
inline constexpr std::array<Link, 1> links_27{{
  {"arXiv", "https://arxiv.org/abs/solv-int/9710002"},
}};
inline constexpr std::array<Link, 1> links_28{{
  {"DOI", "https://doi.org/10.1016/0003-4916(73)90441-7"},
}};
inline constexpr std::array<Link, 1> links_29{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.60.635"},
}};
inline constexpr std::array<Link, 1> links_30{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.60.639"},
}};
inline constexpr std::array<Link, 1> links_31{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.4.2019"},
}};
inline constexpr std::array<Link, 1> links_32{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.45.379"},
}};
inline constexpr std::array<Link, 1> links_33{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.69.2313"},
}};
inline constexpr std::array<Link, 1> links_34{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.46.844"},
}};
inline constexpr std::array<Link, 1> links_35{{
  {"DOI", "https://doi.org/10.1103/PhysRev.150.321"},
}};
inline constexpr std::array<Link, 1> links_36{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.65.243"},
}};
inline constexpr std::array<Link, 2> links_37{{
  {"DOI", "https://doi.org/10.1088/0305-4470/27/18/021"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9306089"},
}};
inline constexpr std::array<Link, 2> links_38{{
  {"arXiv", "https://arxiv.org/abs/2307.12410v1"},
  {"Equations (HTML)", "https://arxiv.org/html/2307.12410v1"},
}};
inline constexpr std::array<Link, 1> links_39{{
  {"DOI", "https://doi.org/10.1063/1.1664947"},
}};
inline constexpr std::array<Link, 2> links_40{{
  {"DOI", "https://doi.org/10.1016/S0550-3213(98)00239-9"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9803118"},
}};
inline constexpr std::array<Link, 2> links_41{{
  {"DOI", "https://doi.org/10.1088/0305-4470/39/5/005"},
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/0511694v2"},
}};
inline constexpr std::array<Link, 2> links_42{{
  {"DOI", "https://doi.org/10.21468/SciPostPhys.7.2.023"},
  {"arXiv v4", "https://arxiv.org/abs/1901.10932v4"},
}};
inline constexpr std::array<Link, 2> links_43{{
  {"DOI", "https://doi.org/10.1088/1742-5468/2015/05/P05037"},
  {"arXiv v2", "https://arxiv.org/abs/1412.8217v2"},
}};

inline constexpr std::array<Reference, 44> references{{
  {"karbach-1998", "Michael Karbach, Kun Hu, and Gerhard Müller", "Introduction to the Bethe ansatz II", "Computers in Physics 12, 565", 1998, links_0},
  {"groha-2017", "Stefan Groha and Fabian H. L. Essler", "Spinon decay in the spin-1/2 Heisenberg chain with weak next nearest neighbour exchange", "J. Phys. A 50, 334002", 2017, links_1},
  {"caux-xxx-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXX spinons", "Online notes", 0, links_2},
  {"caux-xxz-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXZ spinons", "Online notes", 0, links_3},
  {"mei-2017", "Zhongtao Mei and C. J. Bolech", "Derivation of matrix product states for the Heisenberg spin chain with open boundary conditions", "Phys. Rev. E 95, 032127", 2017, links_4},
  {"caux-xxx-states", "Jean-Sébastien Caux", "The Bethe Ansatz: real rapidities, SU(2) descendants, and two-string states", "Online notes", 0, links_5},
  {"vlijm-2016", "R. Vlijm, I. S. Eliëns, and J.-S. Caux", "Correlations of zero-entropy critical states in the XXZ model: integrability and Luttinger theory far from the ground state", "SciPost Phys. 1, 008", 2016, links_6},
  {"kozlowski-2017", "Karol K. Kozlowski", "On condensation properties of Bethe roots associated with the XXZ chain", "arXiv:1508.05741v2", 2017, links_7},
  {"lieb-wu-2003", "Elliott H. Lieb and F. Y. Wu", "The one-dimensional Hubbard model: A reminiscence", "Physica A 321, 1-27", 2003, links_8},
  {"rylands-2022", "Colin Rylands, Bruno Bertini, and Pasquale Calabrese", "Integrable quenches in the Hubbard model", "J. Stat. Mech. 2022, 103103", 2022, links_9},
  {"deguchi-yue-1997", "Tetsuo Deguchi and Ruihong Yue", "Exact solutions of 1-D Hubbard model with open boundary conditions and the conformal dimensions under boundary magnetic fields", "arXiv:cond-mat/9704138", 1997, links_10},
  {"lieb-mattis-1962", "Elliott Lieb and Daniel Mattis", "Theory of Ferromagnetism and the Ordering of Electronic Energy Levels", "Phys. Rev. 125, 164-172", 1962, links_11},
  {"lieb-liniger-1963", "Elliott H. Lieb and Werner Liniger", "Exact Analysis of an Interacting Bose Gas. I. The General Solution and the Ground State", "Phys. Rev. 130, 1605", 1963, links_12},
  {"lieb-1963-excitations", "Elliott H. Lieb", "Exact Analysis of an Interacting Bose Gas. II. The Excitation Spectrum", "Phys. Rev. 130, 1616", 1963, links_13},
  {"gaudin-1971", "M. Gaudin", "Boundary Energy of a Bose Gas in One Dimension", "Phys. Rev. A 4, 386", 1971, links_14},
  {"gaudin-1967", "M. Gaudin", "Un système à une dimension de fermions en interaction", "Phys. Lett. A 24, 55-56", 1967, links_15},
  {"yang-1967", "C. N. Yang", "Some Exact Results for the Many-Body Problem in One Dimension with Repulsive Delta-Function Interaction", "Phys. Rev. Lett. 19, 1312", 1967, links_16},
  {"yang-1968", "C. N. Yang", "S Matrix for the One-Dimensional N-Body Problem with Repulsive or Attractive delta-Function Interaction", "Phys. Rev. 168, 1920", 1968, links_17},
  {"sutherland-1975", "Bill Sutherland", "Model for a multicomponent quantum system", "Phys. Rev. B 12, 3795", 1975, links_18},
  {"essler-korepin-1992", "Fabian H. L. Essler and Vladimir E. Korepin", "Higher conservation laws and algebraic Bethe Ansätze for the supersymmetric t-J model", "Phys. Rev. B 46, 9147", 1992, links_19},
  {"babujian-1982", "H. M. Babujian", "Exact solution of the one-dimensional isotropic Heisenberg chain with arbitrary spins S", "Phys. Lett. A 90, 479-482", 1982, links_20},
  {"vlijm-caux-2014", "Rogier Vlijm and Jean-Sébastien Caux", "Computation of dynamical correlation functions of the spin-1 Babujan-Takhtajan chain", "J. Stat. Mech. 2014, P05009", 2014, links_21},
  {"dukelsky-2004", "J. Dukelsky, S. Pittel, and G. Sierra", "Colloquium: Exactly solvable Richardson-Gaudin models for many-body quantum systems", "Rev. Mod. Phys. 76, 643-662", 2004, links_22},
  {"faribault-schuricht-2013", "Alexandre Faribault and Dirk Schuricht", "Spin decoherence due to a randomly fluctuating spin bath", "Phys. Rev. B 88, 085323", 2013, links_23},
  {"lee-2011", "J. Y. Lee, X. W. Guan, and M. T. Batchelor", "Yang-Yang method for the thermodynamics of one-dimensional multi-component interacting fermions", "J. Phys. A: Math. Theor. 44, 165002", 2011, links_24},
  {"imambekov-demler-2006", "Adilet Imambekov and Eugene Demler", "Exactly solvable case of a one-dimensional Bose-Fermi mixture", "Phys. Rev. A 73, 021602(R)", 2006, links_25},
  {"wang-1999", "Yupeng Wang", "Exact solution of a spin-ladder model", "Phys. Rev. B 60, 9236", 1999, links_26},
  {"bogoliubov-1997", "N. M. Bogoliubov, A. G. Izergin, and N. A. Kitanine", "Correlation functions for a strongly correlated boson system", "arXiv:solv-int/9710002", 1997, links_27},
  {"baxter-1973", "R. J. Baxter", "Eight-vertex model in lattice statistics and one-dimensional anisotropic Heisenberg chain. III. Eigenvectors of the transfer matrix and Hamiltonian", "Ann. Phys. 76, 48-71", 1973, links_28},
  {"haldane-1988", "F. D. M. Haldane", "Exact Jastrow-Gutzwiller resonating-valence-bond ground state of the spin-1/2 antiferromagnetic Heisenberg chain with 1/r^2 exchange", "Phys. Rev. Lett. 60, 635", 1988, links_29},
  {"shastry-1988", "B. Sriram Shastry", "Exact solution of an S=1/2 Heisenberg antiferromagnetic chain with long-ranged interactions", "Phys. Rev. Lett. 60, 639", 1988, links_30},
  {"sutherland-1971", "Bill Sutherland", "Exact Results for a Quantum Many-Body Problem in One Dimension", "Phys. Rev. A 4, 2019", 1971, links_31},
  {"andrei-1980", "N. Andrei", "Diagonalization of the Kondo Hamiltonian", "Phys. Rev. Lett. 45, 379", 1980, links_32},
  {"destri-de-vega-1992", "C. Destri and H. J. de Vega", "New thermodynamic Bethe ansatz equations without strings", "Phys. Rev. Lett. 69, 2313-2317", 1992, links_33},
  {"gwa-spohn-1992", "Leh-Hun Gwa and Herbert Spohn", "Bethe solution for the dynamical-scaling exponent of the noisy Burgers equation", "Phys. Rev. A 46, 844", 1992, links_34},
  {"yang-yang-1966", "C. N. Yang and C. P. Yang", "One-Dimensional Chain of Anisotropic Spin-Spin Interactions. I. Proof of Bethe's Hypothesis for Ground State in a Finite System", "Phys. Rev. 150, 321", 1966, links_35},
  {"shastry-sutherland-1990", "B. Sriram Shastry and Bill Sutherland", "Twisted boundary conditions and effective mass in Heisenberg-Ising and Hubbard rings", "Phys. Rev. Lett. 65, 243", 1990, links_36},
  {"de-vega-gonzalez-ruiz-1994", "H. J. de Vega and A. González-Ruiz", "Boundary K-matrices for the XYZ, XXZ and XXX spin chains", "J. Phys. A: Math. Gen. 27, 6129-6138", 1994, links_37},
  {"essler-de-klerk-2023", "F. H. L. Essler and A. J. J. M. de Klerk", "Statistics of matrix elements of local operators in integrable models", "arXiv:2307.12410v1", 2023, links_38},
  {"yang-yang-1969", "C. N. Yang and C. P. Yang", "Thermodynamics of a One-Dimensional System of Bosons with Repulsive Delta-Function Interaction", "J. Math. Phys. 10, 1115", 1969, links_39},
  {"doikou-nepomechie-1998", "Anastasia Doikou and Rafael I. Nepomechie", "Bulk and Boundary S Matrices for the SU(N) Chain", "Nucl. Phys. B 521, 547-572", 1998, links_40},
  {"oelkers-2006", "N. Oelkers, M. T. Batchelor, M. Bortz, and X.-W. Guan", "Bethe Ansatz study of one-dimensional Bose and Fermi gases with periodic and hard wall boundary conditions", "J. Phys. A 39, 1073-1098", 2006, links_41},
  {"grijalva-2019", "Sebastian Grijalva, Jacopo De Nardis, and Veronique Terras", "Open XXZ chain and boundary modes at zero temperature", "SciPost Phys. 7, 023", 2019, links_42},
  {"dugave-2015", "M. Dugave, F. Göhmann, K. K. Kozlowski, and J. Suzuki", "On form-factor expansions for the XXZ chain in the massive regime", "J. Stat. Mech. (2015) P05037", 2015, links_43},
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
inline constexpr std::array<Use, 4> uses_xxz_pbc{{
  {&references[6], "Massless periodic XXZ equations, energy, and momentum, Eqs. (1), (3)-(5); our energy includes the N*Delta/4 shift."},
  {&references[43], "Massive periodic XXZ Hamiltonian and trigonometric Bethe equations, Eqs. (1.1)-(1.2), with zero twist and field; Pauli exchange is divided by four. Used for Delta>1 sector ground states, not form factors, strings or excitations."},
  {&references[0], "XXX limit at Delta=1 and the conventional real-root quantum-number window."},
  {&references[5], "Background for the XXX window used to restrict the XXZ scan; not a complete XXZ state classification."},
}};
inline constexpr std::array<Use, 2> uses_xxz_obc{{
  {&references[4], "Free-end XXZ equations and boundary reflection phase, Eqs. (11)-(12); our Hamiltonian is divided by four and shifted."},
  {&references[42], "Massive open XXZ boundary root and finite-size deviation, Sec. 4.3.2, specialized to zero boundary fields; used for Delta>1 ground states, not boundary correlations or excitation scans."},
}};
inline constexpr std::array<Use, 2> uses_hubbard_pbc{{
  {&references[8], "Lieb-Wu ground-state equations and quantum-number parity, Eqs. (1), (11), (14)-(18); t=1 and unshifted U*n_up*n_down."},
  {&references[9], "Full particle-hole and partial particle-hole (Shiba) transformations, Sec. II, Eqs. (4)-(7); used for sector mappings, not quench dynamics."},
}};
inline constexpr std::array<Use, 3> uses_hubbard_obc{{
  {&references[10], "Free-end Hubbard energy and nested reflection equations, Eqs. (2.5)-(2.8), (3.1)-(3.2), with all boundary potentials zero; t=1 and unshifted interaction."},
  {&references[11], "Spin ordering for the nearest-neighbour open chain: the sector minimum has S=|Sz|. No periodic-shell spin-branch selection is required."},
  {&references[9], "Bipartite particle-hole and Shiba transformations, Sec. II; used for physical-to-auxiliary sector mappings, not quench dynamics."},
}};
inline constexpr std::array<Use, 3> uses_lieb_liniger_pbc{{
  {&references[12], "Repulsive Bose gas on a ring, H=-sum d_j^2+2c sum delta, and its ground state; c>0."},
  {&references[13], "Excited-state background; we solve finite-volume states in explicit label windows, not a thermodynamic dispersion calculation."},
  {&references[38], "Explicit normalization and finite-ring equations, Eqs. (4), (24)-(27); root-density equation (32)-(33) used for validation, not a thermodynamics API. No matrix elements are implemented."},
}};
inline constexpr std::array<Use, 2> uses_su3_pbc{{
  {&references[18], "Original multicomponent permutation-chain solution; the implementation selects only the fundamental SU(3) periodic balanced ground state."},
  {&references[40], "Nested equations and energy, Eqs. (2.17)-(2.19), logarithmic labels (2.24)-(2.29), and the filled-sea singlet in Sec. 2.3. Our H=sum P gives E=2*E_paper+L. No strings, S matrices or boundary fields are implemented."},
}};
inline constexpr std::array<Use, 3> uses_gaudin_yang_pbc{{
  {&references[15], "Original spin-1/2 continuum fermion solution; we implement repulsive periodic ground states in selected sectors, not attraction."},
  {&references[16], "Original multicomponent delta-gas solution; our implementation has two spin components only."},
  {&references[41], "Hamiltonian and periodic fermion equations (1), (2), (7), (8), (25); even N, odd minority population as in Sec. 5. Weak and strong limits (12), (15)-(16) provide checks. Hard walls, attraction and excitations are not implemented."},
}};

enum class Tool { xxx_pbc, xxx_obc, xxz_pbc, xxz_obc, hubbard_pbc, hubbard_obc, lieb_liniger_pbc, su3_pbc, gaudin_yang_pbc };

[[nodiscard]] constexpr std::span<Use const> for_tool(Tool tool)
{
  switch (tool)
  {
    case Tool::xxx_pbc: return uses_xxx_pbc;
    case Tool::xxx_obc: return uses_xxx_obc;
    case Tool::xxz_pbc: return uses_xxz_pbc;
    case Tool::xxz_obc: return uses_xxz_obc;
    case Tool::hubbard_pbc: return uses_hubbard_pbc;
    case Tool::hubbard_obc: return uses_hubbard_obc;
    case Tool::lieb_liniger_pbc: return uses_lieb_liniger_pbc;
    case Tool::su3_pbc: return uses_su3_pbc;
    case Tool::gaudin_yang_pbc: return uses_gaudin_yang_pbc;
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
