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
inline constexpr std::array<Link, 1> links_0{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.40.4621"},
}};
inline constexpr std::array<Link, 1> links_1{{
  {"arXiv", "https://arxiv.org/abs/cond-mat/0012439"},
}};
inline constexpr std::array<Link, 1> links_2{{
  {"arXiv", "https://arxiv.org/abs/1003.1932"},
}};
inline constexpr std::array<Link, 2> links_3{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.59.1734"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9808018"},
}};
inline constexpr std::array<Link, 2> links_4{{
  {"DOI", "https://doi.org/10.1016/0550-3213(95)00105-2"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9410043"},
}};
inline constexpr std::array<Link, 2> links_5{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.81.205120"},
  {"arXiv", "https://arxiv.org/abs/1002.1671"},
}};
inline constexpr std::array<Link, 2> links_6{{
  {"DOI", "https://doi.org/10.1088/1361-6633/ad7b70"},
  {"arXiv", "https://arxiv.org/abs/2307.00890"},
}};
inline constexpr std::array<Link, 2> links_7{{
  {"DOI", "https://doi.org/10.1063/1.168740"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9809163"},
}};
inline constexpr std::array<Link, 1> links_8{{
  {"arXiv", "https://arxiv.org/abs/1702.06550"},
}};
inline constexpr std::array<Link, 1> links_9{{
  {"XXX spinons", "https://integrability.org/g_h_e.html"},
}};
inline constexpr std::array<Link, 1> links_10{{
  {"XXZ spinons", "https://integrability.org/g_sc_p_e.html"},
}};
inline constexpr std::array<Link, 1> links_11{{
  {"arXiv", "https://arxiv.org/abs/1609.08045"},
}};
inline constexpr std::array<Link, 3> links_12{{
  {"Real rapidities", "https://integrability.org/c_h_e_rr.html"},
  {"SU(2) descendants", "https://integrability.org/c_h_e_rr_10.html"},
  {"Two-string states", "https://integrability.org/c_h_e_s_2.html"},
}};
inline constexpr std::array<Link, 2> links_13{{
  {"DOI", "https://doi.org/10.21468/SciPostPhys.1.1.008"},
  {"arXiv", "https://arxiv.org/abs/1606.09516"},
}};
inline constexpr std::array<Link, 2> links_14{{
  {"arXiv v2", "https://arxiv.org/abs/1910.07805v2"},
  {"HTML", "https://arxiv.org/html/1910.07805v2"},
}};
inline constexpr std::array<Link, 3> links_15{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.104.L081410"},
  {"arXiv v2", "https://arxiv.org/abs/2102.03295v2"},
  {"HTML", "https://arxiv.org/html/2102.03295v2"},
}};
inline constexpr std::array<Link, 2> links_16{{
  {"arXiv v2", "https://arxiv.org/abs/1508.05741v2"},
  {"HTML", "https://arxiv.org/html/1508.05741v2"},
}};
inline constexpr std::array<Link, 1> links_17{{
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/0207529v2"},
}};
inline constexpr std::array<Link, 2> links_18{{
  {"DOI", "https://doi.org/10.1088/1742-5468/ac98be"},
  {"arXiv", "https://arxiv.org/abs/2206.07985"},
}};
inline constexpr std::array<Link, 1> links_19{{
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/9704138v2"},
}};
inline constexpr std::array<Link, 1> links_20{{
  {"DOI", "https://doi.org/10.1103/PhysRev.125.164"},
}};
inline constexpr std::array<Link, 1> links_21{{
  {"DOI", "https://doi.org/10.1103/PhysRev.130.1605"},
}};
inline constexpr std::array<Link, 1> links_22{{
  {"DOI", "https://doi.org/10.1103/PhysRev.130.1616"},
}};
inline constexpr std::array<Link, 1> links_23{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.4.386"},
}};
inline constexpr std::array<Link, 1> links_24{{
  {"DOI", "https://doi.org/10.1016/0375-9601(67)90193-4"},
}};
inline constexpr std::array<Link, 1> links_25{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.19.1312"},
}};
inline constexpr std::array<Link, 1> links_26{{
  {"DOI", "https://doi.org/10.1103/PhysRev.168.1920"},
}};
inline constexpr std::array<Link, 1> links_27{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.12.3795"},
}};
inline constexpr std::array<Link, 2> links_28{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.46.9147"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9207007"},
}};
inline constexpr std::array<Link, 1> links_29{{
  {"DOI", "https://doi.org/10.1016/0375-9601(82)90403-0"},
}};
inline constexpr std::array<Link, 2> links_30{{
  {"DOI", "https://doi.org/10.1088/1742-5468/2014/05/P05009"},
  {"arXiv", "https://arxiv.org/abs/1401.4450"},
}};
inline constexpr std::array<Link, 2> links_31{{
  {"DOI", "https://doi.org/10.1103/RevModPhys.76.643"},
  {"arXiv", "https://arxiv.org/abs/nucl-th/0405011"},
}};
inline constexpr std::array<Link, 2> links_32{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.83.235124"},
  {"arXiv", "https://arxiv.org/abs/1103.0472"},
}};
inline constexpr std::array<Link, 2> links_33{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.88.085323"},
  {"arXiv", "https://arxiv.org/abs/1306.2541"},
}};
inline constexpr std::array<Link, 1> links_34{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.20.98"},
}};
inline constexpr std::array<Link, 1> links_35{{
  {"arXiv", "https://arxiv.org/abs/1011.0128"},
}};
inline constexpr std::array<Link, 2> links_36{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.73.021602"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0505632"},
}};
inline constexpr std::array<Link, 2> links_37{{
  {"DOI", "https://doi.org/10.1016/j.aop.2005.11.017"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0510801"},
}};
inline constexpr std::array<Link, 2> links_38{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.60.9236"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/9901168"},
}};
inline constexpr std::array<Link, 2> links_39{{
  {"DOI", "https://doi.org/10.1016/j.nuclphysb.2004.07.032"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0403587"},
}};
inline constexpr std::array<Link, 1> links_40{{
  {"arXiv", "https://arxiv.org/abs/solv-int/9710002"},
}};
inline constexpr std::array<Link, 2> links_41{{
  {"DOI", "https://doi.org/10.1103/PhysRevB.109.115411"},
  {"arXiv", "https://arxiv.org/abs/2312.00161"},
}};
inline constexpr std::array<Link, 1> links_42{{
  {"DOI", "https://doi.org/10.1016/0003-4916(73)90441-7"},
}};
inline constexpr std::array<Link, 1> links_43{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.60.635"},
}};
inline constexpr std::array<Link, 1> links_44{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.60.639"},
}};
inline constexpr std::array<Link, 1> links_45{{
  {"arXiv", "https://arxiv.org/abs/2606.20168v2"},
}};
inline constexpr std::array<Link, 1> links_46{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.4.2019"},
}};
inline constexpr std::array<Link, 1> links_47{{
  {"arXiv", "https://arxiv.org/abs/hep-th/9908127v2"},
}};
inline constexpr std::array<Link, 1> links_48{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.45.379"},
}};
inline constexpr std::array<Link, 2> links_49{{
  {"DOI", "https://doi.org/10.1016/j.nuclphysb.2026.117385"},
  {"arXiv", "https://arxiv.org/abs/2510.25344"},
}};
inline constexpr std::array<Link, 2> links_50{{
  {"DOI", "https://doi.org/10.1088/1742-5468/abb018"},
  {"arXiv", "https://arxiv.org/abs/2007.06489"},
}};
inline constexpr std::array<Link, 1> links_51{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.69.2313"},
}};
inline constexpr std::array<Link, 1> links_52{{
  {"DOI", "https://doi.org/10.1103/PhysRevA.46.844"},
}};
inline constexpr std::array<Link, 1> links_53{{
  {"DOI", "https://doi.org/10.1103/PhysRev.150.321"},
}};
inline constexpr std::array<Link, 1> links_54{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.65.243"},
}};
inline constexpr std::array<Link, 2> links_55{{
  {"DOI", "https://doi.org/10.1088/0305-4470/27/18/021"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9306089"},
}};
inline constexpr std::array<Link, 2> links_56{{
  {"arXiv", "https://arxiv.org/abs/2307.12410v1"},
  {"Equations (HTML)", "https://arxiv.org/html/2307.12410v1"},
}};
inline constexpr std::array<Link, 1> links_57{{
  {"DOI", "https://doi.org/10.1063/1.1664947"},
}};
inline constexpr std::array<Link, 2> links_58{{
  {"DOI", "https://doi.org/10.1016/S0550-3213(98)00239-9"},
  {"arXiv", "https://arxiv.org/abs/hep-th/9803118"},
}};
inline constexpr std::array<Link, 2> links_59{{
  {"DOI", "https://doi.org/10.1088/0305-4470/39/5/005"},
  {"arXiv v2", "https://arxiv.org/abs/cond-mat/0511694v2"},
}};
inline constexpr std::array<Link, 2> links_60{{
  {"DOI", "https://doi.org/10.21468/SciPostPhys.7.2.023"},
  {"arXiv v4", "https://arxiv.org/abs/1901.10932v4"},
}};
inline constexpr std::array<Link, 2> links_61{{
  {"DOI", "https://doi.org/10.1088/1742-5468/2015/05/P05037"},
  {"arXiv v2", "https://arxiv.org/abs/1412.8217v2"},
}};
inline constexpr std::array<Link, 2> links_62{{
  {"General M", "https://integrability.org/c_h_s_m.html"},
  {"Two-magnon scattering", "https://integrability.org/c_h_s_2.html"},
}};
inline constexpr std::array<Link, 3> links_63{{
  {"arXiv", "https://arxiv.org/abs/0911.1881"},
  {"HTML", "https://arxiv.org/html/0911.1881v1"},
  {"Original 1982 paper", "https://doi.org/10.1007/BF01212176"},
}};
inline constexpr std::array<Link, 1> links_64{{
  {"Author manuscript", "https://web.dm.unipi.it/robol/assets/pdf/secular-paper.pdf"},
}};
inline constexpr std::array<Link, 2> links_65{{
  {"DOI", "https://doi.org/10.1553/etna_vol55s401"},
  {"Open-access article", "https://etna.ricam.oeaw.ac.at/vol.55.2022/pp401-423.dir/pp401-423.pdf"},
}};
inline constexpr std::array<Link, 1> links_66{{
  {"arXiv", "https://arxiv.org/abs/cond-mat/9512120"},
}};
inline constexpr std::array<Link, 2> links_67{{
  {"DOI", "https://doi.org/10.1088/1751-8121/ae05d9"},
  {"arXiv", "https://arxiv.org/abs/2302.13126"},
}};
inline constexpr std::array<Link, 2> links_68{{
  {"DOI", "https://doi.org/10.1007/s00023-006-0304-6"},
  {"arXiv", "https://arxiv.org/abs/math-ph/0508049"},
}};
inline constexpr std::array<Link, 2> links_69{{
  {"DOI", "https://doi.org/10.1103/PhysRevLett.123.250602"},
  {"arXiv", "https://arxiv.org/abs/1908.08172"},
}};
inline constexpr std::array<Link, 1> links_70{{
  {"Open lecture notes", "https://people.sissa.it/~ffranchi/BAnotes.pdf"},
}};
inline constexpr std::array<Link, 1> links_71{{
  {"Lieb equation", "https://integrability.org/g_l_Le.html"},
}};
inline constexpr std::array<Link, 2> links_72{{
  {"DOI", "https://doi.org/10.1088/1742-5468/2014/10/P10045"},
  {"arXiv", "https://arxiv.org/abs/1407.8344"},
}};
inline constexpr std::array<Link, 2> links_73{{
  {"DOI", "https://doi.org/10.1088/1751-8113/44/10/102001"},
  {"arXiv", "https://arxiv.org/abs/1010.4842"},
}};
inline constexpr std::array<Link, 2> links_74{{
  {"DOI", "https://doi.org/10.1088/0305-4470/39/41/S03"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0611701"},
}};
inline constexpr std::array<Link, 2> links_75{{
  {"DOI", "https://doi.org/10.1088/0305-4470/38/7/001"},
  {"arXiv", "https://arxiv.org/abs/cond-mat/0411505"},
}};

inline constexpr std::array<Reference, 76> references{{
  {"barber-batchelor-1989", "Michael N. Barber and Murray T. Batchelor", "Spectrum of the biquadratic spin-1 antiferromagnetic chain", "Phys. Rev. B 40, 4621-4626", 1989, links_0},
  {"albertini-2000", "Giuseppe Albertini", "Is the purely biquadratic spin 1 chain always massive?", "arXiv:cond-mat/0012439", 2000, links_1},
  {"aufgebauer-klumper-2010", "Britta Aufgebauer and Andreas Klümper", "Quantum spin chains of Temperley-Lieb type: periodic boundary conditions, spectral multiplicities and finite temperature", "J. Stat. Mech. 2010, P05018", 2010, links_2},
  {"essler-korepin-1999", "Fabian H. L. Essler and Vladimir E. Korepin", "Form factors in the half-filled Hubbard model", "Phys. Rev. B 59, 1734-1738", 1999, links_3},
  {"melzer-1995", "Ezer Melzer", "On the scaling limit of the 1D Hubbard model at half filling", "Nucl. Phys. B 443, 553-564", 1995, links_4},
  {"essler-2010", "Fabian H. L. Essler", "Threshold singularities in the one-dimensional Hubbard model", "Phys. Rev. B 81, 205120", 2010, links_5},
  {"luo-pu-guan-2024", "Jia-Jia Luo, Han Pu, and Xi-Wen Guan", "Exact results of the one-dimensional repulsive Hubbard model", "Rep. Prog. Phys. 87, 117601", 2024, links_6},
  {"karbach-1998", "Michael Karbach, Kun Hu, and Gerhard Müller", "Introduction to the Bethe ansatz II", "Computers in Physics 12, 565", 1998, links_7},
  {"groha-2017", "Stefan Groha and Fabian H. L. Essler", "Spinon decay in the spin-1/2 Heisenberg chain with weak next nearest neighbour exchange", "J. Phys. A 50, 334002", 2017, links_8},
  {"caux-xxx-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXX spinons", "Online notes", 0, links_9},
  {"caux-xxz-spinons", "Jean-Sébastien Caux", "The Bethe Ansatz: XXZ spinons", "Online notes", 0, links_10},
  {"mei-2017", "Zhongtao Mei and C. J. Bolech", "Derivation of matrix product states for the Heisenberg spin chain with open boundary conditions", "Phys. Rev. E 95, 032127", 2017, links_11},
  {"caux-xxx-states", "Jean-Sébastien Caux", "The Bethe Ansatz: real rapidities, SU(2) descendants, and two-string states", "Online notes", 0, links_12},
  {"vlijm-2016", "R. Vlijm, I. S. Eliëns, and J.-S. Caux", "Correlations of zero-entropy critical states in the XXZ model: integrability and Luttinger theory far from the ground state", "SciPost Phys. 1, 008", 2016, links_13},
  {"bajnok-2020", "Zoltán Bajnok, Etienne Granet, Jesper Lykke Jacobsen, and Rafael I. Nepomechie", "On Generalized Q-systems", "JHEP 03 (2020) 177", 2020, links_14},
  {"popkov-2021", "Vladislav Popkov, Xin Zhang, and Andreas Klümper", "Phantom Bethe excitations and spin helix eigenstates in integrable periodic and open spin chains", "Phys. Rev. B 104, L081410", 2021, links_15},
  {"kozlowski-2017", "Karol K. Kozlowski", "On condensation properties of Bethe roots associated with the XXZ chain", "arXiv:1508.05741v2", 2017, links_16},
  {"lieb-wu-2003", "Elliott H. Lieb and F. Y. Wu", "The one-dimensional Hubbard model: A reminiscence", "Physica A 321, 1-27", 2003, links_17},
  {"rylands-2022", "Colin Rylands, Bruno Bertini, and Pasquale Calabrese", "Integrable quenches in the Hubbard model", "J. Stat. Mech. 2022, 103103", 2022, links_18},
  {"deguchi-yue-1997", "Tetsuo Deguchi and Ruihong Yue", "Exact solutions of 1-D Hubbard model with open boundary conditions and the conformal dimensions under boundary magnetic fields", "arXiv:cond-mat/9704138", 1997, links_19},
  {"lieb-mattis-1962", "Elliott Lieb and Daniel Mattis", "Theory of Ferromagnetism and the Ordering of Electronic Energy Levels", "Phys. Rev. 125, 164-172", 1962, links_20},
  {"lieb-liniger-1963", "Elliott H. Lieb and Werner Liniger", "Exact Analysis of an Interacting Bose Gas. I. The General Solution and the Ground State", "Phys. Rev. 130, 1605", 1963, links_21},
  {"lieb-1963-excitations", "Elliott H. Lieb", "Exact Analysis of an Interacting Bose Gas. II. The Excitation Spectrum", "Phys. Rev. 130, 1616", 1963, links_22},
  {"gaudin-1971", "M. Gaudin", "Boundary Energy of a Bose Gas in One Dimension", "Phys. Rev. A 4, 386", 1971, links_23},
  {"gaudin-1967", "M. Gaudin", "Un système à une dimension de fermions en interaction", "Phys. Lett. A 24, 55-56", 1967, links_24},
  {"yang-1967", "C. N. Yang", "Some Exact Results for the Many-Body Problem in One Dimension with Repulsive Delta-Function Interaction", "Phys. Rev. Lett. 19, 1312", 1967, links_25},
  {"yang-1968", "C. N. Yang", "S Matrix for the One-Dimensional N-Body Problem with Repulsive or Attractive delta-Function Interaction", "Phys. Rev. 168, 1920", 1968, links_26},
  {"sutherland-1975", "Bill Sutherland", "Model for a multicomponent quantum system", "Phys. Rev. B 12, 3795", 1975, links_27},
  {"essler-korepin-1992", "Fabian H. L. Essler and Vladimir E. Korepin", "Higher conservation laws and algebraic Bethe Ansätze for the supersymmetric t-J model", "Phys. Rev. B 46, 9147", 1992, links_28},
  {"babujian-1982", "H. M. Babujian", "Exact solution of the one-dimensional isotropic Heisenberg chain with arbitrary spins S", "Phys. Lett. A 90, 479-482", 1982, links_29},
  {"vlijm-caux-2014", "Rogier Vlijm and Jean-Sébastien Caux", "Computation of dynamical correlation functions of the spin-1 Babujan-Takhtajan chain", "J. Stat. Mech. 2014, P05009", 2014, links_30},
  {"dukelsky-2004", "J. Dukelsky, S. Pittel, and G. Sierra", "Colloquium: Exactly solvable Richardson-Gaudin models for many-body quantum systems", "Rev. Mod. Phys. 76, 643-662", 2004, links_31},
  {"faribault-2011", "Alexandre Faribault, Omar El Araby, Christoph Sträter, and Vladimir Gritsev", "Gaudin models solver based on the correspondence between Bethe ansatz and ordinary differential equations", "Phys. Rev. B 83, 235124", 2011, links_32},
  {"faribault-schuricht-2013", "Alexandre Faribault and Dirk Schuricht", "Spin decoherence due to a randomly fluctuating spin bath", "Phys. Rev. B 88, 085323", 2013, links_33},
  {"sutherland-1968", "Bill Sutherland", "Further Results for the Many-Body Problem in One Dimension", "Phys. Rev. Lett. 20, 98-100", 1968, links_34},
  {"lee-2011", "J. Y. Lee, X. W. Guan, and M. T. Batchelor", "Yang-Yang method for the thermodynamics of one-dimensional multi-component interacting fermions", "J. Phys. A: Math. Theor. 44, 165002", 2011, links_35},
  {"imambekov-demler-2006", "Adilet Imambekov and Eugene Demler", "Exactly solvable case of a one-dimensional Bose-Fermi mixture", "Phys. Rev. A 73, 021602(R)", 2006, links_36},
  {"imambekov-demler-2006-applications", "Adilet Imambekov and Eugene Demler", "Applications of exact solution for strongly interacting one-dimensional Bose-Fermi mixture: Low-temperature correlation functions, density profiles, and collective modes", "Ann. Phys. 321, 2390", 2006, links_37},
  {"wang-1999", "Yupeng Wang", "Exact solution of a spin-ladder model", "Phys. Rev. B 60, 9236", 1999, links_38},
  {"hakobyan-2004", "Tigran Hakobyan", "The ordering of energy levels for SU(n) symmetric antiferromagnetic chains", "Nucl. Phys. B 699, 575-594", 2004, links_39},
  {"bogoliubov-1997", "N. M. Bogoliubov, A. G. Izergin, and N. A. Kitanine", "Correlation functions for a strongly correlated boson system", "arXiv:solv-int/9710002", 1997, links_40},
  {"zhang-klumper-popkov-2024", "Xin Zhang, Andreas Klümper, and Vladislav Popkov", "Pedestrian's way to Baxter's Bethe ansatz for the periodic XYZ chain", "Phys. Rev. B 109, 115411", 2024, links_41},
  {"baxter-1973", "R. J. Baxter", "Eight-vertex model in lattice statistics and one-dimensional anisotropic Heisenberg chain. III. Eigenvectors of the transfer matrix and Hamiltonian", "Ann. Phys. 76, 48-71", 1973, links_42},
  {"haldane-1988", "F. D. M. Haldane", "Exact Jastrow-Gutzwiller resonating-valence-bond ground state of the spin-1/2 antiferromagnetic Heisenberg chain with 1/r^2 exchange", "Phys. Rev. Lett. 60, 635", 1988, links_43},
  {"shastry-1988", "B. Sriram Shastry", "Exact solution of an S=1/2 Heisenberg antiferromagnetic chain with long-ranged interactions", "Phys. Rev. Lett. 60, 639", 1988, links_44},
  {"jiang-lamers-miao-2026", "Yunfeng Jiang, Jules Lamers, and Yuan Miao", "Norms, overlaps and Yangian descendants for the Haldane-Shastry spin chain", "arXiv:2606.20168v2", 2026, links_45},
  {"sutherland-1971", "Bill Sutherland", "Exact Results for a Quantum Many-Body Problem in One Dimension", "Phys. Rev. A 4, 2019", 1971, links_46},
  {"gurappa-panigrahi-1999", "N. Gurappa and Prasanta K. Panigrahi", "Equivalence of the Sutherland Model to Free Particles on a Circle", "arXiv:hep-th/9908127v2", 1999, links_47},
  {"andrei-1980", "N. Andrei", "Diagonalization of the Kondo Hamiltonian", "Phys. Rev. Lett. 45, 379", 1980, links_48},
  {"hegedus-2026", "Arpad Hegedus", "NLIE formulations for the generalized Gibbs ensemble in the sine-Gordon model", "Nucl. Phys. B 1025, 117385", 2026, links_49},
  {"rutkevich-2020", "S. B. Rutkevich", "On the ground-state energy of the finite sine-Gordon ring", "J. Stat. Mech. 2020, 103101", 2020, links_50},
  {"destri-de-vega-1992", "C. Destri and H. J. de Vega", "New thermodynamic Bethe ansatz equations without strings", "Phys. Rev. Lett. 69, 2313-2317", 1992, links_51},
  {"gwa-spohn-1992", "Leh-Hun Gwa and Herbert Spohn", "Bethe solution for the dynamical-scaling exponent of the noisy Burgers equation", "Phys. Rev. A 46, 844", 1992, links_52},
  {"yang-yang-1966", "C. N. Yang and C. P. Yang", "One-Dimensional Chain of Anisotropic Spin-Spin Interactions. I. Proof of Bethe's Hypothesis for Ground State in a Finite System", "Phys. Rev. 150, 321", 1966, links_53},
  {"shastry-sutherland-1990", "B. Sriram Shastry and Bill Sutherland", "Twisted boundary conditions and effective mass in Heisenberg-Ising and Hubbard rings", "Phys. Rev. Lett. 65, 243", 1990, links_54},
  {"de-vega-gonzalez-ruiz-1994", "H. J. de Vega and A. González-Ruiz", "Boundary K-matrices for the XYZ, XXZ and XXX spin chains", "J. Phys. A: Math. Gen. 27, 6129-6138", 1994, links_55},
  {"essler-de-klerk-2023", "F. H. L. Essler and A. J. J. M. de Klerk", "Statistics of matrix elements of local operators in integrable models", "arXiv:2307.12410v1", 2023, links_56},
  {"yang-yang-1969", "C. N. Yang and C. P. Yang", "Thermodynamics of a One-Dimensional System of Bosons with Repulsive Delta-Function Interaction", "J. Math. Phys. 10, 1115", 1969, links_57},
  {"doikou-nepomechie-1998", "Anastasia Doikou and Rafael I. Nepomechie", "Bulk and Boundary S Matrices for the SU(N) Chain", "Nucl. Phys. B 521, 547-572", 1998, links_58},
  {"oelkers-2006", "N. Oelkers, M. T. Batchelor, M. Bortz, and X.-W. Guan", "Bethe Ansatz study of one-dimensional Bose and Fermi gases with periodic and hard wall boundary conditions", "J. Phys. A 39, 1073-1098", 2006, links_59},
  {"grijalva-2019", "Sebastian Grijalva, Jacopo De Nardis, and Veronique Terras", "Open XXZ chain and boundary modes at zero temperature", "SciPost Phys. 7, 023", 2019, links_60},
  {"dugave-2015", "M. Dugave, F. Göhmann, K. K. Kozlowski, and J. Suzuki", "On form-factor expansions for the XXZ chain in the massive regime", "J. Stat. Mech. (2015) P05037", 2015, links_61},
  {"caux-xxz-coordinate", "Jean-Sébastien Caux", "The Bethe Ansatz: coordinate wavefunctions and periodic XXZ equations", "Online notes", 0, links_62},
  {"korepin-2009", "Vladimir E. Korepin", "Norm of Bethe Wave Function as a Determinant", "arXiv:0911.1881 (historical account of the 1982 norm formula)", 2009, links_63},
  {"bini-robol-2013", "Dario A. Bini and Leonardo Robol", "Solving secular and polynomial equations: a multiprecision algorithm", "Author manuscript, May 10, 2013", 2013, links_64},
  {"cameron-graillat-2022", "Thomas R. Cameron and Stef Graillat", "On a compensated Ehrlich-Aberth method for the accurate computation of all polynomial roots", "Electronic Transactions on Numerical Analysis 55, 401-423", 2022, links_65},
  {"koma-nachtergaele-1997", "Tohru Koma and Bruno Nachtergaele", "The spectral gap of the ferromagnetic XXZ chain", "Lett. Math. Phys. 40, 1-16", 1997, links_66},
  {"zhou-2025-biquadratic", "Huan-Qiang Zhou, Qian-Qian Shi, Ian P. McCulloch, and Murray T. Batchelor", "Goldstone modes and the golden spiral in the ferromagnetic spin-1 biquadratic model", "J. Phys. A: Math. Theor. 58, 39LT01", 2025, links_67},
  {"nachtergaele-spitzer-starr-2007", "Bruno Nachtergaele, Wolfgang Spitzer, and Shannon Starr", "Droplet Excitations for the Spin-1/2 XXZ Chain with Kink Boundary Conditions", "Ann. Henri Poincare 8, 165-201", 2007, links_68},
  {"reichert-2019", "Benjamin Reichert, Grigori E. Astrakharchik, Aleksandra Petkovic, and Zoran Ristivojevic", "Exact Results for the Boundary Energy of One-Dimensional Bosons", "Phys. Rev. Lett. 123, 250602", 2019, links_69},
  {"franchini-2011", "Fabio Franchini", "Notes on Bethe Ansatz Techniques", "SISSA lecture notes, May 15, 2011", 2011, links_70},
  {"caux-lieb-liniger", "Jean-Sébastien Caux", "The Bethe Ansatz: The ground state and the Lieb equation", "Online lecture notes, integrability.org", 0, links_71},
  {"pozsgay-2014-q-boson", "B. Pozsgay", "Quantum quenches and Generalized Gibbs Ensemble in a Bethe Ansatz solvable lattice model of interacting bosons", "J. Stat. Mech. (2014) P10045", 2014, links_72},
  {"guan-batchelor-2011", "X.-W. Guan and M. T. Batchelor", "Polylogs, thermodynamics and scaling functions of one-dimensional quantum many-body systems", "J. Phys. A: Math. Theor. 44, 102001", 2011, links_73},
  {"golinelli-mallick-2006", "O. Golinelli and K. Mallick", "The asymmetric simple exclusion process: an integrable model for non-equilibrium statistical mechanics", "J. Phys. A: Math. Gen. 39, 12679-12705", 2006, links_74},
  {"golinelli-mallick-2005", "O. Golinelli and K. Mallick", "Spectral gap of the totally asymmetric exclusion process at arbitrary filling", "J. Phys. A: Math. Gen. 38, 1419-1425", 2005, links_75},
}};

inline constexpr std::array<Use, 3> uses_sine_gordon_vacuum{{
  {&references[51], "Nonlinear integral equation approach to finite-volume sine-Gordon energies."},
  {&references[50], "Bulk-subtracted finite-ring vacuum scaling function and shifted-contour NLIE; soliton-mass convention."},
  {&references[49], "Integer-coupling D-type TBA used for independent p=2 vacuum validation, equations (2.6)-(2.8)."},
}};
inline constexpr std::array<Use, 3> uses_asep_pbc{{
  {&references[74], "Periodic ASEP Bethe equations (66) and Markov eigenvalues (69), with arbitrary nonnegative hopping rates. Leading relaxation branch followed by continuation; no full-spectrum claim."},
  {&references[75], "Finite-size TASEP leading-relaxation branch used as the continuation seed."},
  {&references[52], "Original periodic asymmetric-exclusion relaxation-gap analysis."},
}};
inline constexpr std::array<Use, 2> uses_tasep_pbc{{
  {&references[75], "Periodic TASEP finite-size Bethe equations (2)-(8) and leading relaxation branch (12)-(14) at arbitrary filling. Complex Markov eigenvalues, not energies; no partial asymmetry or full-spectrum claim."},
  {&references[52], "Original periodic asymmetric-exclusion relaxation-gap analysis. This frontend restricts hopping to the totally asymmetric case."},
}};
inline constexpr std::array<Use, 2> uses_xyz_pbc{{
  {&references[41], "Rectangular XYZ coupling convention (2), regular Bethe equations (46)-(48), and energy normalized to S=sigma/2. Even periodic symmetric ground branch only; no singular-pair or excited-spectrum completeness claim."},
  {&references[42], "Eight-vertex/XYZ Bethe-ansatz foundation; this frontend implements only the even-chain regular ground branch."},
}};
inline constexpr std::array<Use, 3> uses_bose_fermi_pbc{{
  {&references[37], "Equal-mass, equal-repulsion Bose-Fermi Hamiltonian (3) and periodic nested equations (28)-(34). Odd-fermion mixed ground shells; no auxiliary-auxiliary scattering. No correlation functions or trapped-gas approximation."},
  {&references[36], "Original equal-coupling integrable mixture and ground-state study."},
  {&references[21], "Pure-boson ground-state reduction on a ring."},
}};
inline constexpr std::array<Use, 2> uses_q_boson_pbc{{
  {&references[72], "Periodic q-boson Hamiltonian (2.1), Bethe equations (2.13), and shifted energy 4 sum sin^2(k/2). Fixed-N ground states and canonical real-root excitation scans; no quenches or correlation functions."},
  {&references[40], "Deformed Fock algebra (1.2)-(1.8), free/phase limits, and exact phase-model momenta (3.4). Our Hamiltonian is twice (1.1), with the +2N shift included."},
}};
inline constexpr std::array<Use, 3> uses_lieb_liniger_thermal{{
  {&references[57], "Finite-temperature thermodynamics of repulsive continuum bosons, at fixed chemical potential or by inversion at fixed density; hbar=2m=k_B=1."},
  {&references[73], "Yang-Yang pseudoenergy, root density and pressure equations (2)-(5). We solve the full integral equations, not the finite-c polylog expansion."},
  {&references[70], "Section 2.11: thermal filling and entropy of Lieb-Liniger Bethe states. No attractive strings, trapped gases or dynamical correlations."},
}};
inline constexpr std::array<Use, 4> uses_lieb_liniger_dispersion{{
  {&references[21], "Repulsive zero-temperature root-density equation and bulk ground-state energy, in units hbar=2m=1."},
  {&references[22], "Thermodynamic type-I particle and type-II hole excitation branches; energies are fixed-N gaps relative to the ground state. No spectral weights or finite-temperature TBA."},
  {&references[70], "Dressed-energy equation and Fermi boundary condition (2.87)-(2.88), related particle-hole energies (2.95); our phase is +2 atan(k/c)."},
  {&references[71], "Pedagogical ground-state integral equation and normalization at fixed density."},
}};
inline constexpr std::array<Use, 7> uses_biquadratic_obc{{
  {&references[0], "Original free-end spin-1 biquadratic/TL spectral correspondence. Even-chain ground states, TL module minima and restricted real-root excitations of H=-sum (S.S)^2, with representation multiplicities."},
  {&references[1], "Hamiltonians and energy shift (2)-(5), real-root Bethe equations (6)-(10), and integer-label window 1<=I<=N-M. Our spin-half reference is half of (3), after a staggered rotation. The odd-chain spinon band is not implemented."},
  {&references[2], "TL loop weight and quantum-group XXZ end fields, Sec. 2.3; open spin-chain module multiplicities, Sec. 3, especially (48)-(55). Periodic twists, physical-spin decomposition and thermodynamics are not implemented."},
  {&references[66], "Ferromagnetic OBC gap and one-defect band: Proposition 2 and Eq. (3.30), transferred through TL equivalence at Delta=3/2 and multiplied by 2*Delta=3. Gap above the entire ground space is 3-2*cos(pi/N), not a zero-mode splitting."},
  {&references[67], "Ferromagnetic convention H=+sum (S.S)^2, Eq. (3), and Fibonacci ground-space multiplicities under free ends. Ground-space zero modes are distinct from the positive-energy TL defect band; entanglement and ground-state wavefunctions are not implemented."},
  {&references[14], "Open quantum-group-invariant XXZ Q-system, Sec. 5: Wronskian (5.17), Bethe equations (5.12), and admissibility conditions. Used for Q-system module searches (validated through N=8; larger sizes experimental), the regularized AF two-string singlet, and empty-sea two-/three-string ferro branches on odd/even long chains. Original equations also underlie odd/even selected real roots and our restricted high-label scattering windows. Two-strings retain signed real deviations; three-strings retain complex deviations. Neither numerical completeness nor long-chain energy ordering is rigorously certified."},
  {&references[68], "Droplet interpretation and thermodynamic module-edge limit, Theorem 2.1; rescaling by 2*Delta=3 gives thresholds 5/3 for two defects and 2 for three. Finite-chain pair/triple modes solve the original Bethe equations with finite deviations, not a bulk/ideal-string substitution. No form factors or physical-spin decomposition are implemented."},
}};
inline constexpr std::array<Use, 4> uses_xxx_pbc{{
  {&references[7], "Periodic XXX equations, energy normalization, and sector quantum numbers; Eqs. (6)-(9), (16), Table I."},
  {&references[8], "Odd-chain one-spinon states and hole labels, Sec. 3; only the unperturbed integrable chain is used."},
  {&references[9], "Thermodynamic XXX spinon dispersion and the bulk reference for --spinons."},
  {&references[12], "Highest-weight and string-classification background; strings and infinite-root descendants are not implemented."},
}};
inline constexpr std::array<Use, 2> uses_xxx_obc{{
  {&references[11], "Free-end XXX equations: rational limit of Eqs. (11)-(12) and footnote 2; our spin-1/2 normalization differs."},
  {&references[12], "SU(2) and string-classification background, not the boundary equations; scans include only the supported finite-real family."},
}};
inline constexpr std::array<Use, 5> uses_xxz_pbc{{
  {&references[13], "Massless periodic XXZ equations, energy, and momentum, Eqs. (1), (3)-(5); our energy includes the N*Delta/4 shift."},
  {&references[16], "Even-ring negative-anisotropy sector ground states: Eqs. (0.4), (0.7), with Pauli exchange divided by four. We use rank-subtracted equations scaled near Delta=-1; the reference explicitly assumes even length and does not justify odd-ring ground-state selection."},
  {&references[61], "Massive periodic XXZ Hamiltonian and trigonometric Bethe equations, Eqs. (1.1)-(1.2), with zero twist and field; Pauli exchange is divided by four. Used for Delta>1 sector ground states, not form factors, strings or excitations."},
  {&references[7], "XXX limit at Delta=1 and the conventional real-root quantum-number window."},
  {&references[12], "Background for the XXX window used to restrict the XXZ scan; not a complete XXZ state classification."},
}};
inline constexpr std::array<Use, 2> uses_xxz_obc{{
  {&references[11], "Free-end XXZ equations and boundary reflection phase, Eqs. (11)-(12); our Hamiltonian is divided by four and shifted. Negative-Delta ground states use an algebraically rank-subtracted, rescaled form of these equations."},
  {&references[60], "Massive open XXZ boundary root and finite-size deviation, Sec. 4.3.2, specialized to zero boundary fields; used for Delta>1 ground states, not boundary correlations or excitation scans."},
}};
inline constexpr std::array<Use, 4> uses_hubbard_dispersion{{
  {&references[3], "Half-filled spinon and holon dispersions, Eqs. (4)-(5), with our U equal to four times the paper's U; no form factors or spectral weights are implemented."},
  {&references[4], "Exact massive-branch K-Bessel series, resummed here into nonoscillatory positive integrals to retain the weak-coupling charge gap; hopping and interaction conventions are converted as documented."},
  {&references[5], "Doped zero-field density and dressed-energy equations (103)-(106), with dressed momenta from (28)-(29); our mu_unshifted equals mu+2u in (105), u=U/4. Elementary lines only, not threshold exponents or continuum minimization."},
  {&references[6], "Real charge particle/hole and spinon excitation interpretation, Sec. II; doped charge-particle is a real-root addition, not the gapped half-filled antiholon or a k-Lambda string."},
}};
inline constexpr std::array<Use, 2> uses_hubbard_pbc{{
  {&references[17], "Lieb-Wu ground-state equations and quantum-number parity, Eqs. (1), (11), (14)-(18); t=1 and unshifted U*n_up*n_down."},
  {&references[18], "Full particle-hole and partial particle-hole (Shiba) transformations, Sec. II, Eqs. (4)-(7); used for sector mappings, not quench dynamics."},
}};
inline constexpr std::array<Use, 3> uses_hubbard_obc{{
  {&references[19], "Free-end Hubbard energy and nested reflection equations, Eqs. (2.5)-(2.8), (3.1)-(3.2), with all boundary potentials zero; t=1 and unshifted interaction."},
  {&references[20], "Spin ordering for the nearest-neighbour open chain: the sector minimum has S=|Sz|. No periodic-shell spin-branch selection is required."},
  {&references[18], "Bipartite particle-hole and Shiba transformations, Sec. II; used for physical-to-auxiliary sector mappings, not quench dynamics."},
}};
inline constexpr std::array<Use, 3> uses_lieb_liniger_pbc{{
  {&references[21], "Repulsive Bose gas on a ring, H=-sum d_j^2+2c sum delta, and its ground state; c>0."},
  {&references[22], "Excited-state background; we solve finite-volume states in explicit label windows, not a thermodynamic dispersion calculation."},
  {&references[56], "Explicit normalization and finite-ring equations, Eqs. (4), (24)-(27); root-density equation (32)-(33) used for validation, not a thermodynamics API. No matrix elements are implemented."},
}};
inline constexpr std::array<Use, 2> uses_lieb_liniger_obc{{
  {&references[23], "Repulsive Bose gas with Dirichlet walls, positive real roots and reflected scattering. We implement finite-volume states and bounded label windows, not general boundary potentials."},
  {&references[69], "Equations (7)-(8) fix the hard-wall ground-state convention and exclusion of self-image scattering; the finite-volume solver extends the logarithmic labels to selected excitations. No boundary-energy integral-equation or thermodynamics API is implemented."},
}};
inline constexpr std::array<Use, 2> uses_su3_pbc{{
  {&references[27], "Original multicomponent permutation-chain solution; the implementation selects only the fundamental SU(3) periodic balanced ground state."},
  {&references[58], "Nested equations and energy, Eqs. (2.17)-(2.19), logarithmic labels (2.24)-(2.29), and the filled-sea singlet in Sec. 2.3. Our H=sum P gives E=2*E_paper+L. No strings, S matrices or boundary fields are implemented."},
}};
inline constexpr std::array<Use, 2> uses_tj_pbc{{
  {&references[28], "Projected t-J Hamiltonian (1.3)-(1.5), Sutherland BFF equations (3.73), and energy (3.75). We remove the shift 2*N_e-L: E=2*N_h-sum 1/(lambda^2+1/4). Doped mixed-spin coverage is restricted to odd N_up and N_down; J=2t=2."},
  {&references[7], "No-hole reduction to the periodic XXX sector solver: H_tJ=2*H_XXX-L/2. Fermionic translation adds the filled-reference phase (-1)^(L-1)."},
}};
inline constexpr std::array<Use, 2> uses_tb_pbc{{
  {&references[29], "Original integrable higher-spin chain family; only the periodic even-length spin-1 singlet ground state is implemented."},
  {&references[30], "Spin-1 Hamiltonian and complex Bethe equations (1.2)-(1.4), filled two-string sea, and finite-deviation equations (3.8)-(3.11). Our bilinear coefficient is 1, i.e. J=4 in this paper. Equation (3.18) supplies only an initial guess; finite-size deviations are solved, not dropped. No dynamical correlations or excitations are implemented."},
}};
inline constexpr std::array<Use, 2> uses_richardson{{
  {&references[31], "Reduced BCS pairing model and blocked-level sectors. We support attractive uniform pairing with distinct doubly degenerate single-particle levels, not arbitrary degeneracies or the full Gaudin family."},
  {&references[32], "Eigenvalue variables (7), quadratic equations in Sec. II, and Richardson equations (32). We use e_i=2*epsilon_i and y_i=g*Lambda(e_i), continue the filled lowest levels from g=0, and include diagonal pair scattering. Energy is sum e_i*y_i-g*M*(L-M+1), plus blocked single-particle energies. No pair-root reconstruction or form factors are implemented."},
}};
inline constexpr std::array<Use, 2> uses_central_spin{{
  {&references[33], "Spin-half central-spin Hamiltonian and all-down-reference Bethe equations, Eqs. (2)-(6), with no bath field: H=B*S0^z+sum A_j*S0.Sj. Only fixed-magnetization ground energies for distinct nonzero couplings are implemented; no dynamics or form factors."},
  {&references[32], "Quadratic eigenvalue-variable method and central-spin realization, Secs. II and IV.3. We compactify the inverse-field continuation to reach B=0; the seed minimizes the central-spin Hamiltonian, not the Richardson energy."},
}};
inline constexpr std::array<Use, 2> uses_sun_fermions_pbc{{
  {&references[34], "Original multicomponent delta-gas solution. We implement equal-mass repulsive fermions on a periodic ring, not arbitrary mixtures or statistics."},
  {&references[35], "Nested finite-size Bethe equations (2)-(3), with c'=c/2, and logarithmic signs from (13)-(19). Our Hamiltonian is -sum d_j^2+2c sum delta with c>=0 and no Zeeman term. Only the centered all-odd occupied-population ground branch is implemented; no strings, TBA or attraction."},
}};
inline constexpr std::array<Use, 2> uses_ladder_pbc{{
  {&references[38], "Permutation-form ladder Hamiltonian (2), rung basis (3), conserved color populations and chemical-potential form (4), and three nested rational equations (5). We use spin-1/2 operators, leg coefficient 1, four-spin coefficient 4, and J_r=2*J in (2)/(4), restoring E=E_perm-L/4+J_r*(L/4-N_s)-h*(N_+-N_-). Longitudinal field enters through the conserved populations; periodic ground states and singlet-count or fixed-magnetization sector minima only."},
  {&references[39], "SU(n) multiplets, Young diagrams and weight spaces. The periodic energy-ordering theorem has a row-parity restriction; we do not assume each population's own highest weight is lowest. Compatible dominant multiplets and displaced packed real seas are compared, not the full excited-state spectrum."},
}};
inline constexpr std::array<Use, 3> uses_haldane_shastry_pbc{{
  {&references[43], "Original periodic spin-1/2 inverse-square chain; we use H=(pi/N)^2 sum_{i<j} S_i.S_j/sin^2(pi*(i-j)/N), J=1."},
  {&references[44], "Independent exact solution of the long-range antiferromagnet; finite even and odd rings, not nearest-neighbor XXX."},
  {&references[45], "Motif spectral rules and Yangian multiplicities, Secs. 2.2-2.3. Convert their H=sum(1-P)/(4*sin^2) by H_ours=E_ferro-2*(pi/N)^2*H_theirs. No wavefunctions, norms or overlaps are implemented."},
}};
inline constexpr std::array<Use, 2> uses_sutherland_pbc{{
  {&references[46], "Original trigonometric inverse-square gas. We select periodic scalar bosons with collision behavior |x_i-x_j|^lambda, lambda>=0; the coupling coefficient alone does not fix this domain."},
  {&references[47], "Hamiltonian (18), ground energy (19), and the partition spectrum after (23), with beta=lambda and hbar=2m=1. Ascending integer labels include common boosts. Exact energies and momenta only; finite label windows do not claim global spectral completeness."},
}};
inline constexpr std::array<Use, 3> uses_gaudin_yang_pbc{{
  {&references[24], "Original spin-1/2 continuum fermion solution; we implement repulsive periodic ground states in selected sectors, not attraction."},
  {&references[25], "Original multicomponent delta-gas solution; our implementation has two spin components only."},
  {&references[59], "Hamiltonian and periodic fermion equations (1), (2), (7), (8), (25); even N, odd minority population as in Sec. 5. Weak and strong limits (12), (15)-(16) provide checks. Hard walls, attraction and excitations are not implemented."},
}};

enum class Tool { sine_gordon_vacuum, asep_pbc, tasep_pbc, xyz_pbc, bose_fermi_pbc, q_boson_pbc, lieb_liniger_thermal, lieb_liniger_dispersion, biquadratic_obc, xxx_pbc, xxx_obc, xxz_pbc, xxz_obc, hubbard_dispersion, hubbard_pbc, hubbard_obc, lieb_liniger_pbc, lieb_liniger_obc, su3_pbc, tj_pbc, tb_pbc, richardson, central_spin, sun_fermions_pbc, ladder_pbc, haldane_shastry_pbc, sutherland_pbc, gaudin_yang_pbc };

[[nodiscard]] constexpr std::span<Use const> for_tool(Tool tool)
{
  switch (tool)
  {
    case Tool::sine_gordon_vacuum: return uses_sine_gordon_vacuum;
    case Tool::asep_pbc: return uses_asep_pbc;
    case Tool::tasep_pbc: return uses_tasep_pbc;
    case Tool::xyz_pbc: return uses_xyz_pbc;
    case Tool::bose_fermi_pbc: return uses_bose_fermi_pbc;
    case Tool::q_boson_pbc: return uses_q_boson_pbc;
    case Tool::lieb_liniger_thermal: return uses_lieb_liniger_thermal;
    case Tool::lieb_liniger_dispersion: return uses_lieb_liniger_dispersion;
    case Tool::biquadratic_obc: return uses_biquadratic_obc;
    case Tool::xxx_pbc: return uses_xxx_pbc;
    case Tool::xxx_obc: return uses_xxx_obc;
    case Tool::xxz_pbc: return uses_xxz_pbc;
    case Tool::xxz_obc: return uses_xxz_obc;
    case Tool::hubbard_dispersion: return uses_hubbard_dispersion;
    case Tool::hubbard_pbc: return uses_hubbard_pbc;
    case Tool::hubbard_obc: return uses_hubbard_obc;
    case Tool::lieb_liniger_pbc: return uses_lieb_liniger_pbc;
    case Tool::lieb_liniger_obc: return uses_lieb_liniger_obc;
    case Tool::su3_pbc: return uses_su3_pbc;
    case Tool::tj_pbc: return uses_tj_pbc;
    case Tool::tb_pbc: return uses_tb_pbc;
    case Tool::richardson: return uses_richardson;
    case Tool::central_spin: return uses_central_spin;
    case Tool::sun_fermions_pbc: return uses_sun_fermions_pbc;
    case Tool::ladder_pbc: return uses_ladder_pbc;
    case Tool::haldane_shastry_pbc: return uses_haldane_shastry_pbc;
    case Tool::sutherland_pbc: return uses_sutherland_pbc;
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
