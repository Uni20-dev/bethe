# GoogleTest is test-only: neither Bethe's public headers nor consumers acquire
# this dependency when BETHE_BUILD_TESTS=OFF.
if(NOT TARGET GTest::gtest_main)
  find_package(GTest 1.12 QUIET CONFIG)
endif()
if(NOT TARGET GTest::gtest_main)
  include(FetchContent)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
  )
  set(INSTALL_GTEST OFF CACHE BOOL "Install GoogleTest")
  set(BUILD_GMOCK OFF CACHE BOOL "Build GoogleMock")
  if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "Use the shared MSVC runtime for GoogleTest")
  endif()
  FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

function(bethe_discover_tests target)
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(${target} PRIVATE GTest::gtest_main)
  # Discovery happens at test time; the directory's configure-time TESTS
  # property does not include these dynamically registered cases.
  gtest_discover_tests(${target}
    DISCOVERY_MODE PRE_TEST
    NO_PRETTY_TYPES
    PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1\;OMP_NUM_THREADS=1"
  )
endfunction()

foreach(module IN ITEMS heisenberg sectors open excitations xxz xxz_open xxz_open_massive xxz_open_negative xxz_periodic_negative xxz_negative xxz_polynomial xxz_open_qsystem xxz_open_two_string xxz_odd_continuation xxz_odd_bounds xxz_wronskian xxz_odd_sectors xxz_regularity xxz_spin_helix xxz_phantom xxz_phantom_wave xxz_phantom_check xxz_coordinate_wave polynomial_roots hubbard hubbard_open hubbard_thermo hubbard_doped lieb_liniger su3 gaudin_yang tj takhtajan_babujian temperley_lieb biquadratic_ferromagnetic biquadratic_bound_pair biquadratic_bound_triple biquadratic_scattering richardson central_spin su_fermions ladder haldane_shastry sutherland support)
  add_executable(bethe_${module}_tests test_${module}.cpp)
  target_link_libraries(bethe_${module}_tests PRIVATE bethe::bethe)
  bethe_discover_tests(bethe_${module}_tests)
endforeach()

add_executable(bethe_lieb_liniger_open_tests test_lieb_liniger_open.cpp)
target_link_libraries(bethe_lieb_liniger_open_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_lieb_liniger_open_tests)

add_executable(bethe_lieb_liniger_thermo_tests test_lieb_liniger_thermo.cpp)
target_link_libraries(bethe_lieb_liniger_thermo_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_lieb_liniger_thermo_tests)

add_executable(bethe_lieb_liniger_thermal_tests test_lieb_liniger_thermal.cpp)
target_link_libraries(bethe_lieb_liniger_thermal_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_lieb_liniger_thermal_tests)

add_executable(bethe_q_boson_tests test_q_boson.cpp)
target_link_libraries(bethe_q_boson_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_q_boson_tests)

add_executable(bethe_bose_fermi_tests test_bose_fermi.cpp)
target_link_libraries(bethe_bose_fermi_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_bose_fermi_tests)

add_executable(bethe_elliptic_theta_tests test_elliptic_theta.cpp)
target_link_libraries(bethe_elliptic_theta_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_elliptic_theta_tests)

add_executable(bethe_xyz_tests test_xyz.cpp)
target_link_libraries(bethe_xyz_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_xyz_tests)

add_executable(bethe_tasep_tests test_tasep.cpp)
target_link_libraries(bethe_tasep_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_tasep_tests)

add_executable(bethe_asep_tests test_asep.cpp)
target_link_libraries(bethe_asep_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_asep_tests)

add_executable(bethe_sine_gordon_tests test_sine_gordon.cpp)
target_link_libraries(bethe_sine_gordon_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_sine_gordon_tests)

add_executable(bethe_kondo_tests test_kondo.cpp)
target_link_libraries(bethe_kondo_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_kondo_tests)

add_executable(bethe_xxz_quantum_group_critical_tests test_xxz_quantum_group_critical.cpp)
target_link_libraries(bethe_xxz_quantum_group_critical_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_xxz_quantum_group_critical_tests)

add_executable(bethe_hubbard_continuum_tests test_hubbard_continuum.cpp)
target_link_libraries(bethe_hubbard_continuum_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_hubbard_continuum_tests)

add_executable(bethe_extrema_tests test_extrema.cpp)
target_link_libraries(bethe_extrema_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_extrema_tests)

add_executable(bethe_hubbard_charge_continuum_tests test_hubbard_charge_continuum.cpp)
target_link_libraries(bethe_hubbard_charge_continuum_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_hubbard_charge_continuum_tests)

add_executable(bethe_biquadratic_pair_defect_tests test_biquadratic_pair_defect.cpp)
target_link_libraries(bethe_biquadratic_pair_defect_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_biquadratic_pair_defect_tests)
add_executable(bethe_biquadratic_two_pairs_tests test_biquadratic_two_pairs.cpp)
target_link_libraries(bethe_biquadratic_two_pairs_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_biquadratic_two_pairs_tests)
add_executable(bethe_biquadratic_triple_defect_tests test_biquadratic_triple_defect.cpp)
target_link_libraries(bethe_biquadratic_triple_defect_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_biquadratic_triple_defect_tests)
add_executable(bethe_biquadratic_bound_four_tests test_biquadratic_bound_four.cpp)
target_link_libraries(bethe_biquadratic_bound_four_tests PRIVATE bethe::bethe)
bethe_discover_tests(bethe_biquadratic_bound_four_tests)

add_executable(bethe_output_tests test_output.cpp)
target_link_libraries(bethe_output_tests PRIVATE bethe_cli_output)
bethe_discover_tests(bethe_output_tests)

add_executable(bethe_run_metadata_tests test_run_metadata.cpp)
target_link_libraries(bethe_run_metadata_tests PRIVATE bethe_cli_output)
bethe_discover_tests(bethe_run_metadata_tests)

add_executable(bethe_citation_tests test_citations.cpp)
# The citation data API is deliberately standard-library-only.
target_compile_features(bethe_citation_tests PRIVATE cxx_std_23)
target_include_directories(bethe_citation_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
bethe_discover_tests(bethe_citation_tests)

find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME Citations.generated
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/scripts/generate_citations.py" --check)
  add_test(NAME Citations.generator
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_citations.py")
else()
  message(STATUS "Python3 unavailable: skipping citation generation/freshness maintainer tests")
endif()

if(BETHE_BUILD_APPS)
  if(Python3_Interpreter_FOUND)
    add_test(NAME BiquadraticTripleDefect.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_triple_defect.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BiquadraticTwoPairs.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_two_pairs.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BiquadraticPairDefect.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_pair_defect.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BiquadraticBoundPair.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_bound_cluster.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}" 2)
    add_test(NAME BiquadraticBoundTriple.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_bound_cluster.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}" 3)
    add_test(NAME BiquadraticBoundQuartet.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_bound_cluster.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}" 4)
    add_test(NAME BiquadraticFerro.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_ferromagnetic.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BiquadraticSpinContent.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_spin_content.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME LiebLinigerOpen.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_lieb_liniger_open.py"
        "$<TARGET_FILE:bethe-lieb-liniger-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME LiebLinigerThermo.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_lieb_liniger_thermo.py"
        "$<TARGET_FILE:bethe-lieb-liniger-dispersion>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME LiebLinigerThermal.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_lieb_liniger_thermal.py"
        "$<TARGET_FILE:bethe-lieb-liniger-thermal>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME QBoson.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_q_boson.py"
        "$<TARGET_FILE:bethe-q-boson-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BoseFermi.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_bose_fermi.py"
        "$<TARGET_FILE:bethe-bose-fermi-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME XYZ.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_xyz.py"
        "$<TARGET_FILE:bethe-xyz-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME HubbardContinuum.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_continuum.py"
        "$<TARGET_FILE:bethe-hubbard-continuum>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME Arguments.bethe-hubbard-continuum.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-hubbard-continuum>")
    add_test(NAME Citations.bethe-hubbard-continuum.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-hubbard-continuum>"
        "-DTOOL=bethe-hubbard-continuum" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Kondo.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_kondo.py"
        "$<TARGET_FILE:bethe-kondo-response>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME Arguments.bethe-kondo-response.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-kondo-response>")
    add_test(NAME Citations.bethe-kondo-response.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-kondo-response>"
        "-DTOOL=bethe-kondo-response" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME SineGordon.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_sine_gordon.py"
        "$<TARGET_FILE:bethe-sine-gordon-vacuum>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME Arguments.bethe-sine-gordon-vacuum.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-sine-gordon-vacuum>")
    add_test(NAME Citations.bethe-sine-gordon-vacuum.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-sine-gordon-vacuum>"
        "-DTOOL=bethe-sine-gordon-vacuum" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME ASEP.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_asep.py"
        "$<TARGET_FILE:bethe-asep-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME Arguments.bethe-asep-pbc.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-asep-pbc>")
    add_test(NAME Citations.bethe-asep-pbc.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-asep-pbc>"
        "-DTOOL=bethe-asep-pbc" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME TASEP.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_tasep.py"
        "$<TARGET_FILE:bethe-tasep-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME Arguments.bethe-tasep-pbc.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-tasep-pbc>")
    add_test(NAME Citations.bethe-tasep-pbc.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-tasep-pbc>"
        "-DTOOL=bethe-tasep-pbc" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Arguments.bethe-xyz-pbc.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-xyz-pbc>")
    add_test(NAME Citations.bethe-xyz-pbc.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xyz-pbc>"
        "-DTOOL=bethe-xyz-pbc" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Arguments.bethe-bose-fermi-pbc.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-bose-fermi-pbc>")
    add_test(NAME Citations.bethe-bose-fermi-pbc.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-bose-fermi-pbc>"
        "-DTOOL=bethe-bose-fermi-pbc" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Arguments.bethe-q-boson-pbc.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-q-boson-pbc>")
    add_test(NAME Citations.bethe-q-boson-pbc.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-q-boson-pbc>"
        "-DTOOL=bethe-q-boson-pbc" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Arguments.bethe-lieb-liniger-thermal.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-lieb-liniger-thermal>")
    add_test(NAME Citations.bethe-lieb-liniger-thermal.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-lieb-liniger-thermal>"
        "-DTOOL=bethe-lieb-liniger-thermal" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME Arguments.bethe-lieb-liniger-dispersion.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
        "$<TARGET_FILE:bethe-lieb-liniger-dispersion>")
    add_test(NAME Citations.bethe-lieb-liniger-dispersion.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-lieb-liniger-dispersion>"
        "-DTOOL=bethe-lieb-liniger-dispersion" -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
    add_test(NAME BiquadraticTwoString.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_two_string.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME BiquadraticQSystem.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_qsystem.py"
        "$<TARGET_FILE:bethe-biquadratic-obc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME DataOutput.models
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_model_tables.py"
        "${UNI20_ENABLE_MPLAPACK}" "$<TARGET_FILE:bethe-su3-pbc>" "$<TARGET_FILE:bethe-tb-pbc>"
        "$<TARGET_FILE:bethe-richardson>" "$<TARGET_FILE:bethe-central-spin>"
        "$<TARGET_FILE:bethe-gaudin-yang-pbc>" "$<TARGET_FILE:bethe-tj-pbc>"
        "$<TARGET_FILE:bethe-sun-fermions-pbc>" "$<TARGET_FILE:bethe-ladder-pbc>"
        "$<TARGET_FILE:bethe-hubbard-pbc>" "$<TARGET_FILE:bethe-hubbard-obc>"
        "$<TARGET_FILE:bethe-xxx-pbc>" "$<TARGET_FILE:bethe-xxx-obc>"
        "$<TARGET_FILE:bethe-xxz-pbc>" "$<TARGET_FILE:bethe-xxz-obc>"
        "$<TARGET_FILE:bethe-lieb-liniger-pbc>" "$<TARGET_FILE:bethe-biquadratic-obc>")
    add_test(NAME DataOutput.spectra
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_spectrum_tables.py"
        "${UNI20_ENABLE_MPLAPACK}"
        "$<TARGET_FILE:bethe-xxx-pbc>" "$<TARGET_FILE:bethe-xxx-obc>"
        "$<TARGET_FILE:bethe-xxz-pbc>" "$<TARGET_FILE:bethe-xxz-obc>"
        "$<TARGET_FILE:bethe-lieb-liniger-pbc>" "$<TARGET_FILE:bethe-biquadratic-obc>")
    add_test(NAME DataOutput.named
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_named_output.py"
        "$<TARGET_FILE:bethe-sutherland-pbc>")
    add_test(NAME Sutherland.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_sutherland_cli.py"
        "$<TARGET_FILE:bethe-sutherland-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME HaldaneShastry.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_haldane_shastry_cli.py"
        "$<TARGET_FILE:bethe-haldane-shastry-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME HubbardThermo.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_thermo_cli.py"
        "$<TARGET_FILE:bethe-hubbard-dispersion>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME HubbardThermo.arguments
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_arguments.py"
        "$<TARGET_FILE:bethe-hubbard-dispersion>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME DataOutput.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_data_output_cli.py"
        "$<TARGET_FILE:bethe-hubbard-dispersion>" "${UNI20_ENABLE_MPLAPACK}")
  endif()
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc bethe-xxz-pbc bethe-xxz-obc bethe-hubbard-pbc bethe-hubbard-obc bethe-hubbard-dispersion bethe-lieb-liniger-pbc bethe-lieb-liniger-obc bethe-su3-pbc bethe-gaudin-yang-pbc bethe-tj-pbc bethe-tb-pbc bethe-biquadratic-obc bethe-richardson bethe-central-spin bethe-sun-fermions-pbc bethe-ladder-pbc bethe-haldane-shastry-pbc bethe-sutherland-pbc)
    if(Python3_Interpreter_FOUND)
      add_test(NAME Arguments.${app}.cli
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_shared_arguments.py"
          "$<TARGET_FILE:${app}>")
    endif()
    add_test(NAME Citations.${app}.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:${app}>" "-DTOOL=${app}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
  endforeach()
  add_test(NAME Biquadratic.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-biquadratic-obc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_biquadratic_cli.cmake")
  add_test(NAME Ladder.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-ladder-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_ladder_cli.cmake")
  add_test(NAME GaudinYang.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-gaudin-yang-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_gaudin_yang_cli.cmake")
  add_test(NAME TJ.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-tj-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_tj_cli.cmake")
  add_test(NAME TB.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-tb-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_tb_cli.cmake")
  add_test(NAME SUFermions.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-sun-fermions-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_su_fermions_cli.cmake")
  add_test(NAME CentralSpin.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-central-spin>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_central_spin_cli.cmake")
  add_test(NAME Richardson.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-richardson>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_richardson_cli.cmake")
  add_test(NAME SU3.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-su3-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_su3_cli.cmake")
  add_test(NAME LiebLiniger.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-lieb-liniger-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_lieb_liniger_cli.cmake")
  add_test(NAME Hubbard.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-hubbard-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_cli.cmake")
  add_test(NAME HubbardOpen.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-hubbard-obc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_open_cli.cmake")
  add_test(NAME XXZOpen.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xxz-obc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_xxz_open_cli.cmake")
  add_test(NAME XXZ.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xxz-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_xxz_cli.cmake")
  add_test(NAME Heisenberg.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xxx-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_cli.cmake")
  add_test(NAME Heisenberg.presentation
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xxx-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_presentation.cmake")
  add_test(NAME Open.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-xxx-obc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_open_cli.cmake")
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc)
    add_test(NAME Excitations.${app}.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:${app}>"
        "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_excitations_cli.cmake")
  endforeach()
endif()

# Small dense systems and short CLI checks do not benefit from BLAS thread pools.
get_property(bethe_test_names DIRECTORY PROPERTY TESTS)
set_tests_properties(${bethe_test_names} PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
