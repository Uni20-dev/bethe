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

foreach(module IN ITEMS heisenberg sectors open excitations xxz xxz_open xxz_open_massive xxz_open_negative xxz_periodic_negative xxz_negative xxz_polynomial xxz_odd_continuation xxz_odd_bounds xxz_wronskian xxz_odd_sectors xxz_regularity xxz_spin_helix xxz_phantom xxz_phantom_wave xxz_phantom_check xxz_coordinate_wave polynomial_roots hubbard hubbard_open hubbard_thermo hubbard_doped lieb_liniger su3 gaudin_yang tj takhtajan_babujian richardson central_spin su_fermions ladder haldane_shastry support)
  add_executable(bethe_${module}_tests test_${module}.cpp)
  target_link_libraries(bethe_${module}_tests PRIVATE bethe::bethe)
  bethe_discover_tests(bethe_${module}_tests)
endforeach()

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
    add_test(NAME HaldaneShastry.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_haldane_shastry_cli.py"
        "$<TARGET_FILE:bethe-haldane-shastry-pbc>" "${UNI20_ENABLE_MPLAPACK}")
    add_test(NAME HubbardThermo.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_thermo_cli.py"
        "$<TARGET_FILE:bethe-hubbard-dispersion>" "${UNI20_ENABLE_MPLAPACK}")
  endif()
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc bethe-xxz-pbc bethe-xxz-obc bethe-hubbard-pbc bethe-hubbard-obc bethe-hubbard-dispersion bethe-lieb-liniger-pbc bethe-su3-pbc bethe-gaudin-yang-pbc bethe-tj-pbc bethe-tb-pbc bethe-richardson bethe-central-spin bethe-sun-fermions-pbc bethe-ladder-pbc bethe-haldane-shastry-pbc)
    add_test(NAME Citations.${app}.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:${app}>" "-DTOOL=${app}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
  endforeach()
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
