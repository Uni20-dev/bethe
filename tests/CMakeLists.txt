add_executable(bethe_heisenberg_tests test_heisenberg.cpp)
target_link_libraries(bethe_heisenberg_tests PRIVATE bethe::bethe)
set_target_properties(bethe_heisenberg_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_sector_tests test_sectors.cpp)
target_link_libraries(bethe_sector_tests PRIVATE bethe::bethe)
set_target_properties(bethe_sector_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_open_tests test_open.cpp)
target_link_libraries(bethe_open_tests PRIVATE bethe::bethe)
set_target_properties(bethe_open_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_excitation_tests test_excitations.cpp)
target_link_libraries(bethe_excitation_tests PRIVATE bethe::bethe)
set_target_properties(bethe_excitation_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_xxz_tests test_xxz.cpp)
target_link_libraries(bethe_xxz_tests PRIVATE bethe::bethe)
set_target_properties(bethe_xxz_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_xxz_open_tests test_xxz_open.cpp)
target_link_libraries(bethe_xxz_open_tests PRIVATE bethe::bethe)
set_target_properties(bethe_xxz_open_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_hubbard_tests test_hubbard.cpp)
target_link_libraries(bethe_hubbard_tests PRIVATE bethe::bethe)
set_target_properties(bethe_hubbard_tests PROPERTIES CXX_EXTENSIONS OFF)
add_executable(bethe_citation_tests test_citations.cpp)
# The citation data API is deliberately standard-library-only.
target_compile_features(bethe_citation_tests PRIVATE cxx_std_23)
target_include_directories(bethe_citation_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
set_target_properties(bethe_citation_tests PROPERTIES CXX_EXTENSIONS OFF)
add_test(NAME Citations.registry COMMAND bethe_citation_tests)

find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME Citations.generated
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/scripts/generate_citations.py" --check)
  add_test(NAME Citations.generator
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_citations.py")
else()
  message(STATUS "Python3 unavailable: skipping citation generation/freshness maintainer tests")
endif()

foreach(precision IN ITEMS fp64 long-double)
  add_test(NAME Heisenberg.${precision} COMMAND bethe_heisenberg_tests ${precision})
  add_test(NAME Sectors.${precision} COMMAND bethe_sector_tests ${precision})
  add_test(NAME Open.${precision} COMMAND bethe_open_tests ${precision})
  add_test(NAME Excitations.${precision} COMMAND bethe_excitation_tests ${precision})
  add_test(NAME XXZ.${precision} COMMAND bethe_xxz_tests ${precision})
  add_test(NAME XXZOpen.${precision} COMMAND bethe_xxz_open_tests ${precision})
  add_test(NAME Hubbard.${precision} COMMAND bethe_hubbard_tests ${precision})
endforeach()
if(UNI20_ENABLE_MPLAPACK)
  add_test(NAME Heisenberg.fp128 COMMAND bethe_heisenberg_tests fp128)
  add_test(NAME Sectors.fp128 COMMAND bethe_sector_tests fp128)
  add_test(NAME Open.fp128 COMMAND bethe_open_tests fp128)
  add_test(NAME Excitations.fp128 COMMAND bethe_excitation_tests fp128)
  add_test(NAME XXZ.fp128 COMMAND bethe_xxz_tests fp128)
  add_test(NAME XXZOpen.fp128 COMMAND bethe_xxz_open_tests fp128)
  add_test(NAME Hubbard.fp128 COMMAND bethe_hubbard_tests fp128)
endif()

if(BETHE_BUILD_APPS)
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc bethe-xxz-pbc bethe-xxz-obc bethe-hubbard-pbc)
    add_test(NAME Citations.${app}.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:${app}>" "-DTOOL=${app}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_citations_cli.cmake")
  endforeach()
  add_test(NAME Hubbard.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:bethe-hubbard-pbc>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_hubbard_cli.cmake")
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
