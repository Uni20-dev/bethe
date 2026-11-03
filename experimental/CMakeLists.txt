# Experimental APIs are opt-in build-tree components, never default installs.
add_library(bethe_experimental INTERFACE)
add_library(bethe::experimental ALIAS bethe_experimental)
target_link_libraries(bethe_experimental INTERFACE bethe::bethe)
target_include_directories(bethe_experimental INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")

if(BETHE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(BETHE_BUILD_EXAMPLES)
  add_executable(bethe-experimental-four-spinon-example examples/xxx_four_spinon.cpp)
  target_link_libraries(bethe-experimental-four-spinon-example PRIVATE bethe::experimental)
  set_target_properties(bethe-experimental-four-spinon-example PROPERTIES CXX_EXTENSIONS OFF)

  add_executable(bethe_example_xxz_vertex examples/xxz_vertex.cpp)
  target_link_libraries(bethe_example_xxz_vertex PRIVATE bethe::experimental TBB::tbb)
  if(BETHE_BUILD_TESTS)
    add_test(NAME Experimental.XXZVertex.example COMMAND bethe_example_xxz_vertex)
    set_tests_properties(Experimental.XXZVertex.example PROPERTIES LABELS "experimental" TIMEOUT 60
      ENVIRONMENT "${bethe_test_environment}")
  endif()

  add_executable(bethe_example_boundary_form_factors examples/xxx_boundary_form_factors.cpp)
  target_link_libraries(bethe_example_boundary_form_factors PRIVATE bethe::experimental)
endif()

if(BETHE_BUILD_APPS)
  add_executable(bethe-experimental-xxx-boundary-dynamics apps/bethe-xxx-boundary-dynamics.cpp)
  target_link_libraries(bethe-experimental-xxx-boundary-dynamics PRIVATE bethe::experimental bethe_cli TBB::tbb)
  set_target_properties(bethe-experimental-xxx-boundary-dynamics PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(BETHE_BUILD_APPS)
  add_executable(bethe-experimental-xxx-boundary-magnetization apps/bethe-xxx-boundary-magnetization.cpp)
  target_link_libraries(bethe-experimental-xxx-boundary-magnetization PRIVATE bethe::experimental bethe_cli TBB::tbb)
  set_target_properties(bethe-experimental-xxx-boundary-magnetization PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(BETHE_BUILD_EXAMPLES)
  add_executable(bethe_example_finite_six_spinon examples/xxx_six_spinon.cpp)
  target_link_libraries(bethe_example_finite_six_spinon PRIVATE bethe::experimental)
  add_executable(bethe_example_six_spinon_scan examples/xxx_six_spinon_scan.cpp)
  target_link_libraries(bethe_example_six_spinon_scan PRIVATE bethe::experimental TBB::tbb)
  set_target_properties(bethe_example_six_spinon_scan PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(BETHE_BUILD_APPS)
  add_executable(bethe-experimental-xxx-four-spinon apps/bethe-xxx-four-spinon.cpp)
  target_link_libraries(bethe-experimental-xxx-four-spinon PRIVATE bethe::experimental bethe_cli TBB::tbb)
  set_target_properties(bethe-experimental-xxx-four-spinon PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(UNIX AND BETHE_BUILD_EXAMPLES)
  add_executable(bethe-experimental-recovery-store-example examples/spinon_recovery_store.cpp)
  target_link_libraries(bethe-experimental-recovery-store-example PRIVATE bethe::experimental)
  set_target_properties(bethe-experimental-recovery-store-example PROPERTIES CXX_EXTENSIONS OFF)
  add_executable(bethe-experimental-observation-store-example examples/spinon_observation_store.cpp)
  target_link_libraries(bethe-experimental-observation-store-example PRIVATE bethe::experimental TBB::tbb)
  set_target_properties(bethe-experimental-observation-store-example PROPERTIES CXX_EXTENSIONS OFF)
endif()
