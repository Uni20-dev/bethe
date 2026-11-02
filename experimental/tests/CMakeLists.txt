add_executable(bethe_experimental_singlet_tests test_xxx_singlet.cpp)
target_link_libraries(bethe_experimental_singlet_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_singlet_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_singlet_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_complex_root_tests test_xxx_complex_roots.cpp)
target_link_libraries(bethe_experimental_complex_root_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_complex_root_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_complex_root_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_boundary_factor_tests test_xxx_boundary_form_factors.cpp)
target_link_libraries(bethe_experimental_boundary_factor_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_boundary_factor_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_boundary_factor_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_boundary_correlation_tests test_xxx_boundary_correlations.cpp)
target_link_libraries(bethe_experimental_boundary_correlation_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_boundary_correlation_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_boundary_correlation_tests TEST_PREFIX "Experimental."
  LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_boundary_dynamics_tests test_xxx_boundary_dynamics.cpp)
target_link_libraries(bethe_experimental_boundary_dynamics_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_boundary_dynamics_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_boundary_dynamics_tests TEST_PREFIX "Experimental."
  LABELS "experimental" TIMEOUT 60)

if(BETHE_BUILD_APPS)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME Experimental.BoundaryDynamics.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_dynamics_cli.py"
        "$<TARGET_FILE:bethe-experimental-xxx-boundary-dynamics>" "${BETHE_HAS_MPREAL}")
    set_tests_properties(Experimental.BoundaryDynamics.cli PROPERTIES
      LABELS "experimental" TIMEOUT 90 ENVIRONMENT "${bethe_test_environment}")
  endif()
endif()

add_executable(bethe_experimental_boundary_magnetization_tests test_xxx_boundary_magnetization.cpp)
target_link_libraries(bethe_experimental_boundary_magnetization_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_boundary_magnetization_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_boundary_magnetization_tests TEST_PREFIX "Experimental."
  LABELS "experimental" TIMEOUT 90)

if(BETHE_BUILD_APPS AND Python3_Interpreter_FOUND)
  add_test(NAME Experimental.BoundaryMagnetization.cli
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_magnetization_cli.py"
      "$<TARGET_FILE:bethe-experimental-xxx-boundary-magnetization>" "${BETHE_HAS_MPREAL}")
  set_tests_properties(Experimental.BoundaryMagnetization.cli PROPERTIES
    LABELS "experimental" TIMEOUT 90 ENVIRONMENT "${bethe_test_environment}")
endif()

add_executable(bethe_experimental_four_spinon_tests test_xxx_four_spinon.cpp)
target_link_libraries(bethe_experimental_four_spinon_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_four_spinon_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_four_spinon_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

if(UNIX)
  add_executable(bethe_experimental_recovery_store_tests test_spinon_recovery_store.cpp)
  target_link_libraries(bethe_experimental_recovery_store_tests PRIVATE bethe::experimental TBB::tbb)
  target_include_directories(bethe_experimental_recovery_store_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
  bethe_discover_tests(bethe_experimental_recovery_store_tests
    TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

  add_executable(bethe_experimental_observation_store_tests test_spinon_observation_store.cpp)
  target_link_libraries(bethe_experimental_observation_store_tests PRIVATE bethe::experimental TBB::tbb)
  target_include_directories(bethe_experimental_observation_store_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
  bethe_discover_tests(bethe_experimental_observation_store_tests
    TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)
endif()

add_executable(bethe_experimental_observation_group_tests test_spinon_observation_group.cpp)
target_link_libraries(bethe_experimental_observation_group_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_observation_group_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_observation_group_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)
