add_executable(bethe_experimental_singlet_tests test_xxx_singlet.cpp)
target_link_libraries(bethe_experimental_singlet_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_singlet_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_singlet_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_boundary_factor_tests test_xxx_boundary_form_factors.cpp)
target_link_libraries(bethe_experimental_boundary_factor_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_boundary_factor_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_boundary_factor_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)
