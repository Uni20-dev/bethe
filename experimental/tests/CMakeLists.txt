add_executable(bethe_experimental_singlet_tests test_xxx_singlet.cpp)
target_link_libraries(bethe_experimental_singlet_tests PRIVATE bethe::experimental GTest::gtest_main)
target_include_directories(bethe_experimental_singlet_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
set_target_properties(bethe_experimental_singlet_tests PROPERTIES CXX_EXTENSIONS OFF)
gtest_discover_tests(bethe_experimental_singlet_tests
  DISCOVERY_MODE PRE_TEST
  NO_PRETTY_TYPES
  TEST_PREFIX "Experimental."
  PROPERTIES ENVIRONMENT "${bethe_discovered_test_environment}" LABELS "experimental" TIMEOUT 60)
