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

add_executable(bethe_experimental_four_spinon_tests test_xxx_four_spinon.cpp)
target_link_libraries(bethe_experimental_four_spinon_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_four_spinon_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_four_spinon_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_four_scan_tests test_xxx_four_spinon_scan.cpp)
target_link_libraries(bethe_experimental_four_scan_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_four_scan_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_four_scan_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

if(BETHE_BUILD_APPS)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME Experimental.FourScan.cli
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_four_scan_cli.py"
        "$<TARGET_FILE:bethe-experimental-xxx-four-spinon>" "${BETHE_HAS_MPREAL}")
    set_tests_properties(Experimental.FourScan.cli PROPERTIES
      LABELS "experimental" TIMEOUT 90 ENVIRONMENT "${bethe_test_environment}")
  endif()
endif()
