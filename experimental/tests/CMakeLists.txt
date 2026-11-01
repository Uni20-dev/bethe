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

add_executable(bethe_experimental_six_spinon_tests test_xxx_six_spinon.cpp)
target_link_libraries(bethe_experimental_six_spinon_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_six_spinon_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_six_spinon_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

add_executable(bethe_experimental_six_spinon_scan_tests test_xxx_six_spinon_scan.cpp)
target_link_libraries(bethe_experimental_six_spinon_scan_tests PRIVATE bethe::experimental TBB::tbb)
target_include_directories(bethe_experimental_six_spinon_scan_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_six_spinon_scan_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)
