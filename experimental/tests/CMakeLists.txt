add_executable(bethe_experimental_singlet_tests test_xxx_singlet.cpp)
target_link_libraries(bethe_experimental_singlet_tests PRIVATE bethe::experimental)
target_include_directories(bethe_experimental_singlet_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
bethe_discover_tests(bethe_experimental_singlet_tests
  TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)

if(UNIX)
  add_executable(bethe_experimental_recovery_store_tests test_spinon_recovery_store.cpp)
  target_link_libraries(bethe_experimental_recovery_store_tests PRIVATE bethe::experimental TBB::tbb)
  target_include_directories(bethe_experimental_recovery_store_tests PRIVATE "${PROJECT_SOURCE_DIR}/tests")
  bethe_discover_tests(bethe_experimental_recovery_store_tests
    TEST_PREFIX "Experimental." LABELS "experimental" TIMEOUT 60)
endif()
