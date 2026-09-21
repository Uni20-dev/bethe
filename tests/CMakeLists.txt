add_executable(bethe_heisenberg_tests test_heisenberg.cpp)
target_link_libraries(bethe_heisenberg_tests PRIVATE bethe::bethe)
set_target_properties(bethe_heisenberg_tests PROPERTIES CXX_EXTENSIONS OFF)

foreach(precision IN ITEMS fp64 long-double)
  add_test(NAME Heisenberg.${precision} COMMAND bethe_heisenberg_tests ${precision})
endforeach()
if(UNI20_ENABLE_MPLAPACK)
  add_test(NAME Heisenberg.fp128 COMMAND bethe_heisenberg_tests fp128)
endif()

if(BETHE_BUILD_APPS)
  add_test(NAME Heisenberg.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:heisenberg-energy>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_cli.cmake")
endif()
