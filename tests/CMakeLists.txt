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

foreach(precision IN ITEMS fp64 long-double)
  add_test(NAME Heisenberg.${precision} COMMAND bethe_heisenberg_tests ${precision})
  add_test(NAME Sectors.${precision} COMMAND bethe_sector_tests ${precision})
  add_test(NAME Open.${precision} COMMAND bethe_open_tests ${precision})
  add_test(NAME Excitations.${precision} COMMAND bethe_excitation_tests ${precision})
endforeach()
if(UNI20_ENABLE_MPLAPACK)
  add_test(NAME Heisenberg.fp128 COMMAND bethe_heisenberg_tests fp128)
  add_test(NAME Sectors.fp128 COMMAND bethe_sector_tests fp128)
  add_test(NAME Open.fp128 COMMAND bethe_open_tests fp128)
  add_test(NAME Excitations.fp128 COMMAND bethe_excitation_tests fp128)
endif()

if(BETHE_BUILD_APPS)
  add_test(NAME Heisenberg.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:heisenberg-energy>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_cli.cmake")
  add_test(NAME Heisenberg.presentation
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:heisenberg-energy>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_presentation.cmake")
  add_test(NAME Open.cli
    COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:heisenberg-open-energy>"
      "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/check_open_cli.cmake")
  foreach(app IN ITEMS heisenberg-energy heisenberg-open-energy)
    add_test(NAME Excitations.${app}.cli
      COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:${app}>"
        "-DFLOAT128=${UNI20_ENABLE_MPLAPACK}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/check_excitations_cli.cmake")
  endforeach()
endif()
