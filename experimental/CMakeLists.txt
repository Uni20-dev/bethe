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
  add_executable(bethe_example_finite_six_spinon examples/xxx_six_spinon.cpp)
  target_link_libraries(bethe_example_finite_six_spinon PRIVATE bethe::experimental)
endif()
