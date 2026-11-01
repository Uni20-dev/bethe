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
endif()

if(BETHE_BUILD_APPS)
  add_executable(bethe-experimental-xxx-four-spinon apps/bethe-xxx-four-spinon.cpp)
  target_link_libraries(bethe-experimental-xxx-four-spinon PRIVATE bethe::experimental bethe_cli TBB::tbb)
  set_target_properties(bethe-experimental-xxx-four-spinon PROPERTIES CXX_EXTENSIONS OFF)
endif()
