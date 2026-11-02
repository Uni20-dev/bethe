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
  add_executable(bethe_example_boundary_form_factors examples/xxx_boundary_form_factors.cpp)
  target_link_libraries(bethe_example_boundary_form_factors PRIVATE bethe::experimental)
endif()

if(BETHE_BUILD_APPS)
  add_executable(bethe-experimental-xxx-boundary-dynamics apps/bethe-xxx-boundary-dynamics.cpp)
  target_link_libraries(bethe-experimental-xxx-boundary-dynamics PRIVATE bethe::experimental bethe_cli TBB::tbb)
  set_target_properties(bethe-experimental-xxx-boundary-dynamics PROPERTIES CXX_EXTENSIONS OFF)
endif()
