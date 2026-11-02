# Experimental APIs are opt-in build-tree components, never default installs.
add_library(bethe_experimental INTERFACE)
add_library(bethe::experimental ALIAS bethe_experimental)
target_link_libraries(bethe_experimental INTERFACE bethe::bethe)
target_include_directories(bethe_experimental INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")

if(BETHE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if(UNIX AND BETHE_BUILD_EXAMPLES)
  add_executable(bethe-experimental-recovery-store-example examples/spinon_recovery_store.cpp)
  target_link_libraries(bethe-experimental-recovery-store-example PRIVATE bethe::experimental)
  set_target_properties(bethe-experimental-recovery-store-example PROPERTIES CXX_EXTENSIONS OFF)
endif()
