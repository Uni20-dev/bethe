# Experimental APIs are opt-in build-tree components, never default installs.
add_library(bethe_experimental INTERFACE)
add_library(bethe::experimental ALIAS bethe_experimental)
target_link_libraries(bethe_experimental INTERFACE bethe::bethe)
target_include_directories(bethe_experimental INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")

if(BETHE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
