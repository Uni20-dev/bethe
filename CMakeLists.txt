cmake_minimum_required(VERSION 3.24)
project(bethe VERSION 0.1.0 LANGUAGES CXX)

option(BETHE_BUILD_APPS "Build Bethe ansatz command-line programs" ${PROJECT_IS_TOP_LEVEL})
option(BETHE_BUILD_TESTS "Build Bethe ansatz regression tests" ${PROJECT_IS_TOP_LEVEL})

# Reuse an enclosing project's Uni20 when present. Otherwise use a reproducible
# source dependency; FETCHCONTENT_SOURCE_DIR_UNI20 provides the local override.
if(NOT TARGET uni20_core)
  include(FetchContent)
  FetchContent_Declare(uni20
    GIT_REPOSITORY https://github.com/Uni20-dev/uni20.git
    GIT_TAG a13cc911b58a670da48fba637ed1d13fd2d4f717
  )
  FetchContent_MakeAvailable(uni20)
  # Keep default builds limited to the components our targets actually link.
  set_property(DIRECTORY "${uni20_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

add_library(bethe INTERFACE)
add_library(bethe::bethe ALIAS bethe)
target_compile_features(bethe INTERFACE cxx_std_23)
target_include_directories(bethe INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(bethe INTERFACE uni20_core)

if(BETHE_BUILD_APPS)
  add_executable(heisenberg-energy apps/heisenberg-energy.cpp)
  target_link_libraries(heisenberg-energy PRIVATE bethe::bethe)
  set_target_properties(heisenberg-energy PROPERTIES CXX_EXTENSIONS OFF)
endif()

if(BETHE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
