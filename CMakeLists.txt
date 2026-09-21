cmake_minimum_required(VERSION 3.24)
project(bethe VERSION 0.1.0 LANGUAGES CXX)

option(BETHE_BUILD_APPS "Build Bethe ansatz command-line programs" ${PROJECT_IS_TOP_LEVEL})
option(BETHE_BUILD_TESTS "Build Bethe ansatz regression tests" ${PROJECT_IS_TOP_LEVEL})

# Reuse an enclosing project's Uni20 when present. Otherwise use a reproducible
# source dependency; FETCHCONTENT_SOURCE_DIR_UNI20 provides the local override.
if(NOT TARGET uni20_core)
  # Uni20's directory-local IPO setting does not reach our final link steps.
  # In particular, Clang LTO archives need a matching LTO-enabled linker.
  # Keep dependency IPO opt-in; CACHE preserves an explicit caller choice.
  set(UNI20_ENABLE_LTO OFF CACHE BOOL "Enable LTO/IPO for optimized Uni20 builds")
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
target_link_libraries(bethe INTERFACE uni20_core uni20_common)

if(BETHE_BUILD_APPS)
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc bethe-xxz-pbc bethe-xxz-obc)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE bethe::bethe)
    set_target_properties(${app} PROPERTIES CXX_EXTENSIONS OFF)
  endforeach()
endif()

if(BETHE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
