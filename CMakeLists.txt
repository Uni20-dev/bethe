cmake_minimum_required(VERSION 3.28)
project(bethe VERSION 0.1.0 LANGUAGES CXX)

option(BETHE_BUILD_APPS "Build Bethe ansatz command-line programs" ${PROJECT_IS_TOP_LEVEL})
option(BETHE_BUILD_TESTS "Build Bethe ansatz regression tests" ${PROJECT_IS_TOP_LEVEL})

# Reuse an enclosing project's Uni20 when present. Otherwise use a reproducible
# source dependency; FETCHCONTENT_SOURCE_DIR_UNI20 provides the local override.
if(NOT TARGET uni20_core)
  # Keep dependency IPO opt-in; enabling it requires coordinating the final
  # application's compiler/linker settings too. Preserve explicit cache choices.
  set(UNI20_ENABLE_LTO OFF CACHE BOOL "Enable LTO/IPO for optimized Uni20 builds")
  if(BETHE_BUILD_APPS)
    set(UNI20_BUILD_CLI ON CACHE BOOL "Build Uni20 command-line application helpers")
  endif()
  include(FetchContent)
  FetchContent_Declare(uni20
    GIT_REPOSITORY https://github.com/Uni20-dev/uni20.git
    GIT_TAG 5aee9d35c39c74f7d4b16e6b976db36836b6acf3
  )
  FetchContent_MakeAvailable(uni20)
  # Keep default builds limited to the components our targets actually link.
  set_property(DIRECTORY "${uni20_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

add_library(bethe INTERFACE)
add_library(bethe::bethe ALIAS bethe)
target_compile_features(bethe INTERFACE cxx_std_23)
target_include_directories(bethe INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(bethe INTERFACE uni20_core uni20_common uni20_linalg)

if(BETHE_BUILD_APPS OR BETHE_BUILD_TESTS)
  # Inspect the actual dependency checkout, including source overrides. Refresh
  # before each build; the script changes the header only when its contents change.
  get_target_property(bethe_uni20_source uni20_common SOURCE_DIR)
  if(NOT EXISTS "${bethe_uni20_source}/run_context.hpp")
    message(FATAL_ERROR
      "Bethe applications and tests require Uni20's run_context API. "
      "Update the supplied Uni20 checkout to include PR56, or use Bethe's pinned dependency.")
  endif()
  set(bethe_generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
  add_custom_target(bethe_build_info
    COMMAND "${CMAKE_COMMAND}" "-DBETHE_SOURCE=${PROJECT_SOURCE_DIR}"
      "-DUNI20_SOURCE=${bethe_uni20_source}" "-DBETHE_VERSION=${PROJECT_VERSION}"
      "-DOUTPUT=${bethe_generated_dir}/bethe-build-info.hpp"
      -P "${PROJECT_SOURCE_DIR}/cmake/WriteBuildInfo.cmake"
    BYPRODUCTS "${bethe_generated_dir}/bethe-build-info.hpp"
    VERBATIM)
  add_library(bethe_cli_output INTERFACE)
  target_link_libraries(bethe_cli_output INTERFACE bethe::bethe)
  target_include_directories(bethe_cli_output INTERFACE "${PROJECT_SOURCE_DIR}/apps" "${bethe_generated_dir}")
  add_dependencies(bethe_cli_output bethe_build_info)
endif()

if(BETHE_BUILD_APPS)
  if(NOT TARGET uni20_cli)
    message(FATAL_ERROR
      "Bethe applications require uni20_cli. Enable UNI20_BUILD_CLI=ON before configuring Uni20, "
      "or use BETHE_BUILD_APPS=OFF for a library-only build.")
  endif()
  add_library(bethe_cli INTERFACE)
  target_link_libraries(bethe_cli INTERFACE bethe_cli_output uni20_cli)
  add_executable(bethe-lieb-liniger-dispersion apps/bethe-lieb-liniger-dispersion.cpp)
  target_link_libraries(bethe-lieb-liniger-dispersion PRIVATE bethe_cli)
  set_target_properties(bethe-lieb-liniger-dispersion PROPERTIES CXX_EXTENSIONS OFF)
  add_executable(bethe-lieb-liniger-thermal apps/bethe-lieb-liniger-thermal.cpp)
  target_link_libraries(bethe-lieb-liniger-thermal PRIVATE bethe_cli)
  set_target_properties(bethe-lieb-liniger-thermal PROPERTIES CXX_EXTENSIONS OFF)
  add_executable(bethe-q-boson-pbc apps/bethe-q-boson-pbc.cpp)
  add_executable(bethe-bose-fermi-pbc apps/bethe-bose-fermi-pbc.cpp)
  target_link_libraries(bethe-bose-fermi-pbc PRIVATE bethe_cli)
  set_target_properties(bethe-bose-fermi-pbc PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(bethe-q-boson-pbc PRIVATE bethe_cli)
  set_target_properties(bethe-q-boson-pbc PROPERTIES CXX_EXTENSIONS OFF)
  foreach(app IN ITEMS bethe-xxx-pbc bethe-xxx-obc bethe-xxz-pbc bethe-xxz-obc bethe-hubbard-pbc bethe-hubbard-obc bethe-hubbard-dispersion bethe-lieb-liniger-pbc bethe-lieb-liniger-obc bethe-su3-pbc bethe-gaudin-yang-pbc bethe-tj-pbc bethe-tb-pbc bethe-biquadratic-obc bethe-richardson bethe-central-spin bethe-sun-fermions-pbc bethe-ladder-pbc bethe-haldane-shastry-pbc bethe-sutherland-pbc)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE bethe_cli)
    set_target_properties(${app} PROPERTIES CXX_EXTENSIONS OFF)
  endforeach()
endif()

if(BETHE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
