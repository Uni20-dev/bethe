# GoogleTest is test-only: neither Bethe's public headers nor consumers acquire
# this dependency when BETHE_BUILD_TESTS=OFF.
# Capabilities are derived in Uni20's directory scope, not cache requests.
# In particular AUTO and scalar-only builds must not imply provider support.
foreach(capability IN ITEMS FLOAT128 MPREAL MPLAPACK_FLOAT128 MPLAPACK_MPREAL)
  get_directory_property(BETHE_HAS_${capability}
    DIRECTORY "${bethe_uni20_core_source}" DEFINITION UNI20_HAS_${capability})
endforeach()
if(NOT TARGET GTest::gtest_main)
  find_package(GTest 1.12 QUIET CONFIG)
endif()
if(NOT TARGET GTest::gtest_main)
  include(FetchContent)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
  )
  set(INSTALL_GTEST OFF CACHE BOOL "Install GoogleTest")
  set(BUILD_GMOCK OFF CACHE BOOL "Build GoogleMock")
  if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "Use the shared MSVC runtime for GoogleTest")
  endif()
  FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

# Dependency lifetime policy belongs to Uni20. Resolve its actual source tree,
# including enclosing projects and FETCHCONTENT_SOURCE_DIR_UNI20 overrides.
# This opt-in affects test processes only, not applications or library targets.
option(BETHE_TEST_UNI20_LSAN_SUPPRESSIONS "Use Uni20's LeakSanitizer policy in Bethe tests" OFF)
set(bethe_test_environment "OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
if(BETHE_TEST_UNI20_LSAN_SUPPRESSIONS)
  get_filename_component(bethe_uni20_lsan_suppressions
    "${bethe_uni20_core_source}/../../../tests/sanitizers/lsan.supp" ABSOLUTE)
  if(NOT EXISTS "${bethe_uni20_lsan_suppressions}")
    message(FATAL_ERROR
      "Uni20's LeakSanitizer policy is missing: ${bethe_uni20_lsan_suppressions}. "
      "Update the supplied Uni20 checkout, or use Bethe's pinned dependency.")
  endif()
  list(APPEND bethe_test_environment
    "LSAN_OPTIONS=suppressions=${bethe_uni20_lsan_suppressions}:print_suppressions=0")
endif()
# Quote generated CTest script values without any intermediate list expansion.
function(bethe_test_literal output value)
  set(equals "=")
  while("${value}" MATCHES "\\]${equals}\\]")
    string(APPEND equals "=")
  endwhile()
  set(${output} "[${equals}[${value}]${equals}]" PARENT_SCOPE)
endfunction()

function(bethe_discover_tests target)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "TEST_PREFIX;TIMEOUT" "LONG_RUNNING_TESTS;LABELS")
  if(arg_UNPARSED_ARGUMENTS OR arg_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "Invalid arguments to bethe_discover_tests(${target})")
  endif()
  if(arg_TEST_PREFIX MATCHES "[][;]")
    message(FATAL_ERROR "Test prefixes cannot contain brackets or list separators")
  endif()
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(${target} PRIVATE GTest::gtest_main)
  set(discovery_options)
  if(DEFINED arg_TEST_PREFIX)
    list(APPEND discovery_options TEST_PREFIX "${arg_TEST_PREFIX}")
  endif()
  set(long_timeout 1800)
  if(DEFINED arg_TIMEOUT)
    set(long_timeout "${arg_TIMEOUT}")
  endif()
  # Keep original GoogleTest identifiers: pretty values can contain brackets,
  # which CMake omits from TEST_LIST. All discovered names must receive policy.
  # Partition by the original GoogleTest name, without registering duplicates.
  set(routine_filter "*")
  if(arg_LONG_RUNNING_TESTS)
    list(JOIN arg_LONG_RUNNING_TESTS ":" long_filter)
    set(routine_filter "*-${long_filter}")
    gtest_discover_tests(${target}
      DISCOVERY_MODE PRE_TEST
      NO_PRETTY_TYPES NO_PRETTY_VALUES
      TEST_FILTER "${long_filter}"
      TEST_LIST ${target}_LONG_RUNNING_TESTS
      ${discovery_options}
    )
  endif()
  gtest_discover_tests(${target}
    DISCOVERY_MODE PRE_TEST
    NO_PRETTY_TYPES NO_PRETTY_VALUES
    TEST_FILTER "${routine_filter}"
    ${discovery_options}
  )
  # PROPERTIES is repeatedly flattened by GoogleTest discovery, differently
  # across supported CMake versions. Apply list-valued properties directly in
  # CTest, after discovery, using its documented TEST_INCLUDE_FILES mechanism.
  bethe_test_literal(environment_literal "${bethe_test_environment}")
  bethe_test_literal(routine_labels_literal "${arg_LABELS}")
  set(long_labels long-running ${arg_LABELS})
  bethe_test_literal(long_labels_literal "${long_labels}")
  bethe_test_literal(routine_timeout_literal "${arg_TIMEOUT}")
  bethe_test_literal(long_timeout_literal "${long_timeout}")
  set(properties_file "${CMAKE_CURRENT_BINARY_DIR}/${target}_properties.cmake")
  configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/BetheTestProperties.cmake.in"
    "${properties_file}" @ONLY)
  set_property(DIRECTORY APPEND PROPERTY TEST_INCLUDE_FILES "${properties_file}")
endfunction()
