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
# GoogleTest discovery forwards PROPERTIES as a list; preserve ENVIRONMENT's
# nested list as one property value for both routine and long-running cases.
string(REPLACE ";" "\\;" bethe_discovered_test_environment "${bethe_test_environment}")

function(bethe_discover_tests target)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "TEST_PREFIX;TIMEOUT" "LONG_RUNNING_TESTS;LABELS")
  if(arg_UNPARSED_ARGUMENTS OR arg_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "Invalid arguments to bethe_discover_tests(${target})")
  endif()
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  target_link_libraries(${target} PRIVATE GTest::gtest_main)
  # Optional family metadata uses the same discovery/environment policy.
  # Existing callers retain their names, labels and timeout defaults.
  set(discovery_options)
  set(routine_properties)
  if(DEFINED arg_TEST_PREFIX)
    list(APPEND discovery_options TEST_PREFIX "${arg_TEST_PREFIX}")
  endif()
  if(arg_LABELS)
    # Preserve the label list through routine_properties and then discovery.
    string(REPLACE ";" "\\\\;" routine_labels "${arg_LABELS}")
    list(APPEND routine_properties LABELS "${routine_labels}")
  endif()
  set(long_timeout 1800)
  if(DEFINED arg_TIMEOUT)
    list(APPEND routine_properties TIMEOUT "${arg_TIMEOUT}")
    set(long_timeout "${arg_TIMEOUT}")
  endif()
  # Discovery happens at test time; the directory's configure-time TESTS
  # property does not include these dynamically registered cases.
  # Partition by GoogleTest name before registration, without duplicating cases
  # or relying on CTest's version-dependent pretty names for typed tests.
  set(routine_filter "*")
  if(arg_LONG_RUNNING_TESTS)
    list(JOIN arg_LONG_RUNNING_TESTS ":" long_filter)
    set(routine_filter "*-${long_filter}")
    set(long_labels "long-running" ${arg_LABELS})
    string(REPLACE ";" "\\;" long_labels "${long_labels}")
    gtest_discover_tests(${target}
      DISCOVERY_MODE PRE_TEST
      NO_PRETTY_TYPES
      TEST_FILTER "${long_filter}"
      TEST_LIST ${target}_LONG_RUNNING_TESTS
      ${discovery_options}
      PROPERTIES
        ENVIRONMENT "${bethe_discovered_test_environment}"
        LABELS "${long_labels}"
        TIMEOUT "${long_timeout}"
    )
  endif()
  gtest_discover_tests(${target}
    DISCOVERY_MODE PRE_TEST
    NO_PRETTY_TYPES
    TEST_FILTER "${routine_filter}"
    ${discovery_options}
    PROPERTIES ENVIRONMENT "${bethe_discovered_test_environment}" ${routine_properties}
  )
endfunction()
