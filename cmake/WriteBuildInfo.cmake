# Build provenance, not runtime inspection of whatever directory the user is in.
find_package(Git QUIET)
function(source_revision source result)
  set(revision "unavailable")
  if(Git_FOUND AND IS_DIRECTORY "${source}")
    # An unpacked source archive inside another checkout is not that checkout.
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source}" ls-files --error-unmatch CMakeLists.txt
      RESULT_VARIABLE tracked OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source}" rev-parse HEAD
      RESULT_VARIABLE status OUTPUT_VARIABLE hash OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(tracked EQUAL 0 AND status EQUAL 0 AND hash MATCHES "^[0-9a-f]+$")
      set(revision "${hash}")
      # The target's SOURCE_DIR can be a subdirectory: inspect the whole checkout.
      execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source}" diff --quiet HEAD -- :/
        RESULT_VARIABLE dirty ERROR_QUIET)
      if(dirty EQUAL 1)
        string(APPEND revision "-dirty")
      elseif(NOT dirty EQUAL 0)
        string(APPEND revision "-status-unknown")
      endif()
    endif()
  endif()
  set(${result} "${revision}" PARENT_SCOPE)
endfunction()
source_revision("${BETHE_SOURCE}" bethe_revision)
source_revision("${UNI20_SOURCE}" uni20_revision)
file(CONFIGURE OUTPUT "${OUTPUT}" CONTENT
"// Generated at build time; do not edit.
#pragma once
namespace bethe::cli::build_info {
inline constexpr char version[] = \"@BETHE_VERSION@\";
inline constexpr char revision[] = \"@bethe_revision@\";
inline constexpr char uni20_revision[] = \"@uni20_revision@\";
}
" @ONLY)
