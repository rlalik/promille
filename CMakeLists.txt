cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    StrawAlignment
    VERSION 0.1.0
    DESCRIPTION "For MillePede with straws"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES NONE
)

include(FetchContent)
FetchContent_Declare(MillePede
  GIT_REPOSITORY /home/dziadu/projects/hades/millepede-ii
  GIT_TAG cmake_integration
)
FetchContent_MakeAvailable(MillePede)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# find ROOT
if (DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()

find_package(ROOT QUIET REQUIRED COMPONENTS GenVector)

# ---- Declare library ----

add_library(StrawAlignment_StrawAlignment INTERFACE)
add_library(StrawAlignment::StrawAlignment ALIAS StrawAlignment_StrawAlignment)

target_link_libraries(StrawAlignment_StrawAlignment INTERFACE ROOT::GenVector)

set_property(
    TARGET StrawAlignment_StrawAlignment PROPERTY
    EXPORT_NAME StrawAlignment
)

target_include_directories(
    StrawAlignment_StrawAlignment ${warning_guard}
    INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(StrawAlignment_StrawAlignment INTERFACE cxx_std_17)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree." "${StrawAlignment_DEVELOPER_MODE}")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

# ---- Developer mode ----

if(NOT StrawAlignment_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of StrawAlignment"
  )
endif()

include(cmake/dev-mode.cmake)
