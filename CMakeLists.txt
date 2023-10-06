cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    promille
    VERSION 0.1.0
    DESCRIPTION "For MillePede with straws"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES NONE
)

include(FetchContent)
FetchContent_Declare(MillePede2
  GIT_REPOSITORY https://gitlab.desy.de/claus.kleinwort/millepede-ii.git
  GIT_TAG main
)
FetchContent_GetProperties(MillePede2)

if(NOT millepede2_POPULATED)
  FetchContent_Populate(MillePede2)

  file(COPY ${PROJECT_SOURCE_DIR}/extras/millepede2/CMakeLists.txt DESTINATION ${millepede2_SOURCE_DIR})

  add_subdirectory(${millepede2_SOURCE_DIR} ${millepede2_BINARY_DIR})
endif()

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(promille_promille INTERFACE)
add_library(promille::promille ALIAS promille_promille)

target_link_libraries(promille_promille INTERFACE MP2::Mille)

set_property(
    TARGET promille_promille PROPERTY
    EXPORT_NAME promille
)

target_include_directories(
    promille_promille ${warning_guard}
    INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(promille_promille INTERFACE cxx_std_17)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES AND CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree." "${promille_DEVELOPER_MODE}")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

# ---- Developer mode ----

if(NOT promille_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of promille"
  )
endif()

include(cmake/dev-mode.cmake)
