cmake_minimum_required(VERSION 3.14)

project(StrawAlignmentExamples CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

if(PROJECT_IS_TOP_LEVEL)
  find_package(StrawAlignment REQUIRED)
endif()

add_custom_target(run-examples)

function(add_example NAME)
  add_executable("${NAME}" "${NAME}.cpp")
  target_link_libraries("${NAME}" PRIVATE StrawAlignment::StrawAlignment)
  target_compile_features("${NAME}" PRIVATE cxx_std_17)
  add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
  add_dependencies("run_${NAME}" "${NAME}")
  add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(hades_basic_example)

target_sources(hades_basic_example PRIVATE)
target_link_libraries(hades_basic_example PRIVATE MP2::Mille)

add_folders(Example)
