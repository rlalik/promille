cmake_minimum_required(VERSION 3.14)

project(promilleTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(promille REQUIRED)
  enable_testing()
endif()

# find ROOT
if (DEFINED ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()

find_package(ROOT QUIET REQUIRED COMPONENTS GenVector)

find_package(GTest)

if(TARGET GTest::gtest_main)
  # new cmake provided
  set(GTEST_TRG GTest::gtest_main)
else()
  # old cmake and pc files
  set(GTEST_TRG gtest gtest_main)
endif()

# ---- Tests ----

add_executable(promille_test
  source/test_euler_angles.cpp
  source/test_parameter_state.cpp
  source/test_measurement_plane.cpp
  source/test_residual_models.cpp
  source/test_promille.cpp
)

target_link_libraries(promille_test PRIVATE ${GTEST_TRG} promille::promille MP2::Mille ROOT::GenVector)
target_compile_features(promille_test PRIVATE cxx_std_17)

target_compile_options(promille_test PUBLIC -Wmaybe-uninitialized;-Wpedantic)

add_test(NAME promille_test COMMAND promille_test)

# ---- End-of-file commands ----

add_folders(Test)
