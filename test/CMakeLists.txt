cmake_minimum_required(VERSION 3.14)

project(StrawAlignmentTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(StrawAlignment REQUIRED)
  enable_testing()
endif()

find_package(GTest)

if(TARGET GTest::gtest_main)
  # new cmake provided
  set(GTEST_TRG GTest::gtest_main)
else()
  # old cmake and pc files
  set(GTEST_TRG gtest gtest_main)
endif()

# ---- Tests ----

add_executable(StrawAlignment_test
  source/test_rotation_matrix.cpp
  source/test_operations.cpp
  source/test_derivatives.cpp
  source/test_mille.cpp
  source/test_euler_angles.cpp
)

target_link_libraries(StrawAlignment_test PRIVATE ${GTEST_TRG} StrawAlignment::StrawAlignment MP2::Mille)
target_compile_features(StrawAlignment_test PRIVATE cxx_std_17)

target_compile_options(StrawAlignment_test PUBLIC -Wmaybe-uninitialized;-Wpedantic)

add_test(NAME StrawAlignment_test COMMAND StrawAlignment_test)

# ---- End-of-file commands ----

add_folders(Test)
