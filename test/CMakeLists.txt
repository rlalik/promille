cmake_minimum_required(VERSION 3.14)

project(mille_builderTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(mille_builder REQUIRED)
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

add_executable(mille_builder_test
  source/test_operations.cpp
  source/test_derivatives.cpp
  source/test_mille.cpp
  source/test_euler_angles.cpp
  source/test_alignment_models.cpp
)

target_link_libraries(mille_builder_test PRIVATE ${GTEST_TRG} mille_builder::mille_builder MP2::Mille)
target_compile_features(mille_builder_test PRIVATE cxx_std_17)

target_compile_options(mille_builder_test PUBLIC -Wmaybe-uninitialized;-Wpedantic)

add_test(NAME mille_builder_test COMMAND mille_builder_test)

# ---- End-of-file commands ----

add_folders(Test)
