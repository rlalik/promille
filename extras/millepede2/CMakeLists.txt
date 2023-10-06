cmake_minimum_required(VERSION 3.14)

project(
  MillePede2
  VERSION 4.13.5
  DESCRIPTION
    "Millepede II is a package for linear least squares fits with a large number of parameters. Developed for the alignment and calibration of tracking detectors."
  HOMEPAGE_URL "https://www.desy.de/~kleinwrt/MP2/doc/html/index.html"
  LANGUAGES C CXX Fortran)

option(MP2_INSTALL "Install MP2 files if used as subproject" OFF)
option(MP2_EXTRA_TOOLS "Install extra MP2 tools" OFF)

include(GNUInstallDirs)

include(FortranCInterface)

add_executable(
  pede
  mpdef.f90
  mpmod.f90
  mpdalc.f90
  mptest1.f90
  mptest2.f90
  minresModule.f90
  minresDataModule.f90
  minresqlpDataModule.f90
  minresqlpBlasModule.f90
  minresqlpModule.f90
  mptext.f90
  mpnum.f90
  mphistab.f90
  mpbits.f90
  mpmon.f90
  mpqldec.f90
  Dbandmatrix.f90
  randoms.f90
  mille.f90
  vertpr.f90
  linesrch.f90
  readc.c
  pede.f90)

target_compile_definitions(pede PRIVATE -DREAD_C_FILES)
target_compile_options(pede PRIVATE -Wall -fopenmp -cpp)
target_link_libraries(pede PRIVATE gomp m z gfortran)

add_library(Mille STATIC Mille.cc Mille.h)
add_library(MP2::Mille ALIAS Mille)

target_include_directories(
  Mille PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(Mille PROPERTIES PUBLIC_HEADER "Mille.h")

include(CMakePackageConfigHelpers)

if(NOT DEFINED MillePede2_INSTALL_CMAKEDIR)
  set(MillePede2_INSTALL_CMAKEDIR
      "${CMAKE_INSTALL_LIBDIR}/cmake/MillePede2"
      CACHE STRING "Path to MillePede2 CMake files")
endif()

write_file(
  ${PROJECT_BINARY_DIR}/Config.cmake.in "include(CMakeFindDependencyMacro)
include(\${CMAKE_CURRENT_LIST_DIR}/MillePede2Targets.cmake)")

configure_package_config_file(
  ${PROJECT_BINARY_DIR}/Config.cmake.in
  ${PROJECT_BINARY_DIR}/MillePede2Config.cmake
  INSTALL_DESTINATION ${MillePede2_INSTALL_CMAKEDIR})

write_basic_package_version_file(MillePede2ConfigVersion.cmake
                                 COMPATIBILITY SameMajorVersion)

string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}"
               PROJECT_IS_TOP_LEVEL)

if(PROJECT_IS_TOP_LEVEL OR MP2_INSTALL)
  install(
    TARGETS pede Mille
    EXPORT MillePede2_Targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            COMPONENT Libraries
            NAMELINK_COMPONENT Development
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)

  install(
    EXPORT MillePede2_Targets
    FILE MillePede2Targets.cmake
    NAMESPACE MP2::
    DESTINATION ${MillePede2_INSTALL_CMAKEDIR}
    COMPONENT Development)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MillePede2Config.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/MillePede2ConfigVersion.cmake
          DESTINATION ${MillePede2_INSTALL_CMAKEDIR})

  if(MP2_EXTRA_TOOLS)
    install(PROGRAMS ${PROJECT_SOURCE_DIR}/tools/readMilleBinary.py
                     ${PROJECT_SOURCE_DIR}/tools/compareResults.py TYPE BIN)
    install(FILES ${PROJECT_SOURCE_DIR}/tools/readPedeHists.C
            DESTINATION ${CMAKE_INSTALL_DATADIR}/millepede2)
  endif()
endif()
