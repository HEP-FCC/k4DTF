cmake_minimum_required(VERSION 3.23)

project(k4DTF
  VERSION 0.1.0
  DESCRIPTION "Key4hep Decay Tree Fitter"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile commands for tooling
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Core EG Physics MathCore Matrix)
find_package(EDM4HEP REQUIRED)
find_package(podio REQUIRED)
find_package(Eigen3 REQUIRED)

find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)

# ---------------------------------------------------------------------------
# Sub-packages
# ---------------------------------------------------------------------------
add_subdirectory(k4DTF)
add_subdirectory(k4DTFGaudi)
add_subdirectory(k4DTFFCCAnalyses)

# ---------------------------------------------------------------------------
# Optional tests
# ---------------------------------------------------------------------------
option(BUILD_TESTING "Build unit tests" ON)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
# CMake package config
# ---------------------------------------------------------------------------
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/k4DTFConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  cmake/k4DTFConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/k4DTFConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/k4DTF
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/k4DTFConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/k4DTFConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/k4DTF
)
