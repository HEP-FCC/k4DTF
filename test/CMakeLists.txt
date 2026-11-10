find_package(Catch2 3 QUIET)
find_package(Python3 COMPONENTS Interpreter QUIET)

if(NOT Catch2_FOUND)
  message(STATUS "Catch2 not found — skipping C++ unit tests (set Catch2_DIR or install Catch2 >= 3)")
else()
  include(Catch)

  add_executable(test_k4DTF test_core.cpp test_toy.cpp)
  target_include_directories(test_k4DTF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(test_k4DTF
    PRIVATE
      k4DTF::k4DTF
      Catch2::Catch2WithMain
  )

  catch_discover_tests(test_k4DTF)
endif()

# ---------------------------------------------------------------------------
# Python / FCCAnalyses pipeline test
# ---------------------------------------------------------------------------
if(NOT Python3_FOUND)
  message(STATUS "Python3 not found — skipping FCCAnalyses pipeline test")
else()
  add_test(NAME test_fccanalyses_pipeline
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_fccanalyses_pipeline.py
  )

  # Prepend the build-tree library directories so the test finds the shared
  # libraries without requiring cmake --install to have been run first.
  set_tests_properties(test_fccanalyses_pipeline PROPERTIES
    ENVIRONMENT
      "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/k4DTFFCCAnalyses:${CMAKE_BINARY_DIR}/k4DTF:$ENV{LD_LIBRARY_PATH}"
  )
endif()
