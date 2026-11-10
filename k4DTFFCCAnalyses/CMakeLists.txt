# Generate ROOT dictionary for use in CLING / RDataFrame
root_generate_dictionary(G__k4DTFFCCAnalyses
  ${CMAKE_CURRENT_SOURCE_DIR}/include/k4DTFFCCAnalyses/DTFFunctions.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/k4DTFFCCAnalyses/LinkDef.h
  MODULE k4DTFFCCAnalyses
)

add_library(k4DTFFCCAnalyses SHARED
  src/DTFFunctions.cpp
  src/DTFInput.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/G__k4DTFFCCAnalyses.cxx
)

add_library(k4DTF::k4DTFFCCAnalyses ALIAS k4DTFFCCAnalyses)

target_include_directories(k4DTFFCCAnalyses
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(k4DTFFCCAnalyses
  PUBLIC
    k4DTF::k4DTF
    EDM4HEP::edm4hep
    ROOT::ROOTVecOps
    ROOT::Core
    ROOT::Physics
)

set_target_properties(k4DTFFCCAnalyses PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

install(TARGETS k4DTFFCCAnalyses
  EXPORT k4DTFTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the ROOT pcm/rootmap files generated alongside the dictionary
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/libk4DTFFCCAnalyses.rootmap"
  "${CMAKE_CURRENT_BINARY_DIR}/libk4DTFFCCAnalyses_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
