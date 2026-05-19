add_library(k4DTF SHARED
  src/DTFParticle.cpp
  src/DecayDescriptor.cpp
  src/DecayTreeFitter.cpp
)

add_library(k4DTF::k4DTF ALIAS k4DTF)

target_include_directories(k4DTF
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(k4DTF
  PUBLIC
    ROOT::Core
    ROOT::EG
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
    EDM4HEP::edm4hep
)

set_target_properties(k4DTF PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

install(TARGETS k4DTF
  EXPORT k4DTFTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT k4DTFTargets
  NAMESPACE k4DTF::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/k4DTF
)
