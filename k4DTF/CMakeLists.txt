add_library(k4DTF SHARED
  src/DTFParticle.cpp
  src/DecayDescriptor.cpp
  src/DecayTreeFitter.cpp
  src/ParticleTable.cpp
  src/TreeInput.cpp
  src/CalorimeterResolution.cpp
  src/FitResult.cpp
  src/core/FitParams.cpp
  src/core/KalmanCalculator.cpp
  src/core/ParticleNode.cpp
  src/core/InternalParticleNode.cpp
  src/core/ResonanceNode.cpp
  src/core/TrackNode.cpp
  src/core/PhotonNode.cpp
  src/core/CompositeNode.cpp
  src/core/MergedPi0Node.cpp
  src/core/InteractionPointNode.cpp
  src/core/DecayChain.cpp
  src/core/Helix.cpp
)

add_library(k4DTF::k4DTF ALIAS k4DTF)

target_include_directories(k4DTF
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(k4DTF
  PUBLIC
    ROOT::Core
    ROOT::EG
    ROOT::MathCore
    ROOT::Matrix
    ROOT::Physics
    EDM4HEP::edm4hep
    Eigen3::Eigen
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
