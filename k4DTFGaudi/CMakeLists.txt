gaudi_add_module(k4DTFGaudi
  SOURCES
    src/DTFAlgorithm.cpp
  LINK
    k4DTF::k4DTF
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    Gaudi::GaudiKernel
)

target_include_directories(k4DTFGaudi
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS k4DTFGaudi
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
