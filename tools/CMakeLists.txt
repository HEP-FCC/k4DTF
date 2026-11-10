# k4dtf-validate: run k4DTF over a stage-1 ntuple of examples/dtf_test_samples
add_executable(k4dtf-validate k4dtf_validate.cpp)
target_link_libraries(k4dtf-validate
  PRIVATE
    k4DTF::k4DTF
    EDM4HEP::edm4hep
    ROOT::RIO
    ROOT::Tree
    ROOT::TreePlayer
)
install(TARGETS k4dtf-validate RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
