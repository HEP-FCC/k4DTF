// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "Gaudi/Property.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/VertexCollection.h"
#include "k4FWCore/Transformer.h"

#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"

namespace k4DTF {

/// Gaudi functional algorithm wrapping DecayTreeFitter.
///
/// Reads an edm4hep::ReconstructedParticleCollection (the final-state
/// daughters, in the order expected by the DecayDescriptor), runs the decay-
/// tree fit, and writes back:
///   - a fitted ReconstructedParticleCollection (daughters + composite nodes)
///   - a VertexCollection carrying the fitted production vertex
///
/// Steering options (set in the job options / Python steering):
/// @code
///   from Configurables import DTFAlgorithm
///   dtf = DTFAlgorithm("DTF",
///       InputParticles = ["MyDaughters"],
///       OutputParticles = ["FittedDaughters"],
///       OutputVertices  = ["FittedVertices"],
///       DecayDescriptor = "B0 -> K+ pi-",
///       MassConstraints = {"511": 5.27963},
///       VertexConstraint = False,
///       MaxIterations = 100)
/// @endcode
struct DTFAlgorithm final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::ReconstructedParticleCollection,
                     edm4hep::VertexCollection>(
              const edm4hep::ReconstructedParticleCollection&)> {
  DTFAlgorithm(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::ReconstructedParticleCollection,
             edm4hep::VertexCollection>
  operator()(const edm4hep::ReconstructedParticleCollection& inputParticles)
      const override;

private:
  Gaudi::Property<std::string> m_decayDescriptor{
      this, "DecayDescriptor", "",
      "Decay descriptor string, e.g. 'B0 -> K+ pi-'"};

  Gaudi::Property<std::map<std::string, double>> m_massConstraints{
      this, "MassConstraints", {},
      "Map of PDG code (as string) -> mass [GeV] for mass constraints; a mass <= 0 means the nominal mass. "
      "A code also applies to its charge conjugate."};

  Gaudi::Property<double> m_Bz{
      this, "Bz", 2.0, "Magnetic field [T] for the track helices"};
  Gaudi::Property<int> m_maxIterations{
      this, "MaxIterations", 30, "Maximum number of fit iterations"};
  Gaudi::Property<double> m_convergenceEps{
      this, "ConvergenceEps", 1e-6,
      "Convergence criterion on delta chi2 (relative to max(1, ndf, chi2))"};

  std::unique_ptr<k4DTF::DecayDescriptor> m_descriptor;
  k4dtf::DecayTreeFitter m_fitter;
};

}  // namespace k4DTF