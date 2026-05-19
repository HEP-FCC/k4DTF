// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTFGaudi/DTFAlgorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/VertexCollection.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/FitResult.h"

namespace k4DTF {

// ---------------------------------------------------------------------------

DTFAlgorithm::DTFAlgorithm(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
          name, svcLoc,
          {KeyValue("InputParticles", "InputParticles")},
          {KeyValue("OutputParticles", "FittedParticles"),
           KeyValue("OutputVertices", "FittedVertices")}) {}

// ---------------------------------------------------------------------------

StatusCode DTFAlgorithm::initialize() {
  const StatusCode sc = MultiTransformer::initialize();
  if (sc.isFailure()) return sc;

  if (m_decayDescriptor.value().empty()) {
    error() << "DecayDescriptor property must be set" << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    auto desc = k4DTF::DecayDescriptor::fromString(m_decayDescriptor);
    m_fitter = std::make_unique<k4DTF::DecayTreeFitter>(std::move(desc));
  } catch (const std::exception& e) {
    error() << "Failed to parse DecayDescriptor: " << e.what() << endmsg;
    return StatusCode::FAILURE;
  }

  m_fitter->setMaxIterations(m_maxIterations)
      .setConvergenceCriterion(m_convergenceEps)
      .setVertexConstraint(m_vertexConstraint);

  for (const auto& [pdgStr, mass] : m_massConstraints.value()) {
    try {
      m_fitter->addMassConstraint(std::stoi(pdgStr), mass);
    } catch (...) {
      warning() << "Could not parse PDG code '" << pdgStr
                << "' in MassConstraints — skipping" << endmsg;
    }
  }

  info() << "Decay descriptor: " << m_decayDescriptor.value() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode DTFAlgorithm::finalize() {
  m_fitter.reset();
  return MultiTransformer::finalize();
}

// ---------------------------------------------------------------------------

std::tuple<edm4hep::ReconstructedParticleCollection,
           edm4hep::VertexCollection>
DTFAlgorithm::operator()(
    const edm4hep::ReconstructedParticleCollection& inputParticles) const {
  edm4hep::ReconstructedParticleCollection fittedParticles;
  edm4hep::VertexCollection fittedVertices;

  // --- convert edm4hep → DTFParticle ---
  // ReconstructedParticle is a proxy; extract ReconstructedParticleData first.
  std::vector<k4DTF::DTFParticle> dtfParticles;
  dtfParticles.reserve(inputParticles.size());
  for (const auto& rp : inputParticles) {
    edm4hep::ReconstructedParticleData data{};
    data.momentum  = rp.getMomentum();
    data.energy    = rp.getEnergy();
    data.charge    = rp.getCharge();
    data.PDG       = rp.getPDG();
    data.covMatrix = rp.getCovMatrix();
    dtfParticles.push_back(k4DTF::DTFParticle::fromEDM4hep(data));
  }

  // --- run the fit ---
  const k4DTF::FitResult result = m_fitter->fit(dtfParticles);

  if (!result.converged()) {
    debug() << "DTF did not converge after " << result.nIter
            << " iterations (chi2=" << result.chi2 << ")" << endmsg;
    // Return empty collections — caller decides how to handle
    return {std::move(fittedParticles), std::move(fittedVertices)};
  }

  // --- convert DTFParticle → edm4hep::ReconstructedParticle ---
  for (const auto& fp : result.fittedParticles)
    fittedParticles.create(fp.toEDM4hep());

  // --- fill vertex ---
  // TODO: store fitted vertex position and covariance once vertex constraints
  //       are implemented in the core fitter.
  auto vtx = fittedVertices.create();
  vtx.setChi2(static_cast<float>(result.chi2));
  vtx.setNdf(result.ndf);

  return {std::move(fittedParticles), std::move(fittedVertices)};
}

// ---------------------------------------------------------------------------

DECLARE_COMPONENT(DTFAlgorithm)

}  // namespace k4DTF