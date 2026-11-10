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
#include "k4DTF/TreeInput.h"
#include "edm4hep/TrackCollection.h"

#include <map>

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
    m_descriptor = std::make_unique<k4DTF::DecayDescriptor>(k4DTF::DecayDescriptor::fromString(m_decayDescriptor));
  } catch (const std::exception& e) {
    error() << "Failed to parse DecayDescriptor: " << e.what() << endmsg;
    return StatusCode::FAILURE;
  }
  for (const auto& [pdgStr, mass] : m_massConstraints.value()) {
    try {
      (void)std::stoi(pdgStr);
    } catch (...) {
      error() << "Could not parse PDG code '" << pdgStr << "' in MassConstraints" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  k4dtf::DecayTreeFitter::Configuration cfg;
  cfg.maxIter = m_maxIterations;
  cfg.convergenceCriterion = m_convergenceEps;
  m_fitter.setConfiguration(cfg);

  info() << "Decay descriptor: " << m_descriptor->toString() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode DTFAlgorithm::finalize() {
  m_descriptor.reset();
  return MultiTransformer::finalize();
}

namespace {

/// Set explicit constraint masses (value > 0) on the nodes with these PDG codes.
void setConstraintMasses(k4dtf::DecayNodeInput& in, const std::map<int, double>& masses) {
  for (const auto& [pdg, m] : masses)
    if (m > 0 && (pdg == in.pdgCode || pdg == -in.pdgCode)) in.massConstraintValue = m;
  for (auto& d : in.daughters) setConstraintMasses(d, masses);
}

}  // namespace

// ---------------------------------------------------------------------------
// The input collection holds the final-state particles of one candidate, in
// the depth-first leaf order of the decay descriptor. Charged particles are
// fitted with the first track state of their first track.
std::tuple<edm4hep::ReconstructedParticleCollection,
           edm4hep::VertexCollection>
DTFAlgorithm::operator()(
    const edm4hep::ReconstructedParticleCollection& inputParticles) const {
  edm4hep::ReconstructedParticleCollection fittedParticles;
  edm4hep::VertexCollection fittedVertices;

  // --- convert edm4hep -> DTFParticle ---
  std::vector<k4DTF::DTFParticle> leaves;
  leaves.reserve(inputParticles.size());
  for (const auto& rp : inputParticles) {
    edm4hep::ReconstructedParticleData data{};
    data.momentum  = rp.getMomentum();
    data.energy    = rp.getEnergy();
    data.mass      = rp.getMass();
    data.charge    = rp.getCharge();
    data.PDG       = rp.getPDG();
    data.covMatrix = rp.getCovMatrix();
    const auto tracks = rp.getTracks();
    if (!tracks.empty() && !tracks[0].getTrackStates().empty())
      leaves.push_back(k4DTF::DTFParticle::fromEDM4hep(data, tracks[0].getTrackStates()[0], m_Bz));
    else
      leaves.push_back(k4DTF::DTFParticle::fromEDM4hep(data));
  }

  std::vector<int> constrained;
  std::map<int, double> masses;
  for (const auto& [pdgStr, m] : m_massConstraints.value()) {
    constrained.push_back(std::stoi(pdgStr));
    masses[std::stoi(pdgStr)] = m;
  }

  // --- run the fit ---
  k4dtf::FitResult result;
  try {
    auto input = k4DTF::makeDecayNodeInput(*m_descriptor, leaves, constrained);
    setConstraintMasses(input, masses);
    result = m_fitter.fit(input);
  } catch (const std::exception& e) {
    debug() << "DTF input error: " << e.what() << endmsg;
    return {std::move(fittedParticles), std::move(fittedVertices)};
  }
  if (!result.converged()) {
    debug() << "DTF failed: " << result.status << " after " << result.nIter << " iterations" << endmsg;
    return {std::move(fittedParticles), std::move(fittedVertices)};
  }

  // --- every node of the tree (head first), and a vertex per decay vertex ---
  std::vector<edm4hep::MutableReconstructedParticle> parts;
  for (const auto& n : result.nodes) {
    auto p = fittedParticles.create();
    p.setPDG(n.pdgCode);
    p.setMomentum({static_cast<float>(n.p4(0)), static_cast<float>(n.p4(1)), static_cast<float>(n.p4(2))});
    p.setEnergy(static_cast<float>(n.p4(3)));
    p.setMass(static_cast<float>(n.mass));
    edm4hep::CovMatrix4f cov;
    for (int i = 0, k = 0; i < 4; ++i)
      for (int j = 0; j <= i; ++j, ++k) cov[k] = static_cast<float>(n.p4Cov(i, j));
    p.setCovMatrix(cov);
    parts.push_back(p);
  }
  for (std::size_t i = 0; i < result.nodes.size(); ++i) {
    const auto& n = result.nodes[i];
    for (int c : n.children) parts[i].addToParticles(parts[c]);
    if (!n.hasDecayVertex || n.children.empty()) continue;
    auto vtx = fittedVertices.create();
    vtx.setPosition({static_cast<float>(n.decayVertex(0)), static_cast<float>(n.decayVertex(1)),
                     static_cast<float>(n.decayVertex(2))});
    edm4hep::CovMatrix3f vcov;
    for (int a = 0, k = 0; a < 3; ++a)
      for (int b = 0; b <= a; ++b, ++k) vcov[k] = static_cast<float>(n.decayVertexCov(a, b));
    vtx.setCovMatrix(vcov);
    vtx.setChi2(static_cast<float>(result.chi2));  // chi2 and ndf of the whole tree
    vtx.setNdf(result.ndf);
    for (int c : n.children) vtx.addToParticles(parts[c]);
    parts[i].setDecayVertex(vtx);
  }

  return {std::move(fittedParticles), std::move(fittedVertices)};
}

// ---------------------------------------------------------------------------
DECLARE_COMPONENT(DTFAlgorithm)

}  // namespace k4DTF