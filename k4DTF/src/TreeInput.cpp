// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/TreeInput.h"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include "k4DTF/ParticleTable.h"

namespace k4DTF {

namespace {

bool constrained(int pdg, const std::vector<int>& list) {
  return std::any_of(list.begin(), list.end(), [&](int c) { return c == pdg || c == ParticleTable::chargeConjugate(pdg); });
}

k4DTF::DecayNodeInput build(const DecayDescriptor& desc, int node, const std::vector<DTFParticle>& leaves,
                            const std::vector<int>& massConstraints) {
  const auto& f = desc.flatten().at(node);
  k4DTF::DecayNodeInput in;
  in.pdgCode = f.pdg;
  in.name = f.name;
  in.mass = f.mass;

  if (f.isLeaf()) {
    const auto& p = leaves.at(f.leafIndex);
    if (in.mass < 0) in.mass = p.mass();
    if (p.hasCluster() && !p.hasTrack()) {
      if (p.cluster->energy <= 0)
        throw std::invalid_argument("makeDecayNodeInput: cluster of '" + f.name + "' has no energy");
      in.measurement = p.cluster->toPhotonMeasurement();
      return in;
    }
    if (p.charge == 0 && !p.hasTrack())
      throw std::invalid_argument("makeDecayNodeInput: neutral final-state particle '" + f.name +
                                  "' needs a calorimeter cluster (DTFParticle::cluster)");
    k4DTF::TrackMeasurement t;
    t.charge = p.charge;
    t.momentum = {p.p4.Px(), p.p4.Py(), p.p4.Pz()};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) t.cov(i, j) = p.cov[i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i];
    if (p.hasTrack()) t.helix = p.track->toHelixMeasurement();
    in.measurement = t;
    return in;
  }

  in.isResonance = ParticleTable::isResonance(f.pdg);
  in.applyMassConstraint = constrained(f.pdg, massConstraints);
  if (in.applyMassConstraint && in.mass <= 0)
    throw std::invalid_argument("makeDecayNodeInput: no nominal mass for constrained particle '" + f.name + "'");
  for (int c : f.children) in.daughters.push_back(build(desc, c, leaves, massConstraints));
  return in;
}

}  // namespace

k4DTF::DecayNodeInput makeDecayNodeInput(const DecayDescriptor& desc, const std::vector<DTFParticle>& leaves,
                                         const std::vector<int>& massConstraints) {
  if (leaves.size() != desc.nLeaves())
    throw std::invalid_argument("makeDecayNodeInput: " + std::to_string(leaves.size()) + " particles for " +
                                std::to_string(desc.nLeaves()) + " final-state particles of '" + desc.toString() +
                                "'");
  return build(desc, 0, leaves, massConstraints);
}

k4DTF::InteractionPointInput makeInteractionPoint(const DTFVertex& v) {
  k4DTF::InteractionPointInput ip;
  ip.position = {v.position.X(), v.position.Y(), v.position.Z()};
  const int idx[3][3] = {{0, 1, 3}, {1, 2, 4}, {3, 4, 5}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) ip.cov(i, j) = v.cov[idx[i][j]];
  return ip;
}

}  // namespace k4DTF
