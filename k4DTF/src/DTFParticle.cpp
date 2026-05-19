// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DTFParticle.h"

namespace k4DTF {

DTFParticle DTFParticle::fromEDM4hep(const edm4hep::ReconstructedParticleData& rp) {
  DTFParticle dp;
  dp.p4.SetPxPyPzE(rp.momentum.x, rp.momentum.y, rp.momentum.z, rp.energy);
  dp.charge = rp.charge;
  dp.pdg    = rp.PDG;
  for (int i = 0; i < 10; ++i) dp.cov[i] = rp.covMatrix[i];
  return dp;
}

edm4hep::ReconstructedParticleData DTFParticle::toEDM4hep() const {
  edm4hep::ReconstructedParticleData rp{};
  rp.momentum  = {static_cast<float>(p4.Px()),
                  static_cast<float>(p4.Py()),
                  static_cast<float>(p4.Pz())};
  rp.energy    = static_cast<float>(p4.E());
  rp.mass      = static_cast<float>(mass());
  rp.charge    = static_cast<float>(charge);
  rp.PDG       = pdg;
  for (int i = 0; i < 10; ++i) rp.covMatrix[i] = static_cast<float>(cov[i]);
  return rp;
}

}  // namespace k4DTF
