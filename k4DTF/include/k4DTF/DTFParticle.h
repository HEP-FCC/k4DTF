// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <string>

#include "Math/Vector4D.h"
#include "edm4hep/ReconstructedParticleData.h"

namespace k4DTF {

/// Internal particle representation used by the fitter.
///
/// Covariance matrix convention: upper triangle of the 4×4 (px,py,pz,E)
/// covariance stored row-major → 10 elements, matching the layout of
/// edm4hep::ReconstructedParticleData::covMatrix.
struct DTFParticle {
  ROOT::Math::PxPyPzEVector p4;
  std::array<double, 10> cov{};  // (px,py,pz,E) covariance, upper triangle
  double charge{0.0};
  int pdg{0};
  std::string name;

  /// Invariant mass from the 4-momentum.
  double mass() const { return p4.M(); }

  /// Construct from an EDM4hep reconstructed particle data struct.
  static DTFParticle fromEDM4hep(const edm4hep::ReconstructedParticleData& rp);

  /// Convert back to an EDM4hep reconstructed particle data struct.
  edm4hep::ReconstructedParticleData toEDM4hep() const;
};

}  // namespace k4DTF
