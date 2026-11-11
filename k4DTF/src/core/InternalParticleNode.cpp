// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/InternalParticleNode.h"

#include <cmath>
#include <iostream>

#include "core/FitParams.h"
#include "core/Projection.h"
#include "core/ResonanceNode.h"
#include "core/TrackNode.h"
#include "k4DTF/core/Helix.h"

#include <vector>

namespace k4DTF {

InternalParticleNode::InternalParticleNode(const DecayNodeInput& input,
                                           const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  for (const auto& dauInput : input.daughters) {
    auto dau = ParticleNode::create(dauInput, this);
    addDaughter(std::move(dau));
  }
  setPdtMass(input.mass);
  if (input.applyMassConstraint) {
    setPdtMass(input.massConstraintValue > 0 ? input.massConstraintValue : input.mass);
    setHasMassConstraint(pdtMass() > 0);
  }
}

ErrCode InternalParticleNode::initPar1(FitParams* fitparams) {
  ErrCode status = ErrCode::success;
  int posindex = posIndex();

  for (int row = 0; row < 3; ++row) {
    fitparams->setPar(posindex + row, 0);
  }

  for (auto* dau : daughters()) {
    status |= dau->initPar1(fitparams);
  }

  seedVertex(fitparams);

  for (int row = 0; row < 3; ++row) {
    if (fitparams->par(posindex + row) != 0) break;
    if (row == 2) {
      if (mother() && mother()->posIndex() >= 0) {
        int posindexmom = mother()->posIndex();
        for (int irow = 0; irow < 3; ++irow) {
          fitparams->setPar(posindex + irow, fitparams->par(posindexmom + irow));
        }
      }
    }
  }

  for (auto* dau : daughters()) {
    status |= dau->initPar2(fitparams);
  }

  initMom(fitparams);

  return status;
}

namespace {

/// Helix tracks produced at the vertex of `node`: its track daughters and
/// those of its resonance daughters, which decay at the same vertex.
void collectHelices(const ParticleNode* node, std::vector<const HelixMeasurement*>& out) {
  for (const auto* dau : node->daughters()) {
    if (const auto* t = dynamic_cast<const TrackNode*>(dau)) {
      if (const auto* h = t->helixMeasurement()) out.push_back(h);
    } else if (dynamic_cast<const ResonanceNode*>(dau)) {
      collectHelices(dau, out);
    }
  }
}

}  // namespace

void InternalParticleNode::seedVertex(FitParams* fitparams) const {
  std::vector<const HelixMeasurement*> helices;
  collectHelices(this, helices);
  if (helices.size() < 2) return;

  // the pair crossing at the largest angle gives the best-defined point
  std::size_t best1 = 0, best2 = 1;
  double bestSin = -1;
  for (std::size_t i = 0; i < helices.size(); ++i)
    for (std::size_t j = i + 1; j < helices.size(); ++j) {
      const double s = std::abs(std::sin(helices[i]->par(helix::Phi) - helices[j]->par(helix::Phi)));
      if (s > bestSin) {
        bestSin = s;
        best1 = i;
        best2 = j;
      }
    }

  // candidates compatible within 5 sigma of the impact parameters
  const auto& h1 = *helices[best1];
  const auto& h2 = *helices[best2];
  const double sigma = std::sqrt(h1.cov(helix::D0, helix::D0) + h1.cov(helix::Z0, helix::Z0) +
                                 h2.cov(helix::D0, helix::D0) + h2.cov(helix::Z0, helix::Z0));
  Eigen::Vector3d vertex;
  if (!helix::intersect(h1.par, h1.referencePoint, h2.par, h2.referencePoint, vertex, 5 * sigma)) return;
  const int posindex = posIndex();
  for (int row = 0; row < 3; ++row) fitparams->setPar(posindex + row, vertex(row));
}

ErrCode InternalParticleNode::initPar2(FitParams* fitparams) {
  int posindex = posIndex();
  if (hasPosition() && mother() && mother()->posIndex() >= 0) {
    bool allZero = true;
    for (int irow = 0; irow < 3; ++irow) {
      if (fitparams->par(posindex + irow) != 0) {
        allZero = false;
        break;
      }
    }
    if (allZero) {
      int posindexmom = mother()->posIndex();
      for (int irow = 0; irow < 3; ++irow) {
        fitparams->setPar(posindex + irow, fitparams->par(posindexmom + irow));
      }
    }
  }
  return initTau(fitparams);
}

ErrCode InternalParticleNode::initMom(FitParams* fitparams) const {
  int momindex = momIndex();
  for (int irow = 0; irow < 4; ++irow) {
    fitparams->setPar(momindex + irow, 0);
  }

  for (auto* dau : daughters()) {
    int daumomindex = dau->momIndex();
    double e2 = 0;
    int maxrow = dau->hasEnergy() ? 4 : 3;
    for (int irow = 0; irow < maxrow; ++irow) {
      double px = fitparams->par(daumomindex + irow);
      e2 += px * px;
      fitparams->setPar(momindex + irow, fitparams->par(momindex + irow) + px);
    }
    if (maxrow == 3) {
      double mass = dau->pdtMass();
      fitparams->setPar(momindex + 3, fitparams->par(momindex + 3) + std::sqrt(e2 + mass * mass));
    }
  }

  return ErrCode::success;
}

ErrCode InternalParticleNode::projectKineConstraint(const FitParams& fitparams,
                                                     Projection& p) const {
  int momindex = momIndex();
  for (int imom = 0; imom < 4; ++imom) {
    p.r(imom) = fitparams.par(momindex + imom);
    p.H(imom, momindex + imom) = 1;
  }

  for (auto* dau : daughters()) {
    int daumomindex = dau->momIndex();
    double mass = dau->pdtMass();
    double e2 = mass * mass;
    int maxrow = dau->hasEnergy() ? 4 : 3;
    for (int imom = 0; imom < maxrow; ++imom) {
      double px = fitparams.par(daumomindex + imom);
      e2 += px * px;
      p.r(imom) += -px;
      p.H(imom, daumomindex + imom) = -1;
    }
    if (maxrow == 3) {
      double energy = std::sqrt(e2);
      p.r(3) += -energy;
      for (int jmom = 0; jmom < 3; ++jmom) {
        double px = fitparams.par(daumomindex + jmom);
        p.H(3, daumomindex + jmom) = -px / energy;
      }
    }
  }
  return ErrCode::success;
}

ErrCode InternalParticleNode::filter(const FitParams& reference, FitParams& params) const {
  for (auto* dau : daughters()) {
    ErrCode status = dau->filter(reference, params);
    if (failure(status)) return status;
  }

  ErrCode status = applyConstraint(reference, params, 4,
                                   &InternalParticleNode::projectKineConstraint);
  if (mother() && lenIndex() >= 0) {
    status |= applyConstraint(reference, params, 3,
                              &ParticleNode::projectGeoConstraint);
  }
  if (hasMassConstraint()) {
    status |= applyConstraint(reference, params, 1,
                              &ParticleNode::projectMassConstraint);
  }
  return status;
}

} // namespace k4DTF
