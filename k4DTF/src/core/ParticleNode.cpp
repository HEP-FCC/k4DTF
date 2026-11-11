// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/ParticleNode.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "core/CompositeNode.h"
#include "core/FitParams.h"
#include "core/InternalParticleNode.h"
#include "core/InteractionPointNode.h"
#include "core/KalmanCalculator.h"
#include "core/MergedPi0Node.h"
#include "core/PhotonNode.h"
#include "core/Projection.h"
#include "core/ResonanceNode.h"
#include "core/TrackNode.h"

namespace k4DTF {

ParticleNode::ParticleNode(const std::string& name, int pdg, const ParticleNode* mother)
    : m_name(name), m_pdgCode(pdg), m_mother(mother) {}

ParticleNode::~ParticleNode() = default;

std::unique_ptr<ParticleNode> ParticleNode::create(const DecayNodeInput& input,
                                                    const ParticleNode* mother) {
  bool isComposite = !input.daughters.empty();

  if (!mother) {
    return std::make_unique<InternalParticleNode>(input, nullptr);
  }

  // Composites with a lifetime get their own decay vertex and a flight length
  // from the mother's vertex (or from the interaction point, for the head);
  // resonances share the vertex of their mother.
  if (isComposite) {
    if (input.isResonance) return std::make_unique<ResonanceNode>(input, mother);
    return std::make_unique<InternalParticleNode>(input, mother);
  }

  return std::visit(
      [&](auto& meas) -> std::unique_ptr<ParticleNode> {
        using T = std::decay_t<decltype(meas)>;
        if constexpr (std::is_same_v<T, TrackMeasurement>) {
          return std::make_unique<TrackNode>(input, mother);
        } else if constexpr (std::is_same_v<T, PhotonMeasurement>) {
          return std::make_unique<PhotonNode>(input, mother);
        } else if constexpr (std::is_same_v<T, CompositeMeasurement>) {
          return std::make_unique<CompositeNode>(input, mother);
        } else if constexpr (std::is_same_v<T, ResonanceMeasurement>) {
          return std::make_unique<ResonanceNode>(input, mother);
        } else if constexpr (std::is_same_v<T, MergedPi0Measurement>) {
          return std::make_unique<MergedPi0Node>(input, mother);
        }
      },
      input.measurement);
}

void ParticleNode::addDaughter(std::unique_ptr<ParticleNode> dau) {
  dau->setMother(this);
  m_daughters.push_back(dau.get());
  m_daughters_owned.push_back(std::move(dau));
}

void ParticleNode::updateIndex(int& offset) {
  for (auto* dau : m_daughters) {
    dau->updateIndex(offset);
  }
  m_index = offset;
  offset += dim();
}

ErrCode ParticleNode::initCov(FitParams* fitparams) const {
  ErrCode status = ErrCode::success;
  for (auto* dau : m_daughters) {
    status |= dau->initCov(fitparams);
  }

  // flight length: wide prior, like LHCb [mm]
  int lenindex = lenIndex();
  if (lenindex >= 0) {
    const double sigl = 1000.0;
    fitparams->cov()(lenindex, lenindex) = sigl * sigl;
  }

  int posindex = posIndex();
  if (posindex >= 0 && hasPosition()) {
    const double sigx = 10.0;
    const double sigy = 10.0;
    const double sigz = 50.0;
    fitparams->cov()(posindex + 0, posindex + 0) = sigx * sigx;
    fitparams->cov()(posindex + 1, posindex + 1) = sigy * sigy;
    fitparams->cov()(posindex + 2, posindex + 2) = sigz * sigz;
  }

  int momindex = momIndex();
  if (momindex >= 0) {
    if (m_daughters.empty()) {
      const double sigmom = 10.0;
      int maxrow = hasEnergy() ? 4 : 3;
      for (int row = 0; row < maxrow; ++row) {
        fitparams->cov()(momindex + row, momindex + row) = sigmom * sigmom;
      }
    } else {
      for (int row = 0; row < 3; ++row) {
        fitparams->cov()(momindex + row, momindex + row) = 0;
      }
      for (auto* dau : m_daughters) {
        int daumomindex = dau->momIndex();
        for (int row = 0; row < 3; ++row) {
          fitparams->cov()(momindex + row, momindex + row) +=
              fitparams->cov()(daumomindex + row, daumomindex + row);
        }
      }
      if (hasEnergy()) {
        fitparams->cov()(momindex + 3, momindex + 3) =
            fitparams->cov()(momindex + 2, momindex + 2);
      }
    }
  }

  return status;
}

ErrCode ParticleNode::projectGeoConstraint(const FitParams& fitparams, Projection& p) const {
  int posindexmother = mother()->posIndex();
  int posindex = posIndex();
  int lenindex = lenIndex();
  int momindex = momIndex();

  double len = fitparams.par(lenindex);

  double px = fitparams.par(momindex + 0);
  double py = fitparams.par(momindex + 1);
  double pz = fitparams.par(momindex + 2);
  double p2 = px * px + py * py + pz * pz;
  double mom = std::sqrt(p2);

  for (int row = 0; row < 3; ++row) {
    double posxmother = fitparams.par(posindexmother + row);
    double posx = fitparams.par(posindex + row);
    double momx = fitparams.par(momindex + row);
    p.r(row) = posxmother - posx + len * momx / mom;
    p.H(row, posindexmother + row) = 1;
    p.H(row, posindex + row) = -1;
    p.H(row, lenindex) = momx / mom;
  }

  p.H(0, momindex + 0) = len / mom * (1 - px * px / p2);
  p.H(0, momindex + 1) = len / mom * (0 - px * py / p2);
  p.H(0, momindex + 2) = len / mom * (0 - px * pz / p2);
  p.H(1, momindex + 0) = len / mom * (0 - py * px / p2);
  p.H(1, momindex + 1) = len / mom * (1 - py * py / p2);
  p.H(1, momindex + 2) = len / mom * (0 - py * pz / p2);
  p.H(2, momindex + 0) = len / mom * (0 - pz * px / p2);
  p.H(2, momindex + 1) = len / mom * (0 - pz * py / p2);
  p.H(2, momindex + 2) = len / mom * (1 - pz * pz / p2);

  return ErrCode::success;
}

ErrCode ParticleNode::projectMassConstraint(const FitParams& fitparams, Projection& p) const {
  double mass = pdtMass();
  double mass2 = mass * mass;
  int momindex = momIndex();

  double px = fitparams.par(momindex + 0);
  double py = fitparams.par(momindex + 1);
  double pz = fitparams.par(momindex + 2);
  double E = fitparams.par(momindex + 3);
  p.r(0) = E * E - px * px - py * py - pz * pz - mass2;

  p.H(0, momindex + 0) = -2.0 * px;
  p.H(0, momindex + 1) = -2.0 * py;
  p.H(0, momindex + 2) = -2.0 * pz;
  p.H(0, momindex + 3) = 2.0 * E;

  return ErrCode::success;
}

ErrCode ParticleNode::initTau(FitParams* fitparams) const {
  int lenindex = lenIndex();
  if (lenindex >= 0 && hasPosition()) {
    const ParticleNode* amother = mother();
    if (amother && amother->posIndex() >= 0) {
      int momposindex = amother->posIndex();
      int posindex = posIndex();
      int momindex = momIndex();
      double momdX = 0;
      double mom2 = 0;
      for (int irow = 0; irow < 3; ++irow) {
        double dx = fitparams->par(posindex + irow) - fitparams->par(momposindex + irow);
        double px = fitparams->par(momindex + irow);
        momdX += px * dx;
        mom2 += px * px;
      }
      double len = momdX / std::sqrt(mom2);
      fitparams->setPar(lenindex, len);
    }
  }
  return ErrCode::success;
}

void ParticleNode::collectVertexDaughters(ParticleVector& particles, int posindex) {
  if (mother() && mother()->posIndex() == posindex) particles.push_back(this);
  for (auto* dau : m_daughters) {
    dau->collectVertexDaughters(particles, posindex);
  }
}

ChiSquare ParticleNode::chiSquare(const FitParams& params) const {
  ChiSquare chi2;
  for (auto* dau : m_daughters) {
    chi2 += dau->chiSquare(params);
  }
  chi2 += params.chiSquare(*this);
  chi2 += ChiSquare(0, -dim());
  return chi2;
}

int ParticleNode::nFinalChargedCandidates() const {
  int rc = 0;
  for (auto* dau : m_daughters) {
    rc += dau->nFinalChargedCandidates();
  }
  return rc;
}

} // namespace k4DTF
