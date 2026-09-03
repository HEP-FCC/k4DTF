// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/PhotonNode.h"

#include <cmath>

#include "core/FitParams.h"
#include "core/Projection.h"

namespace k4dtf {

PhotonNode::PhotonNode(const DecayNodeInput& input, const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  const auto& meas = std::get<PhotonMeasurement>(input.measurement);
  m_clusterPos = meas.clusterPos;
  m_clusterEnergy = meas.energy;
  m_clusterCov = meas.cov;
}

ErrCode PhotonNode::initPar1(FitParams* fitparams) {
  int momindex = momIndex();
  Eigen::Vector3d dx = m_clusterPos;
  double l = dx.norm();
  if (l > 0) {
    fitparams->setPar(momindex + 0, m_clusterEnergy * dx.x() / l);
    fitparams->setPar(momindex + 1, m_clusterEnergy * dx.y() / l);
    fitparams->setPar(momindex + 2, m_clusterEnergy * dx.z() / l);
  }
  return ErrCode::success;
}

ErrCode PhotonNode::initPar2(FitParams* fitparams) {
  if (mother() && mother()->posIndex() >= 0) {
    int posindex = mother()->posIndex();
    Eigen::Vector3d motherpos(fitparams->par(posindex + 0),
                               fitparams->par(posindex + 1),
                               fitparams->par(posindex + 2));
    int momindex = momIndex();
    Eigen::Vector3d dx = m_clusterPos - motherpos;
    double l = dx.norm();
    if (l > 0) {
      fitparams->setPar(momindex + 0, m_clusterEnergy * dx.x() / l);
      fitparams->setPar(momindex + 1, m_clusterEnergy * dx.y() / l);
      fitparams->setPar(momindex + 2, m_clusterEnergy * dx.z() / l);
    }
  }
  return ErrCode::success;
}

ErrCode PhotonNode::initCov(FitParams* fitparams) const {
  int momindex = momIndex();
  double varEnergy = m_clusterCov(2, 2);
  const double factor = 1000;
  for (int row = 0; row < 3; ++row) {
    fitparams->cov()(momindex + row, momindex + row) = factor * varEnergy;
  }
  return ErrCode::success;
}

ErrCode PhotonNode::projectPhotonConstraint(const FitParams& fitparams,
                                             Projection& p) const {
  int momindex = momIndex();
  double px = fitparams.par(momindex + 0);
  double py = fitparams.par(momindex + 1);
  double pz = fitparams.par(momindex + 2);
  double energy = std::sqrt(px * px + py * py + pz * pz);

  double tx = px / pz;
  double ty = py / pz;

  int posindex = mother() ? mother()->posIndex() : -1;
  if (posindex < 0) return ErrCode::badsetup;

  double dz = m_clusterPos.z() - fitparams.par(posindex + 2);
  double dx = m_clusterPos.x() - fitparams.par(posindex + 0) - dz * tx;
  double dy = m_clusterPos.y() - fitparams.par(posindex + 1) - dz * ty;
  double dE = m_clusterEnergy - energy;

  p.r(0) = dx;
  p.r(1) = dy;
  p.r(2) = dE;

  p.H(0, posindex + 0) = -1;
  p.H(1, posindex + 1) = -1;

  if (pz != 0) {
    p.H(0, momindex + 0) = -dz / pz;
    p.H(0, momindex + 2) = dz * px / (pz * pz);
    p.H(1, momindex + 1) = -dz / pz;
    p.H(1, momindex + 2) = dz * py / (pz * pz);
  }

  if (energy > 0) {
    p.H(2, momindex + 0) = -px / energy;
    p.H(2, momindex + 1) = -py / energy;
    p.H(2, momindex + 2) = -pz / energy;
  }

  for (int irow = 0; irow < 3; ++irow) {
    for (int icol = 0; icol <= irow; ++icol) {
      p.V(irow, icol) = m_clusterCov(irow, icol);
    }
  }

  return ErrCode::success;
}

ErrCode PhotonNode::filter(const FitParams& reference, FitParams& params) const {
  return applyConstraint(reference, params, 3, &PhotonNode::projectPhotonConstraint);
}

} // namespace k4dtf
