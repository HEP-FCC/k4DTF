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
  setPdtMass(0);
}

void PhotonNode::initMomentum(FitParams* fitparams, const Eigen::Vector3d& origin) {
  const Eigen::Vector3d dx = m_clusterPos - origin;
  const double l = dx.norm();
  if (l <= 0) return;
  const Eigen::Vector3d p = m_clusterEnergy * dx / l;
  const int momindex = momIndex();
  for (int i = 0; i < 3; ++i) fitparams->setPar(momindex + i, p(i));

  // The running coordinate is fixed here, once per fit, so that all
  // iterations linearise the same constraint.
  p.cwiseAbs().maxCoeff(&m_i1);
  m_i2 = (m_i1 + 1) % 3;
  m_i3 = (m_i1 + 2) % 3;
}

ErrCode PhotonNode::initPar1(FitParams* fitparams) {
  initMomentum(fitparams, Eigen::Vector3d::Zero());
  return m_clusterEnergy > 0 ? ErrCode::success : ErrCode::badsetup;
}

ErrCode PhotonNode::initPar2(FitParams* fitparams) {
  // the mother's vertex is initialised by now: point from it to the cluster
  if (mother() && mother()->posIndex() >= 0) {
    const int posindex = mother()->posIndex();
    initMomentum(fitparams, Eigen::Vector3d(fitparams->par(posindex + 0), fitparams->par(posindex + 1),
                                            fitparams->par(posindex + 2)));
  }
  return ErrCode::success;
}

ErrCode PhotonNode::initCov(FitParams* fitparams) const {
  const int momindex = momIndex();
  const double varEnergy = m_clusterCov(3, 3);
  const double factor = 1000;
  for (int row = 0; row < 3; ++row) {
    fitparams->cov()(momindex + row, momindex + row) = factor * varEnergy;
  }
  return ErrCode::success;
}

ErrCode PhotonNode::projectPhotonConstraint(const FitParams& fitparams, Projection& p) const {
  const int posindex = mother() ? mother()->posIndex() : -1;
  if (posindex < 0) return ErrCode::badsetup;
  const int momindex = momIndex();

  const Eigen::Vector3d x(fitparams.par(posindex), fitparams.par(posindex + 1), fitparams.par(posindex + 2));
  const Eigen::Vector3d mom(fitparams.par(momindex), fitparams.par(momindex + 1), fitparams.par(momindex + 2));
  const double pabs = mom.norm();
  const double p1 = mom(m_i1);
  if (p1 == 0 || pabs == 0) return ErrCode::badsetup;

  const Eigen::Vector3d d = m_clusterPos - x;
  const int other[2] = {m_i2, m_i3};

  // dr/d(cluster x, y, z, E), for the covariance of the residual
  Eigen::Matrix<double, 3, 4> J = Eigen::Matrix<double, 3, 4>::Zero();

  for (int row = 0; row < 2; ++row) {
    const int i = other[row];
    const double slope = mom(i) / p1;
    p.r(row) = d(i) - d(m_i1) * slope;

    // r depends on the mother's vertex and on the photon momentum
    p.H(row, posindex + i) = -1;
    p.H(row, posindex + m_i1) = slope;
    p.H(row, momindex + i) = -d(m_i1) / p1;
    p.H(row, momindex + m_i1) = d(m_i1) * slope / p1;

    J(row, i) = 1;
    J(row, m_i1) = -slope;
  }

  p.r(2) = pabs - m_clusterEnergy;
  for (int k = 0; k < 3; ++k) p.H(2, momindex + k) = mom(k) / pabs;
  J(2, 3) = -1;

  const Eigen::Matrix3d V = J * m_clusterCov * J.transpose();
  for (int row = 0; row < 3; ++row)
    for (int col = 0; col <= row; ++col) p.V(row, col) = V(row, col);

  return ErrCode::success;
}

ErrCode PhotonNode::filter(const FitParams& reference, FitParams& params) const {
  return applyConstraint(reference, params, 3, &PhotonNode::projectPhotonConstraint);
}

} // namespace k4dtf
