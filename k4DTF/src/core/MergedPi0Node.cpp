// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/MergedPi0Node.h"

#include <cmath>

#include "core/FitParams.h"
#include "core/Projection.h"

namespace k4dtf {

MergedPi0Node::MergedPi0Node(const DecayNodeInput& input, const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  const auto& meas = std::get<MergedPi0Measurement>(input.measurement);
  m_p4Cov = meas.cov;
  setHasMassConstraint(input.applyMassConstraint);
  for (const auto& cl : meas.clusters) {
    m_clusterPositions.push_back(cl.pos);
    m_clusterEnergies.push_back(cl.energy);
  }
  m_p4.setZero();
  for (size_t i = 0; i < m_clusterPositions.size(); ++i) {
    m_p4(3) += m_clusterEnergies[i];
  }
}

ErrCode MergedPi0Node::initPar1(FitParams* fitparams) {
  return initPar2(fitparams);
}

ErrCode MergedPi0Node::initPar2(FitParams* fitparams) {
  int momindex = momIndex();
  fitparams->setPar(momindex + 0, m_p4(0));
  fitparams->setPar(momindex + 1, m_p4(1));
  fitparams->setPar(momindex + 2, m_p4(2));
  if (hasEnergy()) {
    fitparams->setPar(momindex + 3, m_p4(3));
  }
  return ErrCode::success;
}

ErrCode MergedPi0Node::initCov(FitParams* fitparams) const {
  int momindex = momIndex();
  int thisdim = dim();
  const double factor = 1000;
  for (int row = 0; row < thisdim; ++row) {
    fitparams->cov()(momindex + row, momindex + row) = factor * m_p4Cov(row, row);
  }
  return ErrCode::success;
}

ErrCode MergedPi0Node::projectPi0Constraint(const FitParams& fitparams,
                                              Projection& p) const {
  const int posindex = mother() ? mother()->posIndex() : -1;
  if (posindex < 0) return ErrCode::badsetup;

  Eigen::Vector3d motherpos(fitparams.par(posindex + 0),
                             fitparams.par(posindex + 1),
                             fitparams.par(posindex + 2));

  Eigen::Vector4d p4calc = Eigen::Vector4d::Zero();
  for (size_t ic = 0; ic < m_clusterPositions.size(); ++ic) {
    Eigen::Vector3d dx = m_clusterPositions[ic] - motherpos;
    double L2 = dx.squaredNorm();
    if (L2 > 0) {
      double L = std::sqrt(L2);
      p4calc(0) += m_clusterEnergies[ic] * dx.x() / L;
      p4calc(1) += m_clusterEnergies[ic] * dx.y() / L;
      p4calc(2) += m_clusterEnergies[ic] * dx.z() / L;
    }
    p4calc(3) += m_clusterEnergies[ic];
  }

  const int momindex = momIndex();
  p.r(0) = p4calc(0) - fitparams.par(momindex + 0);
  p.r(1) = p4calc(1) - fitparams.par(momindex + 1);
  p.r(2) = p4calc(2) - fitparams.par(momindex + 2);
  if (hasEnergy()) {
    p.r(3) = p4calc(3) - fitparams.par(momindex + 3);
  }

  for (size_t ic = 0; ic < m_clusterPositions.size(); ++ic) {
    Eigen::Vector3d dx = m_clusterPositions[ic] - motherpos;
    double dx2 = dx.x() * dx.x();
    double dy2 = dx.y() * dx.y();
    double dz2 = dx.z() * dx.z();
    double L2 = dx2 + dy2 + dz2;
    double EOverL3 = m_clusterEnergies[ic] / (L2 * std::sqrt(L2));
    p.H(0, posindex + 0) += -(dy2 + dz2) * EOverL3;
    p.H(1, posindex + 0) += dx.y() * dx.x() * EOverL3;
    p.H(2, posindex + 0) += dx.z() * dx.x() * EOverL3;
    p.H(0, posindex + 1) += dx.x() * dx.y() * EOverL3;
    p.H(1, posindex + 1) += -(dx2 + dz2) * EOverL3;
    p.H(2, posindex + 1) += dx.z() * dx.y() * EOverL3;
    p.H(0, posindex + 2) += dx.x() * dx.z() * EOverL3;
    p.H(1, posindex + 2) += dx.y() * dx.z() * EOverL3;
    p.H(2, posindex + 2) += -(dx2 + dy2) * EOverL3;
  }

  int maxrow = dimM();
  for (int irow = 0; irow < maxrow; ++irow) {
    p.H(irow, momindex + irow) = -1;
  }

  for (int irow = 0; irow < maxrow; ++irow) {
    for (int icol = 0; icol <= irow; ++icol) {
      p.V(irow, icol) = m_p4Cov(irow, icol);
    }
  }

  return ErrCode::success;
}

ErrCode MergedPi0Node::filter(const FitParams& reference, FitParams& params) const {
  return applyConstraint(reference, params, dimM(),
                         &MergedPi0Node::projectPi0Constraint);
}

} // namespace k4dtf
