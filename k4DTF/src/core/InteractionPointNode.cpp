// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/InteractionPointNode.h"

#include "core/FitParams.h"
#include "core/Projection.h"

namespace k4dtf {

InteractionPointNode::InteractionPointNode(const Eigen::Vector3d& pos,
                                           const Eigen::Matrix3d& cov)
    : ParticleNode("IP", 0, nullptr), m_ipPos(pos), m_ipCov(cov) {}

ErrCode InteractionPointNode::initPar1(FitParams* fitparams) {
  ErrCode status = ErrCode::success;
  int posindex = posIndex();
  for (int row = 0; row < 3; ++row) {
    fitparams->setPar(posindex + row, m_ipPos(row));
  }
  for (auto* dau : daughters()) {
    status |= dau->initPar1(fitparams);
    status |= dau->initPar2(fitparams);
  }
  return status;
}

ErrCode InteractionPointNode::initPar2(FitParams* /*fitparams*/) {
  return ErrCode::success;
}

ErrCode InteractionPointNode::initCov(FitParams* fitparams) const {
  ErrCode status = ErrCode::success;
  int posindex = posIndex();
  for (int row = 0; row < 3; ++row) {
    fitparams->cov()(posindex + row, posindex + row) = 1000 * m_ipCov(row, row);
  }
  for (auto* dau : daughters()) {
    status |= dau->initCov(fitparams);
  }
  return status;
}

ErrCode InteractionPointNode::projectIPConstraint(const FitParams& fitparams,
                                                    Projection& p) const {
  int posindex = posIndex();
  for (int row = 0; row < 3; ++row) {
    p.r(row) = fitparams.par(posindex + row) - m_ipPos(row);
    p.H(row, posindex + row) = 1;
    for (int col = 0; col <= row; ++col) {
      p.V(row, col) = m_ipCov(row, col);
    }
  }
  return ErrCode::success;
}

ErrCode InteractionPointNode::filter(const FitParams& reference, FitParams& params) const {
  for (auto* dau : daughters()) {
    ErrCode status = dau->filter(reference, params);
    if (failure(status)) return status;
  }
  return applyConstraint(reference, params, 3, &InteractionPointNode::projectIPConstraint);
}

} // namespace k4dtf
