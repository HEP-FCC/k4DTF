// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/TrackNode.h"

#include <cmath>

#include "core/FitParams.h"
#include "core/Projection.h"

namespace k4dtf {

TrackNode::TrackNode(const DecayNodeInput& input, const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  const auto& meas = std::get<TrackMeasurement>(input.measurement);
  m_momentum = meas.momentum;
  m_momCov = meas.cov;
  setCharge(static_cast<int>(meas.charge));
}

ErrCode TrackNode::initPar1(FitParams* /*fitparams*/) {
  return ErrCode::success;
}

ErrCode TrackNode::initPar2(FitParams* fitparams) {
  int momindex = momIndex();
  fitparams->setPar(momindex + 0, m_momentum(0));
  fitparams->setPar(momindex + 1, m_momentum(1));
  fitparams->setPar(momindex + 2, m_momentum(2));
  return ErrCode::success;
}

ErrCode TrackNode::initCov(FitParams* fitparams) const {
  int momindex = momIndex();
  const double factor = 1000;
  for (int row = 0; row < 3; ++row) {
    fitparams->cov()(momindex + row, momindex + row) = factor * m_momCov(row, row);
  }
  return ErrCode::success;
}

ErrCode TrackNode::projectTrackConstraint(const FitParams& fitparams,
                                           Projection& p) const {
  int momindex = momIndex();
  for (int i = 0; i < 3; ++i) {
    p.r(i) = fitparams.par(momindex + i) - m_momentum(i);
    p.H(i, momindex + i) = 1;
    for (int j = 0; j <= i; ++j) {
      p.V(i, j) = m_momCov(i, j);
    }
  }
  return ErrCode::success;
}

ErrCode TrackNode::filter(const FitParams& reference, FitParams& params) const {
  return applyConstraint(reference, params, 3, &TrackNode::projectTrackConstraint);
}

} // namespace k4dtf
