// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/TrackNode.h"

#include <algorithm>
#include <cmath>

#include "core/FitParams.h"
#include "core/Projection.h"
#include "k4DTF/core/Helix.h"

namespace k4dtf {

TrackNode::TrackNode(const DecayNodeInput& input, const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  const auto& meas = std::get<TrackMeasurement>(input.measurement);
  m_helix = meas.helix;
  if (m_helix) {
    m_momentum = helix::momentum(m_helix->par, m_helix->Bz);
    m_momCov = helix::momentumCovariance(m_helix->par, m_helix->cov, m_helix->Bz);
    const double q = meas.charge != 0 ? meas.charge : helix::charge(m_helix->par, m_helix->Bz);
    setCharge(static_cast<int>(std::lround(q)));
  } else {
    m_momentum = meas.momentum;
    m_momCov = meas.cov;
    setCharge(static_cast<int>(std::lround(meas.charge)));
  }
  setPdtMass(input.mass);
}

ErrCode TrackNode::initPar1(FitParams* /*fitparams*/) {
  // the energy of a track is computed from its mass hypothesis
  if (pdtMass() < 0) return ErrCode::badsetup;
  if (m_helix && (charge() == 0 || m_helix->Bz == 0)) return ErrCode::badsetup;
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

ErrCode TrackNode::projectHelixConstraint(const FitParams& fitparams, Projection& p) const {
  const int posindex = mother() ? mother()->posIndex() : -1;
  if (posindex < 0) return ErrCode::badsetup;
  const int momindex = momIndex();

  // parameters this constraint depends on: production vertex and momentum
  Eigen::Matrix<double, 6, 1> q;
  for (int i = 0; i < 3; ++i) {
    q(i) = fitparams.par(posindex + i);
    q(3 + i) = fitparams.par(momindex + i);
  }
  const int index[6] = {posindex, posindex + 1, posindex + 2, momindex, momindex + 1, momindex + 2};
  const double Bz = m_helix->Bz;
  const double charge = this->charge();
  const auto& ref = m_helix->referencePoint;
  auto h = [&](const Eigen::Matrix<double, 6, 1>& v) {
    return helix::fromPositionMomentum(v.head<3>(), v.tail<3>(), charge, Bz, ref);
  };

  helix::Vector5d r = h(q) - m_helix->par;
  r(helix::Phi) = helix::wrapPhi(r(helix::Phi));

  // Jacobian by central differences.
  // TODO: analytic derivatives (see the BaBar/LHCb DTF helix constraint).
  const double pmag = q.tail<3>().norm();
  for (int k = 0; k < 6; ++k) {
    const double step = k < 3 ? 1e-4 * std::max(1.0, std::abs(q(k))) : 1e-6 * std::max(1e-3, pmag);
    auto up = q, dn = q;
    up(k) += step;
    dn(k) -= step;
    helix::Vector5d d = h(up) - h(dn);
    d(helix::Phi) = helix::wrapPhi(d(helix::Phi));
    d /= 2 * step;
    for (int row = 0; row < 5; ++row) p.H(row, index[k]) = d(row);
  }
  for (int row = 0; row < 5; ++row) {
    p.r(row) = r(row);
    for (int col = 0; col <= row; ++col) p.V(row, col) = m_helix->cov(row, col);
  }
  return ErrCode::success;
}

ErrCode TrackNode::filter(const FitParams& reference, FitParams& params) const {
  if (m_helix) return applyConstraint(reference, params, 5, &TrackNode::projectHelixConstraint);
  return applyConstraint(reference, params, 3, &TrackNode::projectTrackConstraint);
}

} // namespace k4dtf
