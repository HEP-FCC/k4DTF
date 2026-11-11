// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/CalorimeterResolution.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace k4DTF {

double CalorimeterResolution::angularResolution() const {
  return towerSize / std::sqrt(smearTowerCenter ? 6.0 : 12.0);
}

double CalorimeterResolution::energyVariance(double E) const {
  const double e = std::max(E, 0.0);
  return constant * constant * e * e + stochastic * stochastic * e + noise * noise;
}

Eigen::Matrix3d CalorimeterResolution::positionCovariance(const Eigen::Vector3d& pos) const {
  const double r = pos.norm(), rt = std::hypot(pos.x(), pos.y());
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  if (r <= 0) return cov;
  const double ct = pos.z() / r, st = rt / r;
  const double cp = rt > 0 ? pos.x() / rt : 1.0, sp = rt > 0 ? pos.y() / rt : 0.0;
  const Eigen::Vector3d uTheta(ct * cp, ct * sp, -st), uPhi(-sp, cp, 0.0);
  const double s = angularResolution();
  const double sTheta = r * s * st;  // r sigma(theta), sigma(theta) = sigma(eta) sin(theta)
  const double sPhi = rt * s;        // r sin(theta) sigma(phi)
  cov = sTheta * sTheta * uTheta * uTheta.transpose() + sPhi * sPhi * uPhi * uPhi.transpose();
  return cov;
}

Eigen::Vector3d CalorimeterResolution::frontFace(const Eigen::Vector3d& direction) const {
  const double rt = std::hypot(direction.x(), direction.y());
  const double inf = std::numeric_limits<double>::infinity();
  const double tBarrel = rt > 0 ? radius / rt : inf;
  const double tEndcap = direction.z() != 0 ? halfLength / std::abs(direction.z()) : inf;
  const double t = std::min(tBarrel, tEndcap);
  return std::isfinite(t) ? Eigen::Vector3d(t * direction) : Eigen::Vector3d::Zero();
}

PhotonMeasurement CalorimeterResolution::measurement(const Eigen::Vector3d& position, double E) const {
  PhotonMeasurement m;
  m.clusterPos = position;
  m.energy = E;
  m.cov.setZero();
  m.cov.topLeftCorner<3, 3>() = positionCovariance(position);
  m.cov(3, 3) = energyVariance(E);
  return m;
}

}  // namespace k4DTF
