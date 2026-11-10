// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/core/Helix.h"

#include <cmath>

namespace k4dtf::helix {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

double wrapPhi(double phi) {
  // no loop: stays finite-time for inf/nan input
  const double w = std::remainder(phi, 2 * kPi);
  return w <= -kPi ? w + 2 * kPi : w;
}

Vector5d fromPositionMomentum(const Eigen::Vector3d& xglobal, const Eigen::Vector3d& p, double q, double Bz,
                              const Eigen::Vector3d& reference) {
  // Delphes TrkUtil::XPtoPar, written in mm. a = -q Bz kC is the signed
  // curvature-to-momentum factor: omega = a / pT.
  Vector5d par = Vector5d::Zero();
  const Eigen::Vector3d x = xglobal - reference;
  const double pt = std::hypot(p.x(), p.y());
  const double a = -q * Bz * kC;
  if (pt <= 0.0 || a == 0.0) return par;

  const double r2 = x.x() * x.x() + x.y() * x.y();
  const double cross = x.x() * p.y() - x.y() * p.x();
  const double T = std::sqrt(pt * pt - 2 * a * cross + a * a * r2);
  const double phi0 = std::atan2(p.y() - a * x.x(), p.x() + a * x.y());
  // two algebraically equivalent forms; the second avoids cancellation
  double D = 0.0;
  if (r2 > 0.0) D = (std::abs(a) * std::sqrt(r2) > 0.1 * pt) ? (T - pt) / a : (-2 * cross + a * r2) / (T + pt);
  const double ct = p.z() / pt;
  const double C = a / (2 * pt);  // half curvature [1/mm]
  const double s = wrapPhi(std::atan2(p.y(), p.x()) - phi0);
  const double z0 = x.z() - ct * s / (2 * C);

  par(D0) = D;
  par(Phi) = phi0;
  par(Omega) = -2.0 * C;
  par(Z0) = z0;
  par(TanLambda) = ct;
  return par;
}

Eigen::Vector3d momentum(const Vector5d& par, double Bz) {
  const double omega = std::abs(par(Omega));
  if (omega == 0.0) return Eigen::Vector3d::Zero();
  const double pt = kC * std::abs(Bz) / omega;
  return {pt * std::cos(par(Phi)), pt * std::sin(par(Phi)), pt * par(TanLambda)};
}

Eigen::Vector3d pca(const Vector5d& par, const Eigen::Vector3d& reference) {
  return reference + Eigen::Vector3d(-par(D0) * std::sin(par(Phi)), par(D0) * std::cos(par(Phi)), par(Z0));
}

double charge(const Vector5d& par, double Bz) {
  if (par(Omega) == 0.0 || Bz == 0.0) return 0.0;
  return (par(Omega) > 0) == (Bz > 0) ? 1.0 : -1.0;
}

Eigen::Matrix<double, 3, 5> momentumJacobian(const Vector5d& par, double Bz) {
  Eigen::Matrix<double, 3, 5> J = Eigen::Matrix<double, 3, 5>::Zero();
  const double omega = par(Omega);
  if (omega == 0.0) return J;
  const Eigen::Vector3d p = momentum(par, Bz);
  const double pt = std::hypot(p.x(), p.y());
  J(0, Phi) = -p.y();
  J(1, Phi) = p.x();
  J(0, Omega) = -p.x() / omega;
  J(1, Omega) = -p.y() / omega;
  J(2, Omega) = -p.z() / omega;
  J(2, TanLambda) = pt;
  return J;
}

Eigen::Matrix3d momentumCovariance(const Vector5d& par, const Matrix5d& cov, double Bz) {
  const auto J = momentumJacobian(par, Bz);
  return J * cov * J.transpose();
}

}  // namespace k4dtf::helix
