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

Eigen::Matrix<double, 5, 6> jacobian(const Eigen::Vector3d& xglobal, const Eigen::Vector3d& p, double q, double Bz,
                                     const Eigen::Vector3d& reference) {
  // Derivatives of the formulae in fromPositionMomentum(). With
  //   a = -q Bz kC,  U = px + a y,  W = py - a x,  T^2 = U^2 + W^2
  // the parameters are
  //   phi0 = atan2(W, U)            D0 = (T - pT) / a
  //   omega = -a / pT               tanLambda = pz / pT
  //   Z0 = z - pz s / a,  s = atan2(py, px) - phi0
  enum { X = 0, Y = 1, Zc = 2, PX = 3, PY = 4, PZ = 5 };
  Eigen::Matrix<double, 5, 6> J = Eigen::Matrix<double, 5, 6>::Zero();
  const Eigen::Vector3d x = xglobal - reference;
  const double px = p.x(), py = p.y(), pz = p.z();
  const double pt2 = px * px + py * py, pt = std::sqrt(pt2);
  const double a = -q * Bz * kC;
  if (pt <= 0.0 || a == 0.0) return J;

  const double U = px + a * x.y(), W = py - a * x.x();
  const double T2 = U * U + W * W, T = std::sqrt(T2);
  const double pt3 = pt2 * pt;

  // phi0
  J(Phi, X) = -a * U / T2;
  J(Phi, Y) = -a * W / T2;
  J(Phi, PX) = -W / T2;
  J(Phi, PY) = U / T2;

  // D0 = (T - pT) / a
  J(D0, X) = -W / T;
  J(D0, Y) = U / T;
  J(D0, PX) = (U / T - px / pt) / a;
  J(D0, PY) = (W / T - py / pt) / a;

  // omega = -a / pT
  J(Omega, PX) = a * px / pt3;
  J(Omega, PY) = a * py / pt3;

  // tanLambda = pz / pT
  J(TanLambda, PX) = -pz * px / pt3;
  J(TanLambda, PY) = -pz * py / pt3;
  J(TanLambda, PZ) = 1.0 / pt;

  // Z0 = z - pz s / a, with ds = d(atan2(py, px)) - dphi0
  const double s = wrapPhi(std::atan2(py, px) - std::atan2(W, U));
  J(Z0, Zc) = 1.0;
  J(Z0, X) = pz / a * J(Phi, X);
  J(Z0, Y) = pz / a * J(Phi, Y);
  J(Z0, PX) = -pz / a * (-py / pt2 - J(Phi, PX));
  J(Z0, PY) = -pz / a * (px / pt2 - J(Phi, PY));
  J(Z0, PZ) = -s / a;
  return J;
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
