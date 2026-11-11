// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/core/Helix.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

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

namespace {

/// Centre of the helix circle in the transverse plane.
Eigen::Vector2d centre(const Vector5d& par, const Eigen::Vector3d& reference) {
  const Eigen::Vector3d x0 = pca(par, reference);
  return {x0.x() + std::sin(par(Phi)) / par(Omega), x0.y() - std::cos(par(Phi)) / par(Omega)};
}

/// Arc length of the circle point in direction (point - centre).
double arcLengthTowards(const Vector5d& par, const Eigen::Vector2d& c, const Eigen::Vector2d& point) {
  const double omega = par(Omega);
  const Eigen::Vector2d v = point - c;
  const double phi = std::atan2(-omega * v.x(), omega * v.y());
  return -wrapPhi(phi - par(Phi)) / omega;
}

}  // namespace

Eigen::Vector3d positionAt(const Vector5d& par, double s, const Eigen::Vector3d& reference) {
  const double omega = par(Omega);
  if (omega == 0.0) {
    // straight line
    const Eigen::Vector3d x0 = pca(par, reference);
    return x0 + s * Eigen::Vector3d(std::cos(par(Phi)), std::sin(par(Phi)), par(TanLambda));
  }
  const Eigen::Vector2d c = centre(par, reference);
  const double phi = par(Phi) - omega * s;
  return {c.x() - std::sin(phi) / omega, c.y() + std::cos(phi) / omega,
          reference.z() + par(Z0) + par(TanLambda) * s};
}

Eigen::Vector3d momentumAt(const Vector5d& par, double Bz, double s) {
  Vector5d moved = par;
  moved(Phi) = wrapPhi(par(Phi) - par(Omega) * s);
  return momentum(moved, Bz);
}

double arcLengthNear(const Vector5d& par, const Eigen::Vector3d& point, const Eigen::Vector3d& reference) {
  if (par(Omega) == 0.0) {
    const Eigen::Vector3d x0 = pca(par, reference);
    return (point - x0).head<2>().dot(Eigen::Vector2d(std::cos(par(Phi)), std::sin(par(Phi))));
  }
  return arcLengthTowards(par, centre(par, reference), point.head<2>());
}

bool intersect(const Vector5d& par1, const Eigen::Vector3d& ref1, const Vector5d& par2, const Eigen::Vector3d& ref2,
               Eigen::Vector3d& vertex, double tolerance) {
  if (par1(Omega) == 0.0 || par2(Omega) == 0.0) return false;
  const Eigen::Vector2d c1 = centre(par1, ref1), c2 = centre(par2, ref2);
  const double r1 = 1.0 / std::abs(par1(Omega)), r2 = 1.0 / std::abs(par2(Omega));
  const Eigen::Vector2d d = c2 - c1;
  const double dist = d.norm();
  if (dist < 1e-9 * (r1 + r2)) return false;
  const Eigen::Vector2d u = d / dist;

  // candidate transverse points: (p1 on circle 1, p2 on circle 2)
  std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> cand;
  if (dist <= r1 + r2 && dist >= std::abs(r1 - r2)) {
    const double a = (dist * dist + r1 * r1 - r2 * r2) / (2 * dist);
    const double h = std::sqrt(std::max(0.0, r1 * r1 - a * a));
    const Eigen::Vector2d base = c1 + a * u, perp(-u.y(), u.x());
    cand.emplace_back(base + h * perp, base + h * perp);
    cand.emplace_back(base - h * perp, base - h * perp);
  } else if (dist > r1 + r2) {
    // apart: closest points on the line of centres
    cand.emplace_back(c1 + r1 * u, c2 - r2 * u);
  } else {
    // one circle inside the other: closest points on the same side
    const double side = r1 > r2 ? 1.0 : -1.0;
    cand.emplace_back(c1 + side * r1 * u, c2 + side * r2 * u);
  }
  // Nearly parallel tracks: their circles may not meet (measurement errors)
  // or meet far away. Also try the points closest to the reference point
  // (the beam line), where prompt and short-lived decays are.
  const Eigen::Vector2d r0 = 0.5 * (ref1 + ref2).head<2>();
  auto onCircle = [](const Eigen::Vector2d& c, double r, const Eigen::Vector2d& q) {
    const Eigen::Vector2d v = q - c;
    return v.norm() > 0 ? Eigen::Vector2d(c + r * v.normalized()) : c;
  };
  cand.emplace_back(onCircle(c1, r1, r0), onCircle(c2, r2, r0));

  // The helices should be close at the vertex. Nearly parallel tracks stay
  // close over a long stretch, and their circles cross again far away (tens
  // of cm). So among the candidates where the helices are within
  // `tolerance` of the closest one, take the one nearest the beam line.
  // A boosted pair (a KS of several GeV) can also have its second crossing
  // behind the beam line, as close as the real vertex: candidates more than
  // kBehind behind the beam line, along the pair's transverse momentum, are
  // only used when there is nothing else.
  constexpr double kBehind = 1.0;  // [mm]
  struct Point {
    double gap;
    bool behind;
    Eigen::Vector3d x;
  };
  std::vector<Point> points;
  double minGap = -1;
  for (const auto& [p1, p2] : cand) {
    const double s1 = arcLengthTowards(par1, c1, p1), s2 = arcLengthTowards(par2, c2, p2);
    const Eigen::Vector3d x1 = positionAt(par1, s1, ref1), x2 = positionAt(par2, s2, ref2);
    const Eigen::Vector3d x = 0.5 * (x1 + x2);
    // transverse momentum of the pair at the candidate, up to the common factor kC |Bz|
    const double phi1 = par1(Phi) - par1(Omega) * s1, phi2 = par2(Phi) - par2(Omega) * s2;
    const Eigen::Vector2d pT = Eigen::Vector2d(std::cos(phi1), std::sin(phi1)) / std::abs(par1(Omega)) +
                               Eigen::Vector2d(std::cos(phi2), std::sin(phi2)) / std::abs(par2(Omega));
    const double along = pT.norm() > 0 ? (x.head<2>() - r0).dot(pT.normalized()) : 0.0;
    const double gap = (x1 - x2).norm();
    points.push_back({gap, along < -kBehind, x});
    if (minGap < 0 || gap < minGap) minGap = gap;
  }
  const bool anyAhead = std::any_of(points.begin(), points.end(),
                                    [&](const Point& pt) { return !pt.behind && pt.gap <= minGap + tolerance; });
  double bestR = -1;
  for (const auto& pt : points) {
    if (pt.gap > minGap + tolerance || (anyAhead && pt.behind)) continue;
    const double r = (pt.x.head<2>() - r0).norm();
    if (bestR < 0 || r < bestR) {
      bestR = r;
      vertex = pt.x;
    }
  }
  return true;
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
