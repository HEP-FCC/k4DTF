// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

namespace k4dtf {

/// Helix parameters of a charged track, in the convention of
/// edm4hep::TrackState as written by k4SimDelphes (Delphes TrkUtil):
///
///   par = (D0 [mm], phi [rad], omega [1/mm], Z0 [mm], tanLambda)
///
///   - phi       azimuth of the momentum at the point of closest approach (PCA)
///               to the reference point
///   - omega     signed curvature, sign(omega) = sign(charge * Bz)
///   - PCA       = reference + (-D0 sin(phi), D0 cos(phi), Z0)
///   - pT [GeV]  = 0.299792458e-3 * |Bz [T]| / |omega [1/mm]|
///
/// Units: mm, GeV, T.
namespace helix {

using Vector5d = Eigen::Matrix<double, 5, 1>;
using Matrix5d = Eigen::Matrix<double, 5, 5>;

enum Index { D0 = 0, Phi = 1, Omega = 2, Z0 = 3, TanLambda = 4 };

/// pT [GeV] = kC * Bz [T] * R [mm]
constexpr double kC = 0.299792458e-3;

/// Wrap an angle into (-pi, pi].
double wrapPhi(double phi);

/// Helix parameters of a particle with charge q [e] passing through x [mm]
/// with momentum p [GeV] at that point, with respect to the reference point.
/// Exact (no small-displacement approximation).
Vector5d fromPositionMomentum(const Eigen::Vector3d& x, const Eigen::Vector3d& p, double q, double Bz,
                              const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Momentum [GeV] at the PCA.
Eigen::Vector3d momentum(const Vector5d& par, double Bz);

/// Point of closest approach [mm].
Eigen::Vector3d pca(const Vector5d& par, const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Charge [e] from the sign of omega and Bz (0 if either is 0).
double charge(const Vector5d& par, double Bz);

/// Jacobian d(px, py, pz) / d(par) at the PCA (only phi, omega and tanLambda enter).
Eigen::Matrix<double, 3, 5> momentumJacobian(const Vector5d& par, double Bz);

/// Covariance of (px, py, pz) at the PCA from the helix covariance.
Eigen::Matrix3d momentumCovariance(const Vector5d& par, const Matrix5d& cov, double Bz);

}  // namespace helix
}  // namespace k4dtf
