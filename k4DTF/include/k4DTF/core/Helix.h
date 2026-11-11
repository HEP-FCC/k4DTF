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

/// Jacobian d(par) / d(x, y, z, px, py, pz) of fromPositionMomentum().
Eigen::Matrix<double, 5, 6> jacobian(const Eigen::Vector3d& x, const Eigen::Vector3d& p, double q, double Bz,
                                     const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Momentum [GeV] at the PCA.
Eigen::Vector3d momentum(const Vector5d& par, double Bz);

/// Point of closest approach [mm].
Eigen::Vector3d pca(const Vector5d& par, const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Point on the helix at signed transverse arc length s [mm] from the PCA
/// (s > 0 along the momentum). The helix turns by -omega s.
Eigen::Vector3d positionAt(const Vector5d& par, double s, const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Momentum [GeV] at arc length s.
Eigen::Vector3d momentumAt(const Vector5d& par, double Bz, double s);

/// Arc length of the point of the helix closest to `point` in the transverse
/// plane (within half a turn of the PCA).
double arcLengthNear(const Vector5d& par, const Eigen::Vector3d& point,
                     const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

/// Approximate common point of two helices, to seed a vertex fit. Candidates
/// in the transverse plane: the intersections of the two circles (or, if
/// they do not meet, their closest points), and the points closest to the
/// reference point. Of the candidates where the helices are within
/// `tolerance` [mm] of their smallest 3D distance, the one nearest the beam
/// line (reference point) is returned, as the midpoint between the helices.
/// Candidates more than 1 mm behind the reference point, along the pair's
/// transverse momentum, are only taken if there is no other.
/// (Nearly parallel tracks stay close over a long stretch, and their circles
/// also cross far away.) A tolerance of a few times the impact-parameter
/// resolution separates real ambiguities from boosted displaced decays.
/// Returns false for (near) concentric circles.
bool intersect(const Vector5d& par1, const Eigen::Vector3d& ref1, const Vector5d& par2,
               const Eigen::Vector3d& ref2, Eigen::Vector3d& vertex, double tolerance = 0.1);

/// Charge [e] from the sign of omega and Bz (0 if either is 0).
double charge(const Vector5d& par, double Bz);

/// Jacobian d(px, py, pz) / d(par) at the PCA (only phi, omega and tanLambda enter).
Eigen::Matrix<double, 3, 5> momentumJacobian(const Vector5d& par, double Bz);

/// Covariance of (px, py, pz) at the PCA from the helix covariance.
Eigen::Matrix3d momentumCovariance(const Vector5d& par, const Matrix5d& cov, double Bz);

}  // namespace helix
}  // namespace k4dtf
