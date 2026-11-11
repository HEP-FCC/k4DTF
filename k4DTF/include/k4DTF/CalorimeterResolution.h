// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "k4DTF/core/ParticleMeasurement.h"

namespace k4DTF {

/// Photon measurement resolution of a calorimeter with projective towers, for
/// clusters that come without errors (or without a position).
///
/// The defaults describe the Delphes IDEA card (FCC-config card_IDEA.tcl):
/// crystal ECAL with sigma(E) = sqrt((0.005 E)^2 + 0.03^2 E + 0.002^2) [GeV],
/// towers of 0.02 x 0.02 in (eta, phi), front face at R = 2250 mm and
/// |z| = 2500 mm (the ParticlePropagator cylinder).
///
/// Delphes puts a photon at the (eta, phi) of its tower; with SmearTowerCenter
/// (the IDEA default) the point is drawn uniformly within the tower, otherwise
/// it is the tower centre. Relative to the true impact point (itself anywhere
/// in the tower) the error is then the difference of two uniform
/// distributions, size / sqrt(6), or one, size / sqrt(12).
/// k4SimDelphes v00-09 writes size / sqrt(12) in both cases.
struct CalorimeterResolution {
  double towerSize{0.02};       ///< tower size in eta and in phi
  bool smearTowerCenter{true};  ///< Delphes SmearTowerCenter
  double stochastic{0.03};      ///< [GeV^1/2]
  double noise{0.002};          ///< [GeV]
  double constant{0.005};
  double radius{2250};      ///< [mm] barrel front face
  double halfLength{2500};  ///< [mm] endcap front face

  /// Angular resolution in eta and in phi.
  double angularResolution() const;

  /// Energy variance [GeV^2] at energy E [GeV].
  double energyVariance(double E) const;

  /// Covariance [mm^2] of a cluster position from the angular resolution,
  /// sigma(theta) = sigma(eta) sin(theta): transverse to the line of sight
  /// from the origin, as in k4SimDelphes. (The photon constraint does not
  /// depend on the position along its direction.)
  Eigen::Matrix3d positionCovariance(const Eigen::Vector3d& position) const;

  /// Point where the line from the origin along `direction` reaches the
  /// front face (barrel or endcap).
  Eigen::Vector3d frontFace(const Eigen::Vector3d& direction) const;

  /// Photon measurement for a cluster at `position` with energy E, with
  /// the covariance from this resolution model.
  PhotonMeasurement measurement(const Eigen::Vector3d& position, double E) const;
};

}  // namespace k4DTF
