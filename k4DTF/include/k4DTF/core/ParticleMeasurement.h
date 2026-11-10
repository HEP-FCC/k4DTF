// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <variant>
#include <vector>

#include <Eigen/Dense>

namespace k4dtf {

/// Helix parameters with covariance, see k4DTF/core/Helix.h for the conventions.
struct HelixMeasurement {
  Eigen::Matrix<double, 5, 1> par;
  Eigen::Matrix<double, 5, 5> cov;
  Eigen::Vector3d referencePoint{Eigen::Vector3d::Zero()};  ///< [mm]
  double Bz{0};                                             ///< [T]
};

/// A charged track.
///
/// With a helix the track constrains both the position of its production
/// vertex and its momentum (5 constraints); this is what a decay-tree fit with
/// vertices needs. Without a helix only the 3-momentum is used, and no vertex
/// position can be determined from this track.
struct TrackMeasurement {
  Eigen::Vector3d momentum;  ///< used if there is no helix
  Eigen::Matrix3d cov;       ///< used if there is no helix
  double charge{0};
  std::optional<HelixMeasurement> helix;
};

struct PhotonMeasurement {
  Eigen::Vector3d clusterPos;
  double energy;
  Eigen::Matrix3d cov;
};

struct CompositeMeasurement {
  Eigen::Vector3d vertex;
  Eigen::Vector4d momentum;
  Eigen::Matrix<double, 7, 7> cov;
  bool hasEnergy{true};
};

struct ResonanceMeasurement {
  Eigen::Vector4d momentum;
  Eigen::Matrix4d cov;
};

struct MergedPi0Measurement {
  struct Cluster {
    Eigen::Vector3d pos;
    double energy;
  };
  std::vector<Cluster> clusters;
  Eigen::Matrix4d cov;
};

struct InteractionPointMeasurement {
  Eigen::Vector3d position;
  Eigen::Matrix3d cov;
};

using ParticleMeasurement = std::variant<TrackMeasurement, PhotonMeasurement,
                                         CompositeMeasurement,
                                         ResonanceMeasurement,
                                         MergedPi0Measurement>;

} // namespace k4dtf
