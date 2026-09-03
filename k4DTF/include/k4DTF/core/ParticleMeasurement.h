// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <variant>
#include <vector>

#include <Eigen/Dense>

namespace k4dtf {

struct TrackMeasurement {
  Eigen::Vector3d momentum;
  Eigen::Matrix3d cov;
  double charge{0};
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
