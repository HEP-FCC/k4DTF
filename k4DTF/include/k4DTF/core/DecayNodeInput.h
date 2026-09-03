// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "k4DTF/core/ParticleMeasurement.h"

namespace k4dtf {

struct DecayNodeInput {
  int pdgCode{0};
  std::string name;
  ParticleMeasurement measurement;
  std::vector<DecayNodeInput> daughters;
  bool applyMassConstraint{false};
  double massConstraintValue{0};
};

struct InteractionPointInput {
  Eigen::Vector3d position;
  Eigen::Matrix3d cov;
};

} // namespace k4dtf
