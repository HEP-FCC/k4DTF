// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "k4DTF/core/ParticleMeasurement.h"

namespace k4dtf {

/// One node of the decay tree given to DecayTreeFitter::fit().
///
/// Final-state particles carry a measurement; composite nodes carry daughters
/// (their measurement is ignored).
struct DecayNodeInput {
  int pdgCode{0};
  std::string name;
  ParticleMeasurement measurement;
  std::vector<DecayNodeInput> daughters;

  /// Mass hypothesis [GeV]. Required for tracks (their energy is computed from
  /// it); used as the constraint value of composites if massConstraintValue <= 0.
  double mass{-1};

  bool applyMassConstraint{false};
  double massConstraintValue{0};

  /// Composite nodes only. A resonance (strong or electromagnetic decay: D*,
  /// K*, J/psi, ...) decays where it is produced: it has no vertex or flight
  /// length of its own. Otherwise the node gets a decay vertex and a flight
  /// length from its production vertex.
  bool isResonance{false};
};

struct InteractionPointInput {
  Eigen::Vector3d position;
  Eigen::Matrix3d cov;
};

} // namespace k4dtf
