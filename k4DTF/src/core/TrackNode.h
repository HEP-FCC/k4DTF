// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <Eigen/Dense>

#include "core/ParticleNode.h"

namespace k4dtf {

/// Charged final-state particle. Parameters: its momentum (px, py, pz) at the
/// production vertex; the vertex position belongs to the mother.
///
/// With a helix measurement the constraint is h(x_mother, p) = h_measured
/// (5 constraints), which ties the track to its mother's decay vertex; without
/// one, the momentum is constrained directly (3 constraints).
class TrackNode : public ParticleNode {
public:
  TrackNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return 3; }
  int dimM() const override { return m_helix ? 5 : 3; }

  int  momIndex() const override { return index(); }
  bool hasEnergy() const override { return false; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;
  ErrCode initCov(FitParams* fitparams) const override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

  bool hasHelix() const { return m_helix.has_value(); }
  const HelixMeasurement* helixMeasurement() const { return m_helix ? &*m_helix : nullptr; }

private:
  ErrCode projectTrackConstraint(const FitParams& fitparams, Projection& p) const;
  ErrCode projectHelixConstraint(const FitParams& fitparams, Projection& p) const;

  Eigen::Vector3d m_momentum;
  Eigen::Matrix3d m_momCov;
  std::optional<HelixMeasurement> m_helix;
};

} // namespace k4dtf
