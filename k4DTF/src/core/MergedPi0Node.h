// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "core/ParticleNode.h"

namespace k4DTF {

class MergedPi0Node : public ParticleNode {
public:
  MergedPi0Node(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return hasMassConstraint() ? 3 : 4; }
  int dimM() const override { return dim(); }

  int  momIndex() const override { return index(); }
  bool hasEnergy() const override { return !hasMassConstraint(); }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;
  ErrCode initCov(FitParams* fitparams) const override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

private:
  ErrCode projectPi0Constraint(const FitParams& fitparams, Projection& p) const;

  Eigen::Vector4d m_p4;
  Eigen::Matrix4d m_p4Cov;
  std::vector<Eigen::Vector3d> m_clusterPositions;
  std::vector<double> m_clusterEnergies;
};

} // namespace k4DTF
