// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "core/ParticleNode.h"

namespace k4dtf {

class InteractionPointNode : public ParticleNode {
public:
  InteractionPointNode(const Eigen::Vector3d& pos, const Eigen::Matrix3d& cov);

  int dim() const override { return 3; }
  int dimM() const override { return 3; }

  int posIndex() const override { return index(); }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;
  ErrCode initCov(FitParams* fitparams) const override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

private:
  ErrCode projectIPConstraint(const FitParams& fitparams, Projection& p) const;

  Eigen::Vector3d m_ipPos;
  Eigen::Matrix3d m_ipCov;
};

} // namespace k4dtf
