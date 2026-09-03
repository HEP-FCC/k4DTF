// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ParticleNode.h"

namespace k4dtf {

class PhotonNode : public ParticleNode {
public:
  PhotonNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return 3; }
  int dimM() const override { return 3; }

  int  momIndex() const override { return index(); }
  bool hasEnergy() const override { return false; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;
  ErrCode initCov(FitParams* fitparams) const override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

private:
  ErrCode projectPhotonConstraint(const FitParams& fitparams, Projection& p) const;

  Eigen::Vector3d m_clusterPos;
  double m_clusterEnergy{0};
  Eigen::Matrix3d m_clusterCov;
};

} // namespace k4dtf
