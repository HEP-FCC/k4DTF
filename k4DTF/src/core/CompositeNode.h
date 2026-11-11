// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "core/ParticleNode.h"

namespace k4DTF {

class CompositeNode : public ParticleNode {
public:
  CompositeNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return hasEnergy() ? 8 : 7; }
  int dimM() const override { return hasEnergy() ? 7 : 6; }

  int  posIndex() const override { return index(); }
  int  lenIndex() const override { return index() + 3; }
  int  momIndex() const override { return index() + 4; }
  bool hasEnergy() const override { return m_hasEnergy; }
  bool hasPosition() const override { return true; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

  double chiSquare(const FitParams& fitparams) const;

private:
  ErrCode projectRecoComposite(const FitParams& fitparams, Projection& p) const;

  Eigen::VectorXd m_measurement;
  Eigen::MatrixXd m_measurementCov;
  bool m_hasEnergy{true};
};

} // namespace k4DTF
