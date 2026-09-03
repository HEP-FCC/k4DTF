// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/InternalParticleNode.h"

namespace k4dtf {

class ResonanceNode : public InternalParticleNode {
public:
  ResonanceNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return 4; }
  int dimM() const override { return 0; }

  int  posIndex() const override { return mother()->posIndex(); }
  int  momIndex() const override { return index(); }
  int  lenIndex() const override { return -1; }
  bool hasPosition() const override { return false; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;
};

} // namespace k4dtf
