// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "core/ParticleNode.h"

namespace k4dtf {

class InternalParticleNode : public ParticleNode {
public:
  InternalParticleNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return mother() ? 8 : 7; }
  int dimM() const override { return 0; }

  bool hasEnergy() const override { return true; }
  bool hasPosition() const override { return true; }

  int posIndex() const override { return index(); }
  int lenIndex() const override { return mother() ? index() + 3 : -1; }
  int momIndex() const override { return mother() ? index() + 4 : index() + 3; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

  ErrCode projectKineConstraint(const FitParams& fitparams, Projection& p) const;

protected:
  ErrCode initMom(FitParams* fitparams) const;

  /// Seed the decay vertex from the intersection of two helix tracks made at
  /// it (directly or through resonances); leaves it unchanged if there are
  /// fewer than two.
  void seedVertex(FitParams* fitparams) const;
};

} // namespace k4dtf
