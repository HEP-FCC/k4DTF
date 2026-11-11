// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/ResonanceNode.h"

namespace k4DTF {

ResonanceNode::ResonanceNode(const DecayNodeInput& input, const ParticleNode* mother)
    : InternalParticleNode(input, mother) {}

ErrCode ResonanceNode::initPar1(FitParams* fitparams) {
  ErrCode status = ErrCode::success;
  for (auto* dau : daughters()) {
    status |= dau->initPar1(fitparams);
  }
  return status;
}

ErrCode ResonanceNode::initPar2(FitParams* fitparams) {
  ErrCode status = ErrCode::success;
  for (auto* dau : daughters()) {
    status |= dau->initPar2(fitparams);
  }
  initMom(fitparams);
  return status;
}

ErrCode ResonanceNode::filter(const FitParams& reference, FitParams& params) const {
  for (auto* dau : daughters()) {
    ErrCode status = dau->filter(reference, params);
    if (failure(status)) return status;
  }

  ErrCode status = applyConstraint(reference, params, 4,
                                   &InternalParticleNode::projectKineConstraint);
  if (hasMassConstraint()) {
    status |= applyConstraint(reference, params, 1,
                              &ParticleNode::projectMassConstraint);
  }
  return status;
}

} // namespace k4DTF
