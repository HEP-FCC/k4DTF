// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <optional>
#include <vector>

#include "core/DecayNodeInput.h"
#include "core/InteractionPointNode.h"
#include "core/ParticleNode.h"
#include "k4DTF/ErrCode.h"

namespace k4dtf {

class FitParams;

class DecayChain {
public:
  DecayChain(const DecayNodeInput& input,
             std::optional<InteractionPointInput> pv = std::nullopt);

  int dim() const { return m_dim; }
  int nFinalState() const { return m_nFinalState; }

  const ParticleNode* head() const { return m_mother.get(); }
  ParticleNode*       head()       { return m_mother.get(); }

  std::vector<ParticleNode*> particles() const;

  ErrCode init(FitParams& fitparams);
  ErrCode filter(FitParams& reference, FitParams& fitparams) const;

private:
  std::unique_ptr<ParticleNode> m_mother;
  std::unique_ptr<InteractionPointNode> m_ipNode;
  int m_dim{0};
  int m_nFinalState{0};
};

} // namespace k4dtf
