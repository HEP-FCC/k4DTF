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

namespace k4DTF {

class FitParams;

class DecayChain {
public:
  DecayChain(const DecayNodeInput& input,
             std::optional<InteractionPointInput> pv = std::nullopt);

  int dim() const { return m_dim; }
  int nFinalState() const { return m_nFinalState; }

  /// Top of the internal tree: the interaction point if one was given,
  /// otherwise the head of the decay.
  const ParticleNode* head() const { return m_mother.get(); }
  ParticleNode*       head()       { return m_mother.get(); }

  /// Head of the decay (the node built from the DecayNodeInput root).
  const ParticleNode* decayHead() const { return m_decayHead; }

  /// The interaction point node, or nullptr.
  const ParticleNode* interactionPoint() const { return m_hasIP ? m_mother.get() : nullptr; }

  std::vector<ParticleNode*> particles() const;

  /// False if no measurement in the tree constrains a position (momentum-only
  /// tracks, no interaction point): then vertices and flight lengths are not
  /// determined by the fit.
  bool hasPositionInformation() const { return m_hasPositionInfo; }

  /// Number of vertex and flight-length parameters; these are free when there
  /// is no position information.
  int nGeometricParameters() const { return m_nGeometric; }

  ErrCode init(FitParams& fitparams);
  ErrCode initCov(FitParams& fitparams) const;
  ErrCode filter(FitParams& reference, FitParams& fitparams) const;

private:
  std::unique_ptr<ParticleNode> m_mother;
  std::unique_ptr<InteractionPointNode> m_ipNode;
  const ParticleNode* m_decayHead{nullptr};
  bool m_hasIP{false};
  bool m_hasPositionInfo{false};
  int m_nGeometric{0};
  int m_dim{0};
  int m_nFinalState{0};
};

} // namespace k4DTF
