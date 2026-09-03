// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/DecayChain.h"

#include <iostream>

#include "core/FitParams.h"

namespace k4dtf {

namespace {
int countLeaves(const ParticleNode* node) {
  if (node->daughters().empty()) return 1;
  int n = 0;
  for (auto* dau : node->daughters()) {
    n += countLeaves(dau);
  }
  return n;
}

void collectParticles(ParticleNode* node, std::vector<ParticleNode*>& out) {
  for (auto* dau : node->daughters()) {
    collectParticles(const_cast<ParticleNode*>(dau), out);
  }
  out.push_back(node);
}
} // anonymous namespace

DecayChain::DecayChain(const DecayNodeInput& input,
                       std::optional<InteractionPointInput> pv) {
  if (pv) {
    m_ipNode = std::make_unique<InteractionPointNode>(pv->position, pv->cov);
    m_mother = ParticleNode::create(input, m_ipNode.get());
    m_ipNode->addDaughter(std::move(m_mother));
    m_mother = std::move(m_ipNode);
    m_ipNode = nullptr;
  } else {
    m_mother = ParticleNode::create(input, nullptr);
  }

  m_dim = 0;
  m_mother->updateIndex(m_dim);
  m_nFinalState = countLeaves(m_mother.get());
}

std::vector<ParticleNode*> DecayChain::particles() const {
  std::vector<ParticleNode*> result;
  collectParticles(m_mother.get(), result);
  return result;
}

ErrCode DecayChain::init(FitParams& fitparams) {
  if (fitparams.dim() != m_dim) {
    fitparams = FitParams(m_dim);
  }
  fitparams.resetPar();
  ErrCode status = m_mother->initPar1(&fitparams);
  status |= m_mother->initPar2(&fitparams);
  fitparams.resetCov();
  status |= m_mother->initCov(&fitparams);
  return status;
}

ErrCode DecayChain::filter(FitParams& reference, FitParams& fitparams) const {
  return m_mother->filter(reference, fitparams);
}

} // namespace k4dtf
