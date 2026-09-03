// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <any>
#include <vector>

#include "k4DTF/ChiSquare.h"
#include "k4DTF/ErrCode.h"

namespace k4dtf {

class DecayChain;
class FitParams;

struct FitResult {
  ErrCode status;
  double chi2{0.0};
  int ndf{0};
  int nIter{0};

  bool converged() const { return !failure(status); }
  double prob() const;

  const FitParams* fitParams() const { return m_fitParams; }
  const DecayChain* decayChain() const { return m_decayChain; }

  template <typename T>
  const T* fitContext() const {
    return std::any_cast<T>(&m_fitContext);
  }

  template <typename T>
  void setFitContext(T ctx) {
    m_fitContext = std::move(ctx);
  }

private:
  friend class DecayTreeFitter;

  const FitParams* m_fitParams{nullptr};
  const DecayChain* m_decayChain{nullptr};
  std::any m_fitContext;
};

} // namespace k4dtf
