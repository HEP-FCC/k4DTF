// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include "k4DTF/DTFParticle.h"

namespace k4DTF {

/// Output of a single DecayTreeFitter::fit() call.
struct FitResult {
  enum class Status { Success, MaxIterationsReached, SingularMatrix, Failed };

  Status status{Status::Failed};
  double chi2{0.0};
  int ndf{0};
  int nIter{0};

  /// Fitted particles in the same order as the input (daughters first, then
  /// composite nodes up to the head of the decay tree).
  std::vector<DTFParticle> fittedParticles;

  bool converged() const { return status == Status::Success; }

  /// Returns chi2 probability (TMath::Prob(chi2, ndf)); 0 if ndf <= 0.
  /// TODO: implement using TMath::Prob(chi2, ndf).
  double prob() const { return 0.0; }
};

}  // namespace k4DTF