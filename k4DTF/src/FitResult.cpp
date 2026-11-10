// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/FitResult.h"

#include "TMath.h"
#include "core/DecayChain.h"
#include "core/FitParams.h"

namespace k4dtf {

double FitResult::prob() const {
  if (ndf <= 0 || chi2 < 0) return 0.0;
  return TMath::Prob(chi2, ndf);
}

} // namespace k4dtf
