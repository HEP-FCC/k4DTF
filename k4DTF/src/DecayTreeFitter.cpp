// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DecayTreeFitter.h"

#include <cmath>
#include <iostream>

#include "core/DecayChain.h"
#include "core/FitParams.h"

namespace k4dtf {

DecayTreeFitter::DecayTreeFitter(const Configuration& config)
    : m_config(config) {}

FitResult DecayTreeFitter::fit(const DecayNodeInput& input,
                                std::optional<InteractionPointInput> pv) {
  FitResult result;
  result.status = ErrCode::success;

  DecayChain chain(input, std::move(pv));

  if (chain.dim() == 0) {
    result.status = ErrCode::badsetup;
    return result;
  }

  FitParams fitparams(chain.dim());

  result.status = chain.init(fitparams);
  if (failure(result.status)) return result;

  const int maxndiverging = 3;
  int ndiverging = 0;
  double prevChiSquare = -1;
  double prevAbsDeltaChi2 = -1;

  for (int iter = 0; iter < m_config.maxIter; ++iter) {
    bool firstpass = (iter == 0);

    FitParams reference = fitparams;

    if (firstpass || !fitparams.testCov()) {
      fitparams.resetCov(1000);
      chain.init(fitparams);
    } else {
      fitparams.resetCov(1000);
    }

    ErrCode status = chain.filter(reference, fitparams);

    double chisq = fitparams.chiSquare();
    double deltachisq = chisq - prevChiSquare;
    int ndof = fitparams.nDof();

    if (failure(status)) {
      result.status = status;
      result.chi2 = chisq;
      result.ndf = ndof;
      result.nIter = iter + 1;
      return result;
    }

    if (iter > 0) {
      double dChisqConv =
          m_config.convergenceCriterion *
          std::max(1.0, std::max(static_cast<double>(ndof),
                                 std::min(chisq, prevChiSquare)));
      double dChisqQuit =
          std::max(2.0 * ndof, 2.0 * prevChiSquare);

      if (std::abs(deltachisq) < dChisqConv) {
        result.chi2 = chisq;
        result.ndf = ndof;
        result.nIter = iter + 1;
        result.m_fitParams = &fitparams;
        result.m_decayChain = &chain;
        return result;
      }

      if (iter > 1 && deltachisq > dChisqQuit) {
        result.status = ErrCode::fastdivergingfit;
        result.chi2 = chisq;
        result.ndf = ndof;
        result.nIter = iter + 1;
        return result;
      }

      if (iter > 2 &&
          std::abs(deltachisq) > prevAbsDeltaChi2 &&
          ++ndiverging >= maxndiverging) {
        result.status = ErrCode::slowdivergingfit;
        result.chi2 = chisq;
        result.ndf = ndof;
        result.nIter = iter + 1;
        return result;
      }
    }

    if (deltachisq < 0) ndiverging = 0;
    prevChiSquare = chisq;
    prevAbsDeltaChi2 = std::abs(deltachisq);
  }

  result.status = ErrCode::slowdivergingfit;
  result.chi2 = fitparams.chiSquare();
  result.ndf = fitparams.nDof();
  result.nIter = m_config.maxIter;
  return result;
}

} // namespace k4dtf
