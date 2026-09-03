// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "k4DTF/ErrCode.h"
#include "core/FitParams.h"

namespace k4dtf {

class KalmanCalculator {
public:
  ErrCode init(const Eigen::VectorXd& r, const Eigen::MatrixXd& H,
               const FitParams& fitparams,
               const Eigen::MatrixXd* V = nullptr, int weight = 1);

  void   updatePar(FitParams& fitparams);
  void   updatePar(const Eigen::VectorXd& prediction, FitParams& fitparams);
  void   updateCov(FitParams& fitparams);
  double chisq() const { return m_chisq; }

private:
  int              m_nconstraints{0};
  int              m_nparameters{0};
  Eigen::VectorXd  m_r;
  Eigen::MatrixXd  m_H;
  Eigen::MatrixXd  m_matrixCGT;
  Eigen::MatrixXd  m_matrixR;
  Eigen::MatrixXd  m_matrixRinv;
  Eigen::MatrixXd  m_matrixK;
  double           m_chisq{0};
  int              m_ierr{0};
};

} // namespace k4dtf
