// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/KalmanCalculator.h"

namespace k4dtf {

ErrCode KalmanCalculator::init(const Eigen::VectorXd& r, const Eigen::MatrixXd& H,
                               const FitParams& fitparams,
                               const Eigen::MatrixXd* V, int weight) {
  ErrCode status = ErrCode::success;
  m_nconstraints = r.size();
  m_nparameters  = fitparams.dim();
  m_r = r;
  m_H = H;

  const Eigen::MatrixXd& C = fitparams.cov();

  // C * H^T
  m_matrixCGT = C * H.transpose();

  // R = H * C * H^T + V
  if (V) {
    m_matrixR = H * C * H.transpose() + weight * (*V);
  } else {
    m_matrixR = H * C * H.transpose();
  }

  // Invert R
  m_matrixRinv = m_matrixR.inverse();

  // Check for singularity (Eigen doesn't set an error flag, but we can check)
  // A simple heuristic: if the determinant is extremely small, flag it
  if (std::abs(m_matrixR.determinant()) < 1e-30) {
    status |= ErrCode::inversionerror;
  }

  // Kalman gain: K = C * H^T * R^{-1}
  m_matrixK = m_matrixCGT * m_matrixRinv;

  m_chisq = -1;
  return status;
}

void KalmanCalculator::updatePar(FitParams& fitparams) {
  // par -= K * r
  fitparams.par() -= m_matrixK * m_r;
  // chi2 = r^T * R^{-1} * r
  m_chisq = m_r.transpose() * m_matrixRinv * m_r;
}

void KalmanCalculator::updatePar(const Eigen::VectorXd& prediction, FitParams& fitparams) {
  Eigen::VectorXd rprime = m_r + m_H * (prediction - fitparams.par());
  fitparams.par() = prediction - m_matrixK * rprime;
  m_chisq = rprime.transpose() * m_matrixRinv * rprime;
}

void KalmanCalculator::updateCov(FitParams& fitparams) {
  // deltaCov = -C * H^T * K^T = -K * R * K^T  (Joseph form for stability)
  // Using the simpler form: C -= K * R * K^T
  fitparams.cov() -= m_matrixK * m_matrixR * m_matrixK.transpose();
}

} // namespace k4dtf
