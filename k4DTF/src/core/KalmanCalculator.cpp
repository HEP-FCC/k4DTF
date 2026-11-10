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
    m_V = weight * (*V);
  } else {
    m_V = Eigen::MatrixXd::Zero(m_nconstraints, m_nconstraints);
  }
  m_matrixR = H * C * H.transpose() + m_V;

  // Invert R. R is a covariance matrix, so check positive definiteness with a
  // Cholesky decomposition rather than with the determinant, which scales with
  // the units (a 5x5 helix covariance in mm and rad has det ~ 1e-40).
  Eigen::LLT<Eigen::MatrixXd> llt(m_matrixR);
  if (llt.info() != Eigen::Success) {
    status |= ErrCode::inversionerror;
    return status;
  }
  m_matrixRinv = llt.solve(Eigen::MatrixXd::Identity(m_nconstraints, m_nconstraints));

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
  // Joseph form, C' = (1 - K H) C (1 - K H)^T + K V K^T. Algebraically equal
  // to C - K R K^T, but stays positive semi-definite in floating point: the
  // priors (e.g. 1 m on flight lengths) and the measurements (1e-7 /mm on
  // omega) differ by many orders of magnitude, and the short form loses
  // positive definiteness after a few constraints.
  const Eigen::MatrixXd A =
      Eigen::MatrixXd::Identity(m_nparameters, m_nparameters) - m_matrixK * m_H;
  Eigen::MatrixXd C = A * fitparams.cov() * A.transpose() + m_matrixK * m_V * m_matrixK.transpose();
  fitparams.cov() = 0.5 * (C + C.transpose());
}

} // namespace k4dtf
