// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

namespace k4DTF {

class Projection {
public:
  Projection(int dimP, int dimC)
      : m_H(dimC, dimP), m_r(dimC), m_V(dimC, dimC), m_offset(0) {
    m_H.setZero();
    m_r.setZero();
    m_V.setZero();
  }

  const Eigen::MatrixXd& H() const { return m_H; }
  Eigen::MatrixXd&       H() { return m_H; }

  double& H(int row, int col) { return m_H(m_offset + row, col); }

  const Eigen::VectorXd& r() const { return m_r; }
  Eigen::VectorXd&       r() { return m_r; }
  double& r(int row) { return m_r(m_offset + row); }

  const Eigen::MatrixXd& V() const { return m_V; }
  Eigen::MatrixXd&       V() { return m_V; }
  double& V(int row, int col) { return m_V(m_offset + row, m_offset + col); }

  void reset() {
    m_H.setZero();
    m_r.setZero();
    m_V.setZero();
    m_offset = 0;
  }

  double chiSquare() const {
    Eigen::MatrixXd W = m_V;
    Eigen::VectorXd result = W.inverse() * m_r;
    return m_r.transpose() * result;
  }

  void incrementOffset(unsigned int i) { m_offset += i; }
  unsigned int offset() const { return m_offset; }

private:
  Eigen::MatrixXd m_H;   ///< Jacobian (dimC x dimP)
  Eigen::VectorXd m_r;   ///< constraint residual (dimC)
  Eigen::MatrixXd m_V;   ///< constraint covariance (dimC x dimC)
  unsigned int    m_offset{0};
};

} // namespace k4DTF
