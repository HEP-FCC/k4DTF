// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/FitParams.h"

#include <iomanip>
#include <iostream>

namespace k4dtf {

FitParams::FitParams(int dim)
    : m_dim(dim),
      m_par(dim),
      m_cov(dim, dim),
      m_nConstraintsVec(dim, 0) {
  m_par.setZero();
  m_cov.setZero();
}

void FitParams::resetPar() { m_par.setZero(); }

void FitParams::resetCov(double scale) {
  for (int row = 0; row < m_dim; ++row) {
    for (int col = 0; col < row; ++col) m_cov(row, col) = 0;
    m_cov(row, row) *= scale;
    if (m_cov(row, row) < 0) m_cov(row, row) = -m_cov(row, row);
  }
  m_chiSquare = 0;
  m_nConstraints = 0;
  std::fill(m_nConstraintsVec.begin(), m_nConstraintsVec.end(), 0);
  m_chiSquareMap.clear();
}

bool FitParams::testCov() const {
  for (int row = 0; row < m_dim; ++row) {
    if (m_cov(row, row) <= 0) return false;
  }
  return true;
}

void FitParams::print() const {
  std::cout << std::setw(3) << "idx" << std::setw(15) << "val" << std::setw(15) << "err" << std::endl;
  std::cout << std::setprecision(5);
  for (int row = 0; row < m_dim; ++row) {
    std::cout << std::setw(3) << row << std::setw(15) << m_par(row) << std::setw(15)
              << std::sqrt(m_cov(row, row)) << std::endl;
  }
  std::cout << std::setprecision(10);
}

Eigen::MatrixXd FitParams::cov(const std::vector<int>& indexVec) const {
  int nrow = indexVec.size();
  Eigen::MatrixXd thecov(nrow, nrow);
  for (int row = 0; row < nrow; ++row) {
    for (int col = 0; col <= row; ++col) {
      thecov(row, col) = m_cov(indexVec[row], indexVec[col]);
      thecov(col, row) = thecov(row, col);
    }
  }
  return thecov;
}

Eigen::VectorXd FitParams::par(const std::vector<int>& indexVec) const {
  int nrow = indexVec.size();
  Eigen::VectorXd thepar(nrow);
  for (int row = 0; row < nrow; ++row) {
    thepar(row) = m_par(indexVec[row]);
  }
  return thepar;
}

void FitParams::addChiSquare(double chisq, int nconstraints, const ParticleNode* p) {
  m_chiSquare += chisq;
  m_nConstraints += nconstraints;
  if (p) m_chiSquareMap[p] += ChiSquare(chisq, nconstraints);
}

ChiSquare FitParams::chiSquare(const ParticleNode& p) const {
  auto it = m_chiSquareMap.find(&p);
  return it != m_chiSquareMap.end() ? it->second : ChiSquare();
}

} // namespace k4dtf
