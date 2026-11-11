// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>
#include <map>
#include <vector>

#include "k4DTF/ChiSquare.h"

namespace k4DTF {

class ParticleNode;

class FitParams {
public:
  explicit FitParams(int dim);

  int    dim() const { return m_dim; }
  double chiSquare() const { return m_chiSquare; }
  int    nConstraints() const { return m_nConstraints; }
  int    nDof() const { return m_nConstraints - m_dim; }

  double par(int i) const { return m_par(i); }
  void   setPar(int i, double val) { m_par(i) = val; }

  Eigen::VectorXd&       par() { return m_par; }
  const Eigen::VectorXd& par() const { return m_par; }

  Eigen::MatrixXd&       cov() { return m_cov; }
  const Eigen::MatrixXd& cov() const { return m_cov; }

  double cov(int i, int j) const { return m_cov(i, j); }

  double err(int i) const { return std::sqrt(m_cov(i, i)); }

  Eigen::VectorXd par(const std::vector<int>& indexVec) const;
  Eigen::MatrixXd cov(const std::vector<int>& indexVec) const;

  void resetPar();
  void resetCov(double scale = 100);
  void print() const;
  bool testCov() const;

  void addChiSquare(double chisq, int nconstraints, const ParticleNode* p);
  ChiSquare chiSquare(const ParticleNode& p) const;

private:
  int                                    m_dim;
  Eigen::VectorXd                        m_par;
  Eigen::MatrixXd                        m_cov;
  double                                 m_chiSquare{0};
  int                                    m_nConstraints{0};
  std::vector<int>                       m_nConstraintsVec;
  std::map<const ParticleNode*, ChiSquare> m_chiSquareMap;
};

} // namespace k4DTF
