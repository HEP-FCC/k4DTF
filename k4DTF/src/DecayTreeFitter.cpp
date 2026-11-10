// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DecayTreeFitter.h"

#include <cmath>
#include <iostream>
#include <memory>

#include "core/DecayChain.h"
#include "core/FitParams.h"

namespace k4dtf {

namespace {

/// Covariance sub-matrix for a list of parameter indices.
Eigen::MatrixXd subCov(const FitParams& fp, const std::vector<int>& idx) { return fp.cov(idx); }

/// Fill NodeResults for `node` and its subtree in pre-order.
void extractNode(const ParticleNode* node, const FitParams& fp, int parent, bool geometry,
                 std::vector<NodeResult>& out) {
  const int me = static_cast<int>(out.size());
  out.emplace_back();
  {
    NodeResult& r = out.back();
    r.pdgCode = node->pdgCode();
    r.name = node->name();
    r.parent = parent;

    // 4-momentum
    const int mom = node->momIndex();
    if (mom >= 0) {
      const Eigen::Vector3d p(fp.par(mom), fp.par(mom + 1), fp.par(mom + 2));
      if (node->hasEnergy()) {
        r.p4 << p, fp.par(mom + 3);
        r.p4Cov = subCov(fp, {mom, mom + 1, mom + 2, mom + 3});
      } else {
        const double m = std::max(0.0, node->pdtMass());
        const double E = std::sqrt(p.squaredNorm() + m * m);
        r.p4 << p, E;
        Eigen::Matrix<double, 4, 3> J;
        J.setZero();
        J.topLeftCorner<3, 3>().setIdentity();
        if (E > 0) J.row(3) = p.transpose() / E;
        r.p4Cov = J * subCov(fp, {mom, mom + 1, mom + 2}) * J.transpose();
      }
      const double m2 = r.p4(3) * r.p4(3) - p.squaredNorm();
      r.mass = m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2);
      if (r.mass > 0) {
        Eigen::Vector4d J;
        J << -p / r.mass, r.p4(3) / r.mass;
        r.massErr = std::sqrt(std::max(0.0, double(J.transpose() * r.p4Cov * J)));
      }
    }

    // decay vertex: composites only (tracks have no position of their own)
    const int pos = node->posIndex();
    if (geometry && pos >= 0 && !node->daughters().empty()) {
      r.hasDecayVertex = true;
      r.decayVertex = Eigen::Vector3d(fp.par(pos), fp.par(pos + 1), fp.par(pos + 2));
      r.decayVertexCov = subCov(fp, {pos, pos + 1, pos + 2});
    }

    // flight length and proper decay length
    const int len = node->lenIndex();
    if (geometry && len >= 0 && mom >= 0 && node->hasEnergy()) {
      r.hasLifetime = true;
      const auto C = subCov(fp, {len, mom, mom + 1, mom + 2, mom + 3});
      const double L = fp.par(len);
      r.decayLength = L;
      r.decayLengthErr = std::sqrt(std::max(0.0, C(0, 0)));
      const Eigen::Vector3d p = r.p4.head<3>();
      const double pm = p.norm();
      if (pm > 0 && r.mass > 0) {
        const double m = r.mass, E = r.p4(3);
        r.ctau = L * m / pm;
        // d(ctau)/d(L, px, py, pz, E), with m = sqrt(E^2 - p^2)
        Eigen::Matrix<double, 5, 1> J;
        J(0) = m / pm;
        for (int i = 0; i < 3; ++i) J(1 + i) = L * (-p(i) / (m * pm) - m * p(i) / (pm * pm * pm));
        J(4) = L * E / (m * pm);
        r.ctauErr = std::sqrt(std::max(0.0, double(J.transpose() * C * J)));
      }
    }
  }
  for (const auto* dau : node->daughters()) {
    const int child = static_cast<int>(out.size());
    out[me].children.push_back(child);
    extractNode(dau, fp, me, geometry, out);
  }
}

}  // namespace

DecayTreeFitter::DecayTreeFitter(const Configuration& config)
    : m_config(config) {}

FitResult DecayTreeFitter::fit(const DecayNodeInput& input,
                                std::optional<InteractionPointInput> pv) const {
  FitResult result;
  result.status = ErrCode::success;

  auto chain = std::make_shared<DecayChain>(input, std::move(pv));

  if (chain->dim() == 0) {
    result.status = ErrCode::badsetup;
    return result;
  }

  auto fitparams = std::make_shared<FitParams>(chain->dim());

  result.status = chain->init(*fitparams);
  if (failure(result.status)) return result;

  const int maxndiverging = 3;
  int ndiverging = 0;
  double prevChiSquare = -1;
  double prevAbsDeltaChi2 = -1;

  auto finish = [&](ErrCode status, int iter) {
    result.status = status;
    result.chi2 = fitparams->chiSquare();
    result.ndf = fitparams->nDof();
    // Without any position measurement, vertices and flight lengths are
    // undetermined (they only have their prior): a kinematic fit. Do not
    // count them as fitted parameters and do not report them.
    const bool geometry = chain->hasPositionInformation();
    if (!geometry) result.ndf += chain->nGeometricParameters();
    result.nIter = iter;
    if (!failure(status)) {
      extractNode(chain->decayHead(), *fitparams, -1, geometry, result.nodes);
      result.m_fitParams = fitparams;
      result.m_decayChain = chain;
    }
    return result;
  };

  for (int iter = 0; iter < m_config.maxIter; ++iter) {
    // linearisation point of this iteration
    FitParams reference = *fitparams;

    if (iter == 0 || !fitparams->testCov()) {
      // keep the parameters, restart from the prior covariance
      ErrCode status = chain->initCov(*fitparams);
      if (failure(status)) return finish(status, iter + 1);
    } else {
      // weak prior: the inflated diagonal of the previous result
      fitparams->resetCov(1000);
    }

    ErrCode status = chain->filter(reference, *fitparams);

    double chisq = fitparams->chiSquare();
    double deltachisq = chisq - prevChiSquare;
    int ndof = fitparams->nDof();

    if (failure(status)) return finish(status, iter + 1);

    if (m_config.verbose > 0)
      std::cout << "k4dtf::DecayTreeFitter iteration " << iter << " chi2 = " << chisq << " ndf = " << ndof << "\n";

    if (iter > 0) {
      double dChisqConv =
          m_config.convergenceCriterion *
          std::max(1.0, std::max(static_cast<double>(ndof),
                                 std::min(chisq, prevChiSquare)));
      double dChisqQuit =
          std::max(2.0 * ndof, 2.0 * prevChiSquare);

      if (std::abs(deltachisq) < dChisqConv) return finish(ErrCode::success, iter + 1);

      if (iter > 1 && deltachisq > dChisqQuit) return finish(ErrCode::fastdivergingfit, iter + 1);

      if (iter > 2 &&
          std::abs(deltachisq) > prevAbsDeltaChi2 &&
          ++ndiverging >= maxndiverging) {
        return finish(ErrCode::slowdivergingfit, iter + 1);
      }
    }

    if (deltachisq < 0) ndiverging = 0;
    prevChiSquare = chisq;
    prevAbsDeltaChi2 = std::abs(deltachisq);
  }

  return finish(ErrCode::slowdivergingfit, m_config.maxIter);
}

} // namespace k4dtf
