// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <any>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "k4DTF/ChiSquare.h"
#include "k4DTF/ErrCode.h"

namespace k4dtf {

class DecayChain;
class FitParams;

/// Fitted state of one node of the decay tree.
///
/// FitResult::nodes lists the nodes of the DecayNodeInput tree in depth-first
/// pre-order: node 0 is the head, every composite is followed by its subtree,
/// daughters in input order. The interaction point is not a node.
struct NodeResult {
  int pdgCode{0};
  std::string name;
  int parent{-1};             ///< index of the mother node, -1 for the head
  std::vector<int> children;  ///< indices of the daughter nodes

  /// Fitted (px, py, pz, E) [GeV] and covariance. For tracks E is computed
  /// from the mass hypothesis.
  Eigen::Vector4d p4{Eigen::Vector4d::Zero()};
  Eigen::Matrix4d p4Cov{Eigen::Matrix4d::Zero()};
  double mass{0};
  double massErr{0};

  /// Fitted decay vertex [mm] of composite nodes. For a resonance this is
  /// the vertex it is produced at (and decays at).
  bool hasDecayVertex{false};
  Eigen::Vector3d decayVertex{Eigen::Vector3d::Zero()};
  Eigen::Matrix3d decayVertexCov{Eigen::Matrix3d::Zero()};

  /// Flight length [mm] from the production vertex (mother's decay vertex or
  /// interaction point) and proper decay length c*tau = L m / p [mm]. Only for
  /// composites with a lifetime and a production vertex.
  bool hasLifetime{false};
  double decayLength{0};
  double decayLengthErr{0};
  double ctau{0};
  double ctauErr{0};
};

struct FitResult {
  ErrCode status{ErrCode::success};
  double chi2{0.0};
  int ndf{0};
  int nIter{0};

  /// Per-node results (see NodeResult); filled when the fit converged.
  std::vector<NodeResult> nodes;

  bool converged() const { return !failure(status); }

  /// chi2 probability; 0 if ndf <= 0.
  double prob() const;

  /// Full parameter vector and covariance, and the internal tree; valid for as
  /// long as the FitResult (or a copy of it) exists. Null if the fit failed.
  const FitParams* fitParams() const { return m_fitParams.get(); }
  const DecayChain* decayChain() const { return m_decayChain.get(); }

  template <typename T>
  const T* fitContext() const {
    return std::any_cast<T>(&m_fitContext);
  }

  template <typename T>
  void setFitContext(T ctx) {
    m_fitContext = std::move(ctx);
  }

private:
  friend class DecayTreeFitter;

  std::shared_ptr<const FitParams> m_fitParams;
  std::shared_ptr<const DecayChain> m_decayChain;
  std::any m_fitContext;
};

} // namespace k4dtf
