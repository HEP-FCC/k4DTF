// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>

#include "core/ParticleNode.h"

namespace k4DTF {

/// Photon from a calorimeter cluster. Parameters: its momentum (px, py, pz);
/// it is produced at its mother's vertex x and points at the cluster c.
///
/// Constraint (3 equations), in the form of the BaBar DecayTreeFitter
/// (RecoPhoton), which works in any direction: with i1 the coordinate in
/// which the photon moves most (largest |p_i|) and i2, i3 the other two,
///
///   r0 = (c - x)_i2 - (c - x)_i1 p_i2 / p_i1
///   r1 = (c - x)_i3 - (c - x)_i1 p_i3 / p_i1
///   r2 = |p| - E
///
/// The first two say that the cluster lies on the line from x along p. They
/// do not divide by a small number for any direction, unlike a slope
/// parametrisation tx = px/pz, which diverges for barrel photons.
class PhotonNode : public ParticleNode {
public:
  PhotonNode(const DecayNodeInput& input, const ParticleNode* mother);

  int dim() const override { return 3; }
  int dimM() const override { return 3; }

  int  momIndex() const override { return index(); }
  bool hasEnergy() const override { return false; }

  ErrCode initPar1(FitParams* fitparams) override;
  ErrCode initPar2(FitParams* fitparams) override;
  ErrCode initCov(FitParams* fitparams) const override;

  ErrCode filter(const FitParams& reference, FitParams& params) const override;

private:
  ErrCode projectPhotonConstraint(const FitParams& fitparams, Projection& p) const;
  void initMomentum(FitParams* fitparams, const Eigen::Vector3d& origin);

  Eigen::Vector3d m_clusterPos;
  double m_clusterEnergy{0};
  Eigen::Matrix4d m_clusterCov;  ///< (x, y, z, E)
  int m_i1{2}, m_i2{0}, m_i3{1};  ///< running coordinate and the other two
};

} // namespace k4DTF
