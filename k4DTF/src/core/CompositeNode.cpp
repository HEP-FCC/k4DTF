// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "core/CompositeNode.h"

#include "core/FitParams.h"
#include "core/Projection.h"

namespace k4DTF {

CompositeNode::CompositeNode(const DecayNodeInput& input, const ParticleNode* mother)
    : ParticleNode(input.name, input.pdgCode, mother) {
  const auto& meas = std::get<CompositeMeasurement>(input.measurement);
  m_hasEnergy = meas.hasEnergy;
  int dim = dimM();
  m_measurement.resize(dim);
  m_measurement.head<3>() = meas.vertex;
  m_measurement.segment<4>(3) = meas.momentum;
  if (!m_hasEnergy) {
    m_measurement.head<3>() = meas.vertex;
    m_measurement.segment<3>(3) = meas.momentum.head<3>();
    m_measurementCov = meas.cov.block(0, 0, dim, dim);
  } else {
    m_measurementCov = meas.cov;
  }
}

ErrCode CompositeNode::initPar1(FitParams* fitparams) {
  int posindex = posIndex();
  int momindex = momIndex();
  int indexmap[7];
  for (int i = 0; i < 3; ++i) indexmap[i] = posindex + i;
  for (int i = 0; i < 4; ++i) indexmap[i + 3] = momindex + i;

  for (int row = 0; row < dimM(); ++row) {
    fitparams->setPar(indexmap[row], m_measurement(row));
  }
  return ErrCode::success;
}

ErrCode CompositeNode::initPar2(FitParams* fitparams) {
  return initTau(fitparams);
}

ErrCode CompositeNode::projectRecoComposite(const FitParams& fitparams,
                                             Projection& p) const {
  int posindex = posIndex();
  int momindex = momIndex();
  int indexmap[7];
  for (int i = 0; i < 3; ++i) indexmap[i] = posindex + i;
  for (int i = 0; i < 4; ++i) indexmap[i + 3] = momindex + i;

  for (int row = 0; row < dimM(); ++row) {
    p.r(row) = fitparams.par(indexmap[row]) - m_measurement(row);
    p.H(row, indexmap[row]) = 1;
    for (int col = 0; col <= row; ++col) {
      p.V(row, col) = m_measurementCov(row, col);
    }
  }
  return ErrCode::success;
}

ErrCode CompositeNode::filter(const FitParams& reference, FitParams& params) const {
  ErrCode status = applyConstraint(reference, params, dimM(),
                                   &CompositeNode::projectRecoComposite);
  if (mother() && lenIndex() >= 0) {
    status |= applyConstraint(reference, params, 3,
                              &ParticleNode::projectGeoConstraint);
  }
  return status;
}

double CompositeNode::chiSquare(const FitParams& fitparams) const {
  Projection p(fitparams.dim(), dimM());
  projectRecoComposite(fitparams, p);
  return p.chiSquare();
}

} // namespace k4DTF
