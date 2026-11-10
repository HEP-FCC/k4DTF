// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "k4DTF/ChiSquare.h"
#include "k4DTF/ErrCode.h"
#include "core/DecayNodeInput.h"
#include "core/FitParams.h"
#include "core/KalmanCalculator.h"
#include "core/Projection.h"

namespace k4dtf {

class ParticleNode {
public:
  using ParticleVector = std::vector<ParticleNode*>;

  ParticleNode(const std::string& name, int pdg, const ParticleNode* mother);
  virtual ~ParticleNode();

  static std::unique_ptr<ParticleNode> create(const DecayNodeInput& input,
                                              const ParticleNode* mother);

  virtual int dim() const = 0;
  virtual int dimM() const = 0;

  virtual bool hasEnergy() const { return false; }
  virtual bool hasPosition() const { return false; }
  bool hasMassConstraint() const { return m_hasMassConstraint; }
  void setHasMassConstraint(bool v) { m_hasMassConstraint = v; }

  int index() const { return m_index; }
  virtual int posIndex() const { return m_posIndex; }
  virtual int momIndex() const { return m_momIndex; }
  virtual int lenIndex() const { return m_lenIndex; }

  const std::string& name() const { return m_name; }
  int pdgCode() const { return m_pdgCode; }
  double pdtMass() const { return m_pdtMass; }
  void setPdtMass(double m) { m_pdtMass = m; }
  int charge() const { return m_charge; }
  void setCharge(int c) { m_charge = c; }

  const ParticleNode* mother() const { return m_mother; }
  const ParticleVector& daughters() const { return m_daughters; }

  void setIndex(int i) { m_index = i; }
  void setPosIndex(int i) { m_posIndex = i; }
  void setMomIndex(int i) { m_momIndex = i; }
  void setLenIndex(int i) { m_lenIndex = i; }
  void setMother(const ParticleNode* m) { m_mother = m; }

  void addDaughter(std::unique_ptr<ParticleNode> dau);
  void updateIndex(int& offset);

  virtual ErrCode initPar1(FitParams* fitparams) = 0;
  virtual ErrCode initPar2(FitParams* fitparams) = 0;
  virtual ErrCode initCov(FitParams* fitparams) const;

  virtual ErrCode filter(const FitParams& reference, FitParams& params) const = 0;

  ChiSquare chiSquare(const FitParams& params) const;
  int nFinalChargedCandidates() const;

protected:
  ErrCode initTau(FitParams* fitparams) const;

  template <typename Derived>
  ErrCode applyConstraint(const FitParams& reference, FitParams& params,
                          int dimC,
                          ErrCode (Derived::*proj)(const FitParams &,
                                                   Projection &) const) const {
    Projection p(params.dim(), dimC);
    p.reset();
    const Derived *self = static_cast<const Derived *>(this);
    // Linearise the constraint at the reference (the result of the previous
    // iteration), so that the result does not depend on the order in which
    // constraints are applied, and evaluate it at the current parameters:
    //   r(x) ~ r(x_ref) + H (x - x_ref)
    ErrCode status = (self->*proj)(reference, p);
    if (!failure(status)) {
      p.r() += p.H() * (params.par() - reference.par());
      // projections fill the lower triangle of V only
      p.V() = p.V().template selfadjointView<Eigen::Lower>();
      // The gain uses the current covariance, which already contains the
      // constraints applied before this one.
      KalmanCalculator kalman;
      status |= kalman.init(p.r(), p.H(), params, &p.V());
      if (!failure(status)) {
        kalman.updatePar(params);
        kalman.updateCov(params);
        params.addChiSquare(kalman.chisq(), dimC, this);
      }
    }
    return status;
  }

public:
  ErrCode projectGeoConstraint(const FitParams &fitparams, Projection &p) const;
  ErrCode projectMassConstraint(const FitParams &fitparams, Projection &p) const;

  void collectVertexDaughters(ParticleVector& particles, int posindex);

private:
  std::string m_name;
  int m_pdgCode{0};
  int m_index{-1};
  int m_posIndex{-1};
  int m_momIndex{-1};
  int m_lenIndex{-1};
  double m_pdtMass{0};
  int m_charge{0};
  bool m_hasMassConstraint{false};

  const ParticleNode* m_mother{nullptr};
  ParticleVector m_daughters;
  std::vector<std::unique_ptr<ParticleNode>> m_daughters_owned;
};

} // namespace k4dtf
