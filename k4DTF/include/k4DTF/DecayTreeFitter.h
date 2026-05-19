// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <map>
#include <vector>

#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/FitResult.h"

namespace k4DTF {

/// Kinematic decay-tree fitter.
///
/// Implements a Lagrange-multiplier least-squares fit (following the Moser &
/// Roussarie approach) with optional per-node mass constraints and optional
/// vertex constraints.  The implementation in DecayTreeFitter.cpp is a
/// skeleton — fill in buildConstraintMatrix(), buildJacobian(), and iterate()
/// for a complete implementation.
///
/// Typical usage:
/// @code
///   DecayTreeFitter dtf(DecayDescriptor::fromString("B0 -> K+ pi-"));
///   dtf.addMassConstraint(511, 5.27963);  // B0 mass in GeV
///   FitResult r = dtf.fit(inputParticles);
///   if (r.converged()) { ... }
/// @endcode
class DecayTreeFitter {
public:
  explicit DecayTreeFitter(DecayDescriptor decayDescriptor);

  // -------------------------------------------------------------------------
  // Configuration (call before fit())
  // -------------------------------------------------------------------------

  /// Add a mass constraint for a node identified by its PDG code.
  DecayTreeFitter& addMassConstraint(int pdg, double massGeV);

  /// Enable/disable vertex constraint at the production vertex of each node.
  DecayTreeFitter& setVertexConstraint(bool enable);

  DecayTreeFitter& setMaxIterations(int n);
  DecayTreeFitter& setConvergenceCriterion(double eps);

  // -------------------------------------------------------------------------
  // Fit
  // -------------------------------------------------------------------------

  /// Fit the decay tree given the measured final-state particles.
  /// @param particles  Must match the leaf order of the DecayDescriptor.
  FitResult fit(const std::vector<DTFParticle>& particles) const;

private:
  DecayDescriptor m_desc;
  std::map<int, double> m_massConstraints;
  bool m_vertexConstraint{false};
  int m_maxIter{100};
  double m_eps{1e-6};

  // --- internal helpers (implement in DecayTreeFitter.cpp) -----------------

  struct FitState;  // forward declaration of internal state structure

  /// Initialise the state vector and covariance from input particles.
  FitState initialiseState(const std::vector<DTFParticle>& particles) const;

  /// Build the constraint vector h(x) and its Jacobian H = dh/dx.
  bool buildConstraints(FitState& state) const;

  /// Perform one Lagrange-multiplier Newton step; return delta-chi2.
  double iterate(FitState& state) const;

  /// Package the converged state into a FitResult.
  FitResult extractResult(const FitState& state, int nIter,
                          FitResult::Status status) const;
};

}  // namespace k4DTF