// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DecayTreeFitter.h"

#include <stdexcept>

#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

namespace k4DTF {

// ---------------------------------------------------------------------------
// Internal state
// ---------------------------------------------------------------------------

struct DecayTreeFitter::FitState {
  TVectorD x;          // state vector: (px,py,pz,E) × nParticles
  TMatrixDSym C;       // covariance of x
  TVectorD x0;         // initial (measured) state
  TMatrixDSym C0;      // initial covariance
  TVectorD lambda;     // Lagrange multipliers
  TVectorD h;          // constraint residuals h(x) = 0
  TMatrixD H;          // Jacobian dh/dx
  double chi2{0.0};
  int nConstraints{0};
  int nParams{0};
};

// ---------------------------------------------------------------------------

DecayTreeFitter::DecayTreeFitter(DecayDescriptor decayDescriptor)
    : m_desc(std::move(decayDescriptor)) {}

DecayTreeFitter& DecayTreeFitter::addMassConstraint(int pdg, double massGeV) {
  m_massConstraints[pdg] = massGeV;
  return *this;
}

DecayTreeFitter& DecayTreeFitter::setVertexConstraint(bool enable) {
  m_vertexConstraint = enable;
  return *this;
}

DecayTreeFitter& DecayTreeFitter::setMaxIterations(int n) {
  m_maxIter = n;
  return *this;
}

DecayTreeFitter& DecayTreeFitter::setConvergenceCriterion(double eps) {
  m_eps = eps;
  return *this;
}

// ---------------------------------------------------------------------------
// Main fit loop
// ---------------------------------------------------------------------------

FitResult DecayTreeFitter::fit(const std::vector<DTFParticle>& particles) const {
  const auto leafPdgs = m_desc.finalStatePdgs();
  if (particles.size() != leafPdgs.size()) {
    throw std::invalid_argument(
        "DecayTreeFitter::fit: number of input particles (" +
        std::to_string(particles.size()) +
        ") does not match the number of final-state nodes in the descriptor (" +
        std::to_string(leafPdgs.size()) + ")");
  }

  // TODO: implement the fit loop:
  //
  //   FitState state = initialiseState(particles);
  //   for (int iter = 0; iter < m_maxIter; ++iter) {
  //     if (!buildConstraints(state)) return extractResult(...Failed);
  //     const double deltaChi2 = iterate(state);
  //     if (std::abs(deltaChi2) < m_eps) return extractResult(...Success);
  //   }
  //   return extractResult(...MaxIterationsReached);

  return FitResult{};
}

// ---------------------------------------------------------------------------
// Internal helpers — implement these to complete the fitter
// ---------------------------------------------------------------------------

DecayTreeFitter::FitState DecayTreeFitter::initialiseState(
    const std::vector<DTFParticle>& particles) const {
  // TODO: build the state vector x0 and covariance C0 from the input particles.
  //
  // State vector layout: (px, py, pz, E) × nParticles (row-major).
  // Covariance layout:   block-diagonal, each 4×4 block from DTFParticle::cov
  //                      (10-element upper triangle, same order as
  //                       edm4hep::ReconstructedParticle::covMatrix).

  (void)particles;
  return FitState{};
}

bool DecayTreeFitter::buildConstraints(FitState& state) const {
  // TODO: populate state.h and state.H for each non-leaf node in m_desc.
  //
  // For each internal node:
  //   - 4-momentum conservation: sum(daughters p4) - mother p4 = 0
  //     (3 independent components; energy fixed by mass-shell if constrained)
  //   - mass constraint (if PDG is in m_massConstraints):
  //     E^2 - |p|^2 - m^2 = 0
  //   - vertex constraint (if m_vertexConstraint):
  //     requires extending the state vector with track/vertex parameters
  //
  // Set state.nConstraints to the total number of constraints.
  // Return false if the Jacobian H is singular.

  (void)state;
  return true;
}

double DecayTreeFitter::iterate(FitState& state) const {
  // TODO: perform one Lagrange-multiplier Newton step (Billoir / Moser-Roussarie):
  //
  //   W      = H C H^T                    (nConstraints × nConstraints)
  //   lambda = W^{-1} (H (x - x0) + h)   (Lagrange multipliers)
  //   x_new  = x0 - C H^T lambda          (updated state)
  //   chi2   = lambda^T W lambda
  //
  // Update state.x, state.C, state.chi2.
  // Return |chi2_new - chi2_old| for convergence testing.

  (void)state;
  return 0.0;
}

FitResult DecayTreeFitter::extractResult(const FitState& state, int nIter,
                                         FitResult::Status status) const {
  // TODO: unpack state.x and state.C back into a FitResult.
  //
  // Fitted particles should follow the same (px,py,pz,E) × nParticles layout
  // used in initialiseState(), with composite (non-leaf) nodes appended after
  // the final-state daughters.

  (void)state;
  (void)nIter;
  (void)status;
  return FitResult{};
}

}  // namespace k4DTF
