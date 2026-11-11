// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/FitResult.h"
#include "k4DTFFCCAnalyses/DTFInput.h"

/// RDataFrame-compatible functions for the decay-tree fitter.
///
/// Two levels:
///  - DTFRunner + fit* extractors: candidates from k4DTF::CandidateBuilder,
///    charged leaves with their tracks (helix) and an optional primary vertex.
///    This is the full decay-tree fit with vertices.
///  - DTF_*: the earlier one-call functions on a list of reconstructed
///    particles. They use the 4-momentum covariance of the particles only (no
///    tracks, so no vertices), which Delphes leaves empty.
namespace k4DTF {

using RVecRCP = ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>;
using RVecF = ROOT::VecOps::RVec<float>;

// ---------------------------------------------------------------------------
// One-call functions (momentum-only)
// ---------------------------------------------------------------------------

/// chi2 of the fit of `particles` (leaves in descriptor order); with
/// motherPdg/motherMass > 0 the head is mass-constrained. {-1} on failure.
RVecF DTF_chi2(const RVecRCP& particles, const std::string& decayDesc, int motherPdg = -1, float motherMass = -1.0f);
RVecF DTF_prob(const RVecRCP& particles, const std::string& decayDesc, int motherPdg = -1, float motherMass = -1.0f);
/// Fitted mass of the head.
RVecF DTF_mass(const RVecRCP& particles, const std::string& decayDesc, int motherPdg = -1, float motherMass = -1.0f);
/// Fitted particles: every node in descriptor pre-order (head first).
RVecRCP DTF_fittedParticles(const RVecRCP& particles, const std::string& decayDesc, int motherPdg = -1,
                            float motherMass = -1.0f);
float DTF_chi2_single(const RVecRCP& daughters, const std::string& decayDesc, int motherPdg = -1,
                      float motherMass = -1.0f);

// ---------------------------------------------------------------------------
// Decay-tree fit of candidates built with k4DTF::CandidateBuilder
// ---------------------------------------------------------------------------

using RVecFit = ROOT::VecOps::RVec<k4DTF::FitResult>;

/// Fits every candidate once; extract quantities with the fit* functions
/// below instead of refitting per output column.
///
/// Leaves are passed with the descriptor's mass hypotheses and, for charged
/// leaves, their track state and Bz; photons with their cluster
/// (photonClusterPerReco; the overloads without clusters fail on photons). Mass constraints use the nominal masses
/// (ParticleTable) and apply to the listed PDG codes and their charge
/// conjugates. Declare one instance per configuration at global scope in the
/// analysis script; operator() is const and thread safe.
class DTFRunner {
public:
  DTFRunner(const CandidateBuilder& builder, std::vector<int> massConstraintPdgs = {},
            bool primaryVertexConstraint = false, int maxIterations = 30);

  /// Fit; pv is used if the primary-vertex constraint is enabled.
  RVecFit operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, float Bz,
                     const edm4hep::VertexData& pv) const;

  /// Fit without a primary vertex.
  RVecFit operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, float Bz) const;

  /// Fit with photons: clusterPerReco from photonClusterPerReco().
  RVecFit operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, const Vec_cl& clusterPerReco,
                     float Bz, const edm4hep::VertexData& pv) const;
  RVecFit operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, const Vec_cl& clusterPerReco,
                     float Bz) const;

  /// Leaves of candidate j as DTFParticles (mass hypotheses, tracks and
  /// clusters attached).
  std::vector<DTFParticle> inputs(const Candidates& c, std::size_t j, const Vec_rp& recos,
                                  const Vec_ts& tsPerReco, float Bz, const Vec_cl& clusterPerReco = {}) const;

  /// Fitter input of candidate j.
  k4DTF::DecayNodeInput treeInput(const Candidates& c, std::size_t j, const Vec_rp& recos, const Vec_ts& tsPerReco,
                                  float Bz, const Vec_cl& clusterPerReco = {}) const;

private:
  DecayDescriptor m_desc;
  DecayDescriptor m_cc;
  std::vector<int> m_massConstraints;
  bool m_pv;
  k4DTF::DecayTreeFitter m_fitter;

  RVecFit run(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, const Vec_cl& clusterPerReco,
              float Bz, const std::optional<k4DTF::InteractionPointInput>& pv) const;
};

/// Fit status as the k4DTF::ErrCode bit pattern: 0 = success.
Vec_i fitStatus(const RVecFit& fits);
Vec_f fitChi2(const RVecFit& fits);  ///< kInvalid unless converged
Vec_i fitNdf(const RVecFit& fits);
Vec_i fitNIter(const RVecFit& fits);
Vec_f fitProb(const RVecFit& fits);  ///< kInvalid unless converged

/// Per-node fitted quantities (node index as in DecayDescriptor::flatten());
/// kInvalid unless the fit converged (and, for vertex/ctau, the node has them).
Vec_f fitNodeMass(const RVecFit& fits, int node);
Vec_f fitNodeMassErr(const RVecFit& fits, int node);
Vec_f fitNodeP(const RVecFit& fits, int node);
Vec_f fitNodeMomentum(const RVecFit& fits, int node, int coord);  ///< coord 0,1,2 = px,py,pz
Vec_f fitNodeDecayVertex(const RVecFit& fits, int node, int coord);
Vec_f fitNodeDecayVertexErr(const RVecFit& fits, int node, int coord);
Vec_f fitNodeCtau(const RVecFit& fits, int node);
Vec_f fitNodeCtauErr(const RVecFit& fits, int node);
Vec_f fitNodeDecayLength(const RVecFit& fits, int node);
Vec_f fitNodeDecayLengthErr(const RVecFit& fits, int node);

}  // namespace k4DTF
