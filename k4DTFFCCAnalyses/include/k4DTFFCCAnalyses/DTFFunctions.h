// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"

/// RDataFrame-compatible free functions for the decay-tree fitter.
///
/// All functions are designed to be used as RDataFrame Define/Filter lambdas:
/// @code
///   df.Define("DTF_chi2",
///       k4DTF::DTF_chi2,
///       {"B0_daughters", "B0 -> K+ pi-", 511, 5.27963})
/// @endcode
///
/// Because RDataFrame passes column values positionally, the recommended
/// pattern is to wrap each call in a lambda capturing the configuration:
/// @code
///   const std::string desc = "B0 -> K+ pi-";
///   df.Define("B0_DTF_chi2", [desc](const RVec<ReconstructedParticleData>& v){
///       return k4DTF::DTF_chi2(v, desc, 511, 5.27963);
///   }, {"B0_daughters"});
/// @endcode
namespace k4DTF {

using RVecRCP = ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>;
using RVecF = ROOT::VecOps::RVec<float>;

// ---------------------------------------------------------------------------
// Per-combination outputs (one entry per candidate)
// ---------------------------------------------------------------------------

/// Run the DTF and return chi2 for each candidate.
/// @param particles    Final-state daughters (order must match descriptor).
/// @param decayDesc    Decay descriptor string, e.g. "B0 -> K+ pi-".
/// @param motherPdg    PDG code of the head node for mass constraint (-1 = none).
/// @param motherMass   Mass [GeV] for the head-node constraint (ignored if motherPdg < 0).
RVecF DTF_chi2(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg = -1,
               float motherMass = -1.0f);

/// Run the DTF and return chi2 probability (TMath::Prob) for each candidate.
RVecF DTF_prob(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg = -1,
               float motherMass = -1.0f);

/// Return the fitted invariant mass of the mother particle for each candidate.
RVecF DTF_mass(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg = -1,
               float motherMass = -1.0f);

/// Return the fitted daughter particles after the DTF for each candidate.
RVecRCP DTF_fittedParticles(const RVecRCP& particles,
                             const std::string& decayDesc,
                             int motherPdg = -1,
                             float motherMass = -1.0f);

// ---------------------------------------------------------------------------
// Convenience single-candidate wrappers (no RVec loop)
// ---------------------------------------------------------------------------

/// Same as DTF_chi2 but for a single set of daughters.
float DTF_chi2_single(const RVecRCP& daughters,
                      const std::string& decayDesc,
                      int motherPdg = -1,
                      float motherMass = -1.0f);

}  // namespace k4DTF