// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTFFCCAnalyses/DTFFunctions.h"

#include "TMath.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/FitResult.h"

namespace k4DTF {

namespace {

k4DTF::DecayTreeFitter makeFitter(const std::string& decayDesc,
                                  int motherPdg, float motherMass) {
  k4DTF::DecayTreeFitter dtf(k4DTF::DecayDescriptor::fromString(decayDesc));
  if (motherPdg > 0 && motherMass > 0.0f)
    dtf.addMassConstraint(motherPdg, static_cast<double>(motherMass));
  return dtf;
}

std::vector<k4DTF::DTFParticle> toDTF(const RVecRCP& rps) {
  std::vector<k4DTF::DTFParticle> out;
  out.reserve(rps.size());
  for (const auto& rp : rps) out.push_back(k4DTF::DTFParticle::fromEDM4hep(rp));
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------

RVecF DTF_chi2(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg,
               float motherMass) {
  auto dtf = makeFitter(decayDesc, motherPdg, motherMass);
  const auto result = dtf.fit(toDTF(particles));
  if (!result.converged()) return {-1.0f};
  return {static_cast<float>(result.chi2)};
}

RVecF DTF_prob(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg,
               float motherMass) {
  auto dtf = makeFitter(decayDesc, motherPdg, motherMass);
  const auto result = dtf.fit(toDTF(particles));
  if (!result.converged() || result.ndf <= 0) return {-1.0f};
  return {static_cast<float>(TMath::Prob(result.chi2, result.ndf))};
}

RVecF DTF_mass(const RVecRCP& particles,
               const std::string& decayDesc,
               int motherPdg,
               float motherMass) {
  auto dtf = makeFitter(decayDesc, motherPdg, motherMass);
  const auto result = dtf.fit(toDTF(particles));
  if (!result.converged() || result.fittedParticles.empty()) return {-1.0f};
  return {static_cast<float>(result.fittedParticles.back().mass())};
}

RVecRCP DTF_fittedParticles(const RVecRCP& particles,
                             const std::string& decayDesc,
                             int motherPdg,
                             float motherMass) {
  auto dtf = makeFitter(decayDesc, motherPdg, motherMass);
  const auto result = dtf.fit(toDTF(particles));
  if (!result.converged()) return {};
  RVecRCP out;
  out.reserve(result.fittedParticles.size());
  for (const auto& fp : result.fittedParticles) out.push_back(fp.toEDM4hep());
  return out;
}

float DTF_chi2_single(const RVecRCP& daughters,
                      const std::string& decayDesc,
                      int motherPdg,
                      float motherMass) {
  const auto v = DTF_chi2(daughters, decayDesc, motherPdg, motherMass);
  return v.empty() ? -1.0f : v[0];
}

}  // namespace k4DTF
