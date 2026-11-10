// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTFFCCAnalyses/DTFFunctions.h"

#include <cmath>
#include <exception>

#include "k4DTF/ParticleTable.h"
#include "k4DTF/TreeInput.h"

namespace k4DTF {

// ---------------------------------------------------------------------------
// One-call functions (momentum-only)
// ---------------------------------------------------------------------------

namespace {

k4dtf::FitResult fitParticles(const RVecRCP& particles, const std::string& decayDesc, int motherPdg,
                              float motherMass) {
  const auto desc = DecayDescriptor::fromString(decayDesc);
  std::vector<DTFParticle> leaves;
  for (const auto& rp : particles) leaves.push_back(DTFParticle::fromEDM4hep(rp));
  auto input = makeDecayNodeInput(desc, leaves);
  if (motherPdg > 0 && motherMass > 0.0f) {
    input.applyMassConstraint = true;
    input.massConstraintValue = motherMass;
  }
  return k4dtf::DecayTreeFitter().fit(input);
}

edm4hep::ReconstructedParticleData toEDM4hep(const k4dtf::NodeResult& n, float charge) {
  edm4hep::ReconstructedParticleData rp{};
  rp.PDG = n.pdgCode;
  rp.momentum = {static_cast<float>(n.p4(0)), static_cast<float>(n.p4(1)), static_cast<float>(n.p4(2))};
  rp.energy = static_cast<float>(n.p4(3));
  rp.mass = static_cast<float>(n.mass);
  rp.charge = charge;
  for (int i = 0, k = 0; i < 4; ++i)
    for (int j = 0; j <= i; ++j, ++k) rp.covMatrix[k] = static_cast<float>(n.p4Cov(i, j));
  return rp;
}

}  // namespace

RVecF DTF_chi2(const RVecRCP& particles, const std::string& decayDesc, int motherPdg, float motherMass) {
  try {
    const auto r = fitParticles(particles, decayDesc, motherPdg, motherMass);
    if (r.converged()) return {static_cast<float>(r.chi2)};
  } catch (const std::exception&) {
  }
  return {-1.0f};
}

RVecF DTF_prob(const RVecRCP& particles, const std::string& decayDesc, int motherPdg, float motherMass) {
  try {
    const auto r = fitParticles(particles, decayDesc, motherPdg, motherMass);
    if (r.converged() && r.ndf > 0) return {static_cast<float>(r.prob())};
  } catch (const std::exception&) {
  }
  return {-1.0f};
}

RVecF DTF_mass(const RVecRCP& particles, const std::string& decayDesc, int motherPdg, float motherMass) {
  try {
    const auto r = fitParticles(particles, decayDesc, motherPdg, motherMass);
    if (r.converged() && !r.nodes.empty()) return {static_cast<float>(r.nodes.front().mass)};
  } catch (const std::exception&) {
  }
  return {-1.0f};
}

RVecRCP DTF_fittedParticles(const RVecRCP& particles, const std::string& decayDesc, int motherPdg,
                            float motherMass) {
  RVecRCP out;
  try {
    const auto r = fitParticles(particles, decayDesc, motherPdg, motherMass);
    if (!r.converged()) return out;
    for (const auto& n : r.nodes) out.push_back(toEDM4hep(n, static_cast<float>(ParticleTable::charge(n.pdgCode))));
  } catch (const std::exception&) {
  }
  return out;
}

float DTF_chi2_single(const RVecRCP& daughters, const std::string& decayDesc, int motherPdg, float motherMass) {
  const auto v = DTF_chi2(daughters, decayDesc, motherPdg, motherMass);
  return v.empty() ? -1.0f : v[0];
}

// ---------------------------------------------------------------------------
// Decay-tree fit of candidates built with k4DTF::CandidateBuilder
// ---------------------------------------------------------------------------

DTFRunner::DTFRunner(const CandidateBuilder& builder, std::vector<int> massConstraintPdgs,
                     bool primaryVertexConstraint, int maxIterations)
    : m_desc(builder.descriptor()),
      m_cc(builder.ccDescriptor()),
      m_massConstraints(std::move(massConstraintPdgs)),
      m_pv(primaryVertexConstraint) {
  k4dtf::DecayTreeFitter::Configuration cfg;
  cfg.maxIter = maxIterations;
  m_fitter.setConfiguration(cfg);
}

std::vector<DTFParticle> DTFRunner::inputs(const Candidates& c, std::size_t j, const Vec_rp& recos,
                                           const Vec_ts& tsPerReco, float Bz) const {
  const auto& d = c.isCC.at(j) ? m_cc : m_desc;
  const auto& flat = d.flatten();
  const auto leafNodes = d.leafNodes();
  const auto& cand = c.leaves.at(j);

  std::vector<DTFParticle> out;
  out.reserve(cand.size());
  for (std::size_t l = 0; l < cand.size(); ++l) {
    const auto& node = flat[leafNodes[l]];
    const int idx = cand[l];
    const auto& rp = recos.at(idx);
    const bool charged = static_cast<std::size_t>(idx) < tsPerReco.size() && tsPerReco[idx].location != kNoTrack;
    DTFParticle p = charged ? DTFParticle::fromEDM4hep(rp, tsPerReco[idx], Bz, node.mass)
                            : DTFParticle::fromEDM4hep(rp);
    if (!charged && node.mass >= 0) p.setMassHypothesis(node.mass);
    p.pdg = node.pdg;
    p.name = node.name;
    out.push_back(std::move(p));
  }
  return out;
}

k4dtf::DecayNodeInput DTFRunner::treeInput(const Candidates& c, std::size_t j, const Vec_rp& recos,
                                           const Vec_ts& tsPerReco, float Bz) const {
  return makeDecayNodeInput(c.isCC.at(j) ? m_cc : m_desc, inputs(c, j, recos, tsPerReco, Bz), m_massConstraints);
}

RVecFit DTFRunner::run(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, float Bz,
                       const std::optional<k4dtf::InteractionPointInput>& pv) const {
  RVecFit out;
  out.reserve(c.leaves.size());
  for (std::size_t j = 0; j < c.leaves.size(); ++j) {
    try {
      out.push_back(m_fitter.fit(treeInput(c, j, recos, tsPerReco, Bz), m_pv ? pv : std::nullopt));
    } catch (const std::exception&) {
      k4dtf::FitResult failed;
      failed.status = k4dtf::ErrCode::badsetup;
      out.push_back(std::move(failed));
    }
  }
  return out;
}

RVecFit DTFRunner::operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco, float Bz,
                              const edm4hep::VertexData& pv) const {
  return run(c, recos, tsPerReco, Bz, makeInteractionPoint(DTFVertex::fromEDM4hep(pv)));
}

RVecFit DTFRunner::operator()(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco,
                              float Bz) const {
  return run(c, recos, tsPerReco, Bz, std::nullopt);
}

// --- extraction ------------------------------------------------------------

namespace {

template <typename F>
Vec_f perFit(const RVecFit& fits, F&& f) {
  Vec_f out;
  out.reserve(fits.size());
  for (const auto& r : fits) out.push_back(r.converged() ? static_cast<float>(f(r)) : kInvalid);
  return out;
}

enum class Needs { Nothing, Vertex, Lifetime };

template <typename F>
Vec_f perNode(const RVecFit& fits, int node, Needs needs, F&& f) {
  Vec_f out;
  out.reserve(fits.size());
  for (const auto& r : fits) {
    bool ok = r.converged() && node >= 0 && static_cast<std::size_t>(node) < r.nodes.size();
    if (ok && needs == Needs::Vertex) ok = r.nodes[node].hasDecayVertex;
    if (ok && needs == Needs::Lifetime) ok = r.nodes[node].hasLifetime;
    out.push_back(ok ? static_cast<float>(f(r.nodes[node])) : kInvalid);
  }
  return out;
}

}  // namespace

Vec_i fitStatus(const RVecFit& fits) {
  Vec_i out;
  out.reserve(fits.size());
  for (const auto& r : fits) out.push_back(static_cast<int>(r.status));
  return out;
}

Vec_i fitNdf(const RVecFit& fits) {
  Vec_i out;
  out.reserve(fits.size());
  for (const auto& r : fits) out.push_back(r.ndf);
  return out;
}

Vec_i fitNIter(const RVecFit& fits) {
  Vec_i out;
  out.reserve(fits.size());
  for (const auto& r : fits) out.push_back(r.nIter);
  return out;
}

Vec_f fitChi2(const RVecFit& fits) {
  return perFit(fits, [](const k4dtf::FitResult& r) { return r.chi2; });
}

Vec_f fitProb(const RVecFit& fits) {
  return perFit(fits, [](const k4dtf::FitResult& r) { return r.prob(); });
}

Vec_f fitNodeMass(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Nothing, [](const k4dtf::NodeResult& n) { return n.mass; });
}

Vec_f fitNodeMassErr(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Nothing, [](const k4dtf::NodeResult& n) { return n.massErr; });
}

Vec_f fitNodeP(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Nothing, [](const k4dtf::NodeResult& n) { return n.p4.head<3>().norm(); });
}

Vec_f fitNodeMomentum(const RVecFit& fits, int node, int coord) {
  return perNode(fits, node, Needs::Nothing, [coord](const k4dtf::NodeResult& n) { return n.p4(coord); });
}

Vec_f fitNodeDecayVertex(const RVecFit& fits, int node, int coord) {
  return perNode(fits, node, Needs::Vertex, [coord](const k4dtf::NodeResult& n) { return n.decayVertex(coord); });
}

Vec_f fitNodeDecayVertexErr(const RVecFit& fits, int node, int coord) {
  return perNode(fits, node, Needs::Vertex,
                 [coord](const k4dtf::NodeResult& n) { return std::sqrt(n.decayVertexCov(coord, coord)); });
}

Vec_f fitNodeCtau(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Lifetime, [](const k4dtf::NodeResult& n) { return n.ctau; });
}

Vec_f fitNodeCtauErr(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Lifetime, [](const k4dtf::NodeResult& n) { return n.ctauErr; });
}

Vec_f fitNodeDecayLength(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Lifetime, [](const k4dtf::NodeResult& n) { return n.decayLength; });
}

Vec_f fitNodeDecayLengthErr(const RVecFit& fits, int node) {
  return perNode(fits, node, Needs::Lifetime, [](const k4dtf::NodeResult& n) { return n.decayLengthErr; });
}

}  // namespace k4DTF
