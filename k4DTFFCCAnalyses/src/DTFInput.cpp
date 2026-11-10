// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTFFCCAnalyses/DTFInput.h"

#include <algorithm>
#include <cmath>
#include <functional>

#include "k4DTF/DTFParticle.h"
#include "k4DTF/ParticleTable.h"

namespace k4DTF {

namespace {

double p2(const edm4hep::ReconstructedParticleData& rp) {
  return double(rp.momentum.x) * rp.momentum.x + double(rp.momentum.y) * rp.momentum.y +
         double(rp.momentum.z) * rp.momentum.z;
}

bool hasTrackState(const edm4hep::TrackState& ts) { return ts.location != kNoTrack; }

bool sameTrack(const edm4hep::TrackState& a, const edm4hep::TrackState& b) {
  return a.D0 == b.D0 && a.Z0 == b.Z0 && a.phi == b.phi && a.omega == b.omega && a.tanLambda == b.tanLambda;
}

/// Mass of each leaf hypothesis (nominal mass from the descriptor).
std::vector<double> leafMasses(const DecayDescriptor& d) {
  std::vector<double> m;
  for (int n : d.leafNodes()) m.push_back(d.flatten()[n].mass);
  return m;
}

/// Invariant mass of leaves [first, first + n) of a candidate.
double massOfLeaves(const Vec_i& cand, int first, int n, const Vec_rp& recos, const std::vector<double>& masses) {
  double E = 0, px = 0, py = 0, pz = 0;
  for (int l = first; l < first + n; ++l) {
    const auto& rp = recos[cand[l]];
    const double m = masses[l] >= 0 ? masses[l] : rp.mass;
    E += std::sqrt(p2(rp) + m * m);
    px += rp.momentum.x;
    py += rp.momentum.y;
    pz += rp.momentum.z;
  }
  const double m2 = E * E - px * px - py * py - pz * pz;
  return m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2);
}

double coordOf(const edm4hep::Vector3d& v, int coord) { return coord == 0 ? v.x : (coord == 1 ? v.y : v.z); }

/// Range of related indices [begin, end) of an MC particle, clipped to the column.
template <typename F>
void forEachRelated(const Vec_i& rel, unsigned begin, unsigned end, F&& f) {
  for (unsigned k = begin; k < end && k < rel.size(); ++k) f(rel[k]);
}

std::vector<int> related(const Vec_i& rel, unsigned begin, unsigned end) {
  std::vector<int> out;
  forEachRelated(rel, begin, end, [&](int i) { out.push_back(i); });
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Per reconstructed particle
// ---------------------------------------------------------------------------

Vec_i recoToMC(const Vec_rp& recos, const Vec_i& linkFrom, const Vec_i& linkTo) {
  Vec_i out(recos.size(), -1);
  const auto n = std::min(linkFrom.size(), linkTo.size());
  for (std::size_t k = 0; k < n; ++k) {
    const int from = linkFrom[k];
    if (from >= 0 && static_cast<std::size_t>(from) < out.size() && out[from] < 0) out[from] = linkTo[k];
  }
  return out;
}

Vec_i truthPDG(const Vec_i& recoMC, const Vec_mc& mc) {
  Vec_i out(recoMC.size(), 0);
  for (std::size_t i = 0; i < recoMC.size(); ++i) {
    const int m = recoMC[i];
    if (m >= 0 && static_cast<std::size_t>(m) < mc.size()) out[i] = mc[m].PDG;
  }
  return out;
}

Vec_ts trackStatePerReco(const Vec_rp& recos, const Vec_i& rpTrackIdx, const Vec_trk& tracks,
                         const Vec_ts& trackStates) {
  Vec_ts out;
  out.reserve(recos.size());
  for (const auto& rp : recos) {
    edm4hep::TrackState ts{};
    ts.location = kNoTrack;
    if (rp.tracks_end > rp.tracks_begin && rp.tracks_begin < rpTrackIdx.size()) {
      const int t = rpTrackIdx[rp.tracks_begin];
      if (t >= 0 && static_cast<std::size_t>(t) < tracks.size()) {
        const auto& trk = tracks[t];
        if (trk.trackStates_end > trk.trackStates_begin && trk.trackStates_begin < trackStates.size())
          ts = trackStates[trk.trackStates_begin];
      }
    }
    out.push_back(ts);
  }
  return out;
}

float firstOr(const Vec_f& v, float def) { return v.empty() ? def : v[0]; }

// ---------------------------------------------------------------------------
// Candidates
// ---------------------------------------------------------------------------

CandidateBuilder::CandidateBuilder(const std::string& descriptor, MassWindows massWindows, bool includeCC,
                                   std::size_t maxCandidates)
    : m_desc(DecayDescriptor::fromString(descriptor)),
      m_cc(m_desc.chargeConjugate()),
      m_windows(std::move(massWindows)),
      m_includeCC(includeCC && !m_desc.isSelfConjugate()),
      m_max(maxCandidates) {}

Vec_cand CandidateBuilder::build(const DecayDescriptor& d, const Vec_rp& recos, const Vec_i& pid) const {
  const auto& flat = d.flatten();
  const auto masses = leafMasses(d);
  std::vector<Vec_cand> sub(flat.size());

  // Children have larger pre-order indices than their mother: walk backwards.
  for (int n = static_cast<int>(flat.size()) - 1; n >= 0; --n) {
    const auto& node = flat[n];
    auto& out = sub[n];

    if (node.isLeaf()) {
      const long q = std::lround(ParticleTable::charge(node.pdg));
      for (std::size_t i = 0; i < recos.size() && i < pid.size(); ++i) {
        if (pid[i] == node.pdg && std::lround(recos[i].charge) == q) out.push_back(Vec_i{static_cast<int>(i)});
      }
      continue;
    }

    // identical sibling subtrees: combine each set once (ordered by first reco index)
    std::vector<std::pair<int, int>> identical;
    for (std::size_t a = 0; a < node.children.size(); ++a)
      for (std::size_t b = a + 1; b < node.children.size(); ++b)
        if (d.toString(node.children[a]) == d.toString(node.children[b])) identical.emplace_back(a, b);

    std::vector<int> chosenFirst(node.children.size(), -1);
    Vec_i current;
    std::function<void(std::size_t)> combine = [&](std::size_t k) {
      if (out.size() >= m_max) return;
      if (k == node.children.size()) {
        if (const auto w = m_windows.find(std::abs(node.pdg)); w != m_windows.end()) {
          // current holds the leaves of this node; shift masses to the node's leaf range
          std::vector<double> m(masses.begin() + node.firstLeaf, masses.begin() + node.firstLeaf + node.nLeaves);
          const double mass = massOfLeaves(current, 0, node.nLeaves, recos, m);
          if (mass < w->second.first || mass > w->second.second) return;
        }
        out.push_back(current);
        return;
      }
      for (const auto& s : sub[node.children[k]]) {
        bool ok = true;
        for (int idx : s)
          if (std::find(current.begin(), current.end(), idx) != current.end()) ok = false;
        for (const auto& [a, b] : identical)
          if (static_cast<std::size_t>(b) == k && s.front() <= chosenFirst[a]) ok = false;
        if (!ok) continue;
        const auto size = current.size();
        current.insert(current.end(), s.begin(), s.end());
        chosenFirst[k] = s.front();
        combine(k + 1);
        current.erase(current.begin() + size, current.end());
      }
    };
    combine(0);
  }
  return sub[0];
}

Candidates CandidateBuilder::operator()(const Vec_rp& recos, const Vec_i& pid) const {
  Candidates c;
  for (const auto& cand : build(m_desc, recos, pid)) {
    c.leaves.push_back(cand);
    c.isCC.push_back(0);
  }
  if (m_includeCC) {
    for (const auto& cand : build(m_cc, recos, pid)) {
      c.leaves.push_back(cand);
      c.isCC.push_back(1);
    }
  }
  return c;
}

int nCandidates(const Candidates& c) { return static_cast<int>(c.leaves.size()); }

Vec_i leafRecoIndex(const Candidates& c, int leaf) {
  Vec_i out;
  out.reserve(c.leaves.size());
  for (const auto& cand : c.leaves) out.push_back(cand.at(leaf));
  return out;
}

Vec_rp leafParticles(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco,
                     const CandidateBuilder& b, int leaf, float Bz) {
  Vec_rp out;
  out.reserve(c.leaves.size());
  for (std::size_t j = 0; j < c.leaves.size(); ++j) {
    const auto& d = c.isCC[j] ? b.ccDescriptor() : b.descriptor();
    const auto& node = d.flatten()[d.leafNodes().at(leaf)];
    const int idx = c.leaves[j].at(leaf);
    auto rp = recos[idx];
    const double m = node.mass >= 0 ? node.mass : rp.mass;
    rp.PDG = node.pdg;
    rp.mass = static_cast<float>(m);
    rp.energy = static_cast<float>(std::sqrt(p2(rp) + m * m));
    const bool noCov = std::all_of(rp.covMatrix.values.begin(), rp.covMatrix.values.end(),
                                   [](float v) { return v == 0.f; });
    if (noCov && static_cast<std::size_t>(idx) < tsPerReco.size() && hasTrackState(tsPerReco[idx])) {
      const auto cov = DTFTrack::fromEDM4hep(tsPerReco[idx], Bz).p4Covariance(m);
      for (int i = 0; i < 10; ++i) rp.covMatrix[i] = static_cast<float>(cov[i]);
    }
    out.push_back(rp);
  }
  return out;
}

Vec_ts leafTrackStates(const Candidates& c, const Vec_ts& tsPerReco, int leaf) {
  Vec_ts out;
  out.reserve(c.leaves.size());
  for (const auto& cand : c.leaves) out.push_back(tsPerReco.at(cand.at(leaf)));
  return out;
}

Vec_f nodeMass(const Candidates& c, const Vec_rp& recos, const CandidateBuilder& b, int node) {
  const auto& f = b.descriptor().flatten().at(node);
  const auto masses = leafMasses(b.descriptor());  // identical for the CC descriptor
  Vec_f out;
  out.reserve(c.leaves.size());
  for (const auto& cand : c.leaves)
    out.push_back(static_cast<float>(massOfLeaves(cand, f.firstLeaf, f.nLeaves, recos, masses)));
  return out;
}

Vec_i nLeavesInPV(const Candidates& c, const Vec_ts& tsPerReco, const Vec_ts& pvTracks) {
  Vec_i out;
  out.reserve(c.leaves.size());
  for (const auto& cand : c.leaves) {
    int n = 0;
    for (int idx : cand) {
      const auto& ts = tsPerReco.at(idx);
      if (!hasTrackState(ts)) continue;
      if (std::any_of(pvTracks.begin(), pvTracks.end(), [&](const auto& pv) { return sameTrack(ts, pv); })) ++n;
    }
    out.push_back(n);
  }
  return out;
}

// ---------------------------------------------------------------------------
// MC truth
// ---------------------------------------------------------------------------

TruthMatcher::TruthMatcher(const CandidateBuilder& builder)
    : m_desc(builder.descriptor()), m_cc(builder.ccDescriptor()) {}

TruthMatch TruthMatcher::operator()(const Candidates& c, const Vec_i& recoMC, const Vec_mc& mc,
                                    const Vec_i& mcParents, const Vec_i& mcDaughters) const {
  TruthMatch t;
  const auto nMC = static_cast<int>(mc.size());
  auto valid = [&](int i) { return i >= 0 && i < nMC; };
  auto parentsOf = [&](int i) { return related(mcParents, mc[i].parents_begin, mc[i].parents_end); };
  auto daughtersOf = [&](int i) { return related(mcDaughters, mc[i].daughters_begin, mc[i].daughters_end); };

  for (std::size_t j = 0; j < c.leaves.size(); ++j) {
    const auto& d = c.isCC[j] ? m_cc : m_desc;
    const auto& flat = d.flatten();
    Vec_i nm(flat.size(), -1);

    for (int n = static_cast<int>(flat.size()) - 1; n >= 0; --n) {
      const auto& node = flat[n];
      if (node.isLeaf()) {
        const int idx = c.leaves[j].at(node.leafIndex);
        const int m = static_cast<std::size_t>(idx) < recoMC.size() ? recoMC[idx] : -1;
        if (valid(m) && mc[m].PDG == node.pdg) nm[n] = m;
        continue;
      }
      std::vector<int> kids;
      for (int ch : node.children) kids.push_back(nm[ch]);
      if (std::any_of(kids.begin(), kids.end(), [&](int k) { return !valid(k); })) continue;

      for (int q : parentsOf(kids.front())) {
        if (!valid(q) || mc[q].PDG != node.pdg) continue;
        const auto dau = daughtersOf(q);
        const bool allKids = std::all_of(kids.begin(), kids.end(), [&](int k) {
          return std::find(dau.begin(), dau.end(), k) != dau.end();
        });
        const bool onlyPhotonsElse = std::all_of(dau.begin(), dau.end(), [&](int x) {
          return std::find(kids.begin(), kids.end(), x) != kids.end() || (valid(x) && mc[x].PDG == 22);
        });
        if (allKids && onlyPhotonsElse) {
          nm[n] = q;
          break;
        }
      }
    }

    // follow mixing / copies of the head back to the produced particle
    int origin = nm[0];
    for (int guard = 0; valid(origin) && guard < 10; ++guard) {
      const auto par = parentsOf(origin);
      if (par.size() != 1 || !valid(par[0])) break;
      const int q = par[0];
      if (std::abs(mc[q].PDG) != std::abs(mc[origin].PDG) || daughtersOf(q).size() != 1) break;
      origin = q;
    }

    const int nMatched = static_cast<int>(std::count_if(nm.begin(), nm.end(), [](int x) { return x >= 0; }));
    t.nodeMC.push_back(nm);
    t.nMatched.push_back(nMatched);
    t.isSignal.push_back(nMatched == static_cast<int>(flat.size()) ? 1 : 0);
    t.headOrigin.push_back(origin);
  }
  return t;
}

Vec_mc nodeMCParticle(const TruthMatch& t, const Vec_mc& mc, int node) {
  Vec_mc out;
  out.reserve(t.nodeMC.size());
  for (const auto& nm : t.nodeMC) {
    const int m = nm.at(node);
    out.push_back(m >= 0 && static_cast<std::size_t>(m) < mc.size() ? mc[m] : edm4hep::MCParticleData{});
  }
  return out;
}

Vec_d nodeTrueProductionVertex(const TruthMatch& t, const Vec_mc& mc, int node, int coord) {
  Vec_d out;
  out.reserve(t.nodeMC.size());
  for (std::size_t j = 0; j < t.nodeMC.size(); ++j) {
    const int m = node == 0 ? t.headOrigin[j] : t.nodeMC[j].at(node);
    out.push_back(m >= 0 && static_cast<std::size_t>(m) < mc.size() ? coordOf(mc[m].vertex, coord) : kInvalid);
  }
  return out;
}

Vec_d nodeTrueDecayVertex(const TruthMatch& t, const Vec_mc& mc, const Vec_i& mcDaughters, int node, int coord) {
  Vec_d out;
  out.reserve(t.nodeMC.size());
  for (const auto& nm : t.nodeMC) {
    double v = kInvalid;
    const int m = nm.at(node);
    if (m >= 0 && static_cast<std::size_t>(m) < mc.size()) {
      const auto dau = related(mcDaughters, mc[m].daughters_begin, mc[m].daughters_end);
      if (!dau.empty() && dau.front() >= 0 && static_cast<std::size_t>(dau.front()) < mc.size())
        v = coordOf(mc[dau.front()].vertex, coord);
    }
    out.push_back(v);
  }
  return out;
}

Vec_d nodeTrueCtau(const TruthMatch& t, const Vec_mc& mc, const Vec_i& mcDaughters, int node) {
  Vec_d out;
  out.reserve(t.nodeMC.size());
  Vec_d dv[3], pv[3];
  for (int k = 0; k < 3; ++k) {
    dv[k] = nodeTrueDecayVertex(t, mc, mcDaughters, node, k);
    pv[k] = nodeTrueProductionVertex(t, mc, node, k);
  }
  for (std::size_t j = 0; j < t.nodeMC.size(); ++j) {
    const int m = t.nodeMC[j].at(node);
    if (m < 0 || dv[0][j] == kInvalid || pv[0][j] == kInvalid) {
      out.push_back(kInvalid);
      continue;
    }
    const auto& p = mc[m].momentum;
    const double pmag = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    const double L = std::sqrt(std::pow(dv[0][j] - pv[0][j], 2) + std::pow(dv[1][j] - pv[1][j], 2) +
                               std::pow(dv[2][j] - pv[2][j], 2));
    out.push_back(pmag > 0 ? L * mc[m].mass / pmag : kInvalid);
  }
  return out;
}

Vec_ts leafTrueTrackState(const TruthMatch& t, const Vec_mc& mc, const CandidateBuilder& b, int leaf, float Bz) {
  const int node = b.descriptor().leafNodes().at(leaf);
  Vec_ts out;
  out.reserve(t.nodeMC.size());
  for (const auto& nm : t.nodeMC) {
    edm4hep::TrackState ts{};
    ts.location = kNoTrack;
    const int m = nm.at(node);
    if (m >= 0 && static_cast<std::size_t>(m) < mc.size() && mc[m].charge != 0) {
      const auto& p = mc[m];
      ts = DTFTrack::fromPositionMomentum({p.vertex.x, p.vertex.y, p.vertex.z},
                                          {p.momentum.x, p.momentum.y, p.momentum.z}, p.charge, Bz)
               .toEDM4hep();
    }
    out.push_back(ts);
  }
  return out;
}

double truePrimaryVertex(const Vec_mc& mc, int coord) {
  for (const auto& p : mc)
    if (p.PDG == 23) return coordOf(p.vertex, coord);
  for (const auto& p : mc)
    if (p.parents_begin == p.parents_end) return coordOf(p.vertex, coord);
  return 0.0;
}

}  // namespace k4DTF
