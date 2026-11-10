// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

// Candidate building and truth matching (k4DTFFCCAnalyses/DTFInput.h) on
// hand-made events in the k4SimDelphes EDM4hep layout.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

#include "SyntheticEvent.h"
#include "k4DTFFCCAnalyses/DTFFunctions.h"
#include "k4DTFFCCAnalyses/DTFInput.h"

using Catch::Matchers::WithinAbs;
using namespace k4DTF;
using namespace k4DTF::test;
using k4dtf::FitResult;

namespace {

const char* kBd = "B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+";
const CandidateBuilder::MassWindows kWindows = {{511, {4.9, 5.7}}, {413, {1.95, 2.07}}, {421, {1.80, 1.93}}};

struct Prepared {
  SyntheticEvent ev;
  BdChain chain{};
  Vec_i recoMC, pid;
  Vec_ts ts;
  std::vector<int> recoOf;  // MC index -> reco index
};

/// Event with one signal chain plus prompt background tracks and a photon.
Prepared makeEvent(bool cc, bool mixed) {
  Prepared p;
  const ROOT::Math::XYZPoint pv(0.003, -1e-5, 0.12);
  const int z = addZ(p.ev, pv);
  p.chain = addBdChain(p.ev, z, pv, cc, mixed);
  // prompt background
  const int bkg[4] = {p.ev.addMC(211, {1.0, 0.5, 2.0, std::sqrt(5.25 + 0.0195)}, pv, z),
                      p.ev.addMC(-211, {-2.0, 0.3, 1.0, std::sqrt(5.09 + 0.0195)}, pv, z),
                      p.ev.addMC(321, {0.4, -1.5, -0.5, std::sqrt(2.66 + 0.2437)}, pv, z),
                      p.ev.addMC(22, {0.5, 0.5, 0.5, std::sqrt(0.75)}, pv, z)};
  p.ev.finalise();

  p.recoOf.assign(p.ev.mc.size(), -1);
  for (int m : {p.chain.K, p.chain.piD0, p.chain.piSoft, p.chain.piB, bkg[0], bkg[1], bkg[2], bkg[3]})
    p.recoOf[m] = p.ev.addReco(m);

  p.recoMC = recoToMC(p.ev.recos, p.ev.linkFrom, p.ev.linkTo);
  p.pid = truthPDG(p.recoMC, p.ev.mc);
  p.ts = trackStatePerReco(p.ev.recos, p.ev.rpTracks, p.ev.tracks, p.ev.trackStates);
  return p;
}

/// RVec == RVec is element-wise (returns an RVec); compare as std::vector.
std::vector<int> asVector(const Vec_i& v) { return {v.begin(), v.end()}; }

int signalIndex(const TruthMatch& t) {
  for (std::size_t j = 0; j < t.isSignal.size(); ++j)
    if (t.isSignal[j]) return static_cast<int>(j);
  return -1;
}

}  // namespace

TEST_CASE("Per-reco columns are aligned with ReconstructedParticles", "[input]") {
  const auto p = makeEvent(false, false);
  REQUIRE(p.recoMC.size() == p.ev.recos.size());
  REQUIRE(p.ts.size() == p.ev.recos.size());
  CHECK(p.pid[p.recoOf[p.chain.K]] == 321);
  CHECK(p.pid[p.recoOf[p.chain.piB]] == 211);
  // the photon (last reco) has no track
  CHECK(p.ts.back().location == kNoTrack);
  CHECK(p.ts[p.recoOf[p.chain.K]].location == edm4hep::TrackState::AtIP);
}

TEST_CASE("CandidateBuilder finds the signal chain and TruthMatcher flags it", "[input]") {
  const auto p = makeEvent(false, false);
  const CandidateBuilder b(kBd, kWindows);
  const TruthMatcher matcher(b);

  const auto cands = b(p.ev.recos, p.pid);
  REQUIRE(nCandidates(cands) >= 1);
  const auto truth = matcher(cands, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  CHECK(std::count(truth.isSignal.begin(), truth.isSignal.end(), 1) == 1);

  const int j = signalIndex(truth);
  REQUIRE(j >= 0);
  CHECK(cands.isCC[j] == 0);
  CHECK(asVector(cands.leaves[j]) == std::vector<int>{p.recoOf[p.chain.K], p.recoOf[p.chain.piD0], p.recoOf[p.chain.piSoft],
                                 p.recoOf[p.chain.piB]});
  CHECK(asVector(truth.nodeMC[j]) ==
        std::vector<int>{p.chain.B, p.chain.Dst, p.chain.D0, p.chain.K, p.chain.piD0, p.chain.piSoft, p.chain.piB});
  CHECK(truth.headOrigin[j] == p.chain.B);

  // pre-fit masses with the descriptor's mass hypotheses
  CHECK_THAT(nodeMass(cands, p.ev.recos, b, 2)[j], WithinAbs(ParticleTable::mass(421), 2e-4));
  CHECK_THAT(nodeMass(cands, p.ev.recos, b, 1)[j], WithinAbs(ParticleTable::mass(413), 2e-4));
  CHECK_THAT(nodeMass(cands, p.ev.recos, b, 0)[j], WithinAbs(ParticleTable::mass(511), 5e-4));
}

TEST_CASE("Truth vertices, ctau and track parameters of the signal", "[input]") {
  const auto p = makeEvent(false, false);
  const CandidateBuilder b(kBd, kWindows);
  const auto cands = b(p.ev.recos, p.pid);
  const auto truth = TruthMatcher(b)(cands, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  const int j = signalIndex(truth);
  REQUIRE(j >= 0);

  const auto& mc = p.ev.mc;
  const auto& dau = p.ev.mcDaughters;
  CHECK_THAT(nodeTrueCtau(truth, mc, dau, 0)[j], WithinAbs(0.455, 1e-9));
  CHECK_THAT(nodeTrueCtau(truth, mc, dau, 1)[j], WithinAbs(0.0, 1e-12));  // D* is a resonance
  CHECK_THAT(nodeTrueCtau(truth, mc, dau, 2)[j], WithinAbs(0.123, 1e-9));
  CHECK_THAT(nodeTrueDecayVertex(truth, mc, dau, 2, 0)[j], WithinAbs(p.chain.dVertex.X(), 1e-12));
  CHECK_THAT(nodeTrueProductionVertex(truth, mc, 0, 2)[j], WithinAbs(p.chain.pv.Z(), 1e-12));
  CHECK_THAT(truePrimaryVertex(mc, 0), WithinAbs(p.chain.pv.X(), 1e-12));

  // reco tracks were made from the truth helix: identical parameters
  for (int leaf = 0; leaf < 4; ++leaf) {
    const auto reco = leafTrackStates(cands, p.ts, leaf)[j];
    const auto tru = leafTrueTrackState(truth, mc, b, leaf, p.ev.Bz)[j];
    CHECK(tru.location != kNoTrack);
    CHECK_THAT(reco.D0, WithinAbs(tru.D0, 1e-6));
    CHECK_THAT(reco.phi, WithinAbs(tru.phi, 1e-6));
    CHECK_THAT(reco.omega, WithinAbs(tru.omega, 1e-9));
    CHECK_THAT(reco.Z0, WithinAbs(tru.Z0, 1e-6));
    CHECK_THAT(reco.tanLambda, WithinAbs(tru.tanLambda, 1e-6));
  }
  // displaced D0 daughters have a non-zero impact parameter
  CHECK(std::abs(leafTrackStates(cands, p.ts, 0)[j].D0) > 1e-3);

  // leaf particles carry the mass hypothesis and a propagated covariance
  const auto K = leafParticles(cands, p.ev.recos, p.ts, b, 0, p.ev.Bz)[j];
  CHECK(K.PDG == 321);
  CHECK_THAT(K.mass, WithinAbs(ParticleTable::mass(321), 1e-6));
  CHECK(K.covMatrix[0] > 0.f);
  CHECK(nLeavesInPV(cands, p.ts, Vec_ts{})[j] == 0);
  CHECK(nLeavesInPV(cands, p.ts, Vec_ts{p.ts[p.recoOf[p.chain.piB]]})[j] == 1);
}

TEST_CASE("Charge-conjugate candidates", "[input]") {
  const auto p = makeEvent(true, false);
  const CandidateBuilder b(kBd, kWindows);
  const auto cands = b(p.ev.recos, p.pid);
  const auto truth = TruthMatcher(b)(cands, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  const int j = signalIndex(truth);
  REQUIRE(j >= 0);
  CHECK(cands.isCC[j] == 1);
  CHECK(leafParticles(cands, p.ev.recos, p.ts, b, 0, p.ev.Bz)[j].PDG == -321);

  // without CC the signal is not found
  const CandidateBuilder noCC(kBd, kWindows, false);
  const auto c2 = noCC(p.ev.recos, p.pid);
  const auto t2 = TruthMatcher(noCC)(c2, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  CHECK(signalIndex(t2) == -1);
}

TEST_CASE("Mixed B0: head origin is the produced particle", "[input]") {
  const auto p = makeEvent(false, true);
  const CandidateBuilder b(kBd, kWindows);
  const auto cands = b(p.ev.recos, p.pid);
  const auto truth = TruthMatcher(b)(cands, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  const int j = signalIndex(truth);
  REQUIRE(j >= 0);
  CHECK(truth.nodeMC[j][0] == p.chain.B);
  CHECK(truth.headOrigin[j] == p.chain.Bproduced);
  CHECK(p.chain.Bproduced != p.chain.B);
  CHECK_THAT(nodeTrueCtau(truth, p.ev.mc, p.ev.mcDaughters, 0)[j], WithinAbs(0.455, 1e-9));
}

TEST_CASE("Mass windows and identical siblings", "[input]") {
  const auto p = makeEvent(false, false);
  // a D0 window that excludes the D0 removes the signal
  const CandidateBuilder tight(kBd, {{421, {1.90, 1.95}}});
  const auto c = tight(p.ev.recos, p.pid);
  const auto t = TruthMatcher(tight)(c, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  CHECK(signalIndex(t) == -1);

  // two identical pi- leaves: each pair is used once
  SyntheticEvent ev;
  const ROOT::Math::XYZPoint o(0, 0, 0);
  const int z = addZ(ev, o);
  for (int pdg : {321, -211, -211, -211}) ev.addMC(pdg, {0.3 * pdg / 100.0, 1.0, 2.0, 4.0}, o, z);
  ev.finalise();
  for (int m = 3; m < 7; ++m) ev.addReco(m);
  const auto pid = truthPDG(recoToMC(ev.recos, ev.linkFrom, ev.linkTo), ev.mc);
  const CandidateBuilder dplus("D- -> K+ pi- pi-", {}, false);
  CHECK(nCandidates(dplus(ev.recos, pid)) == 3);  // C(3,2) pairs of pi-
}

TEST_CASE("DTFRunner fits the signal candidate", "[input]") {
  // tracks carry the exact helix of their MC particle, so the signal fits
  // with chi2 ~ 0 and reproduces the true vertices
  const auto p = makeEvent(false, false);
  const CandidateBuilder b(kBd, kWindows);
  const auto cands = b(p.ev.recos, p.pid);
  const auto truth = TruthMatcher(b)(cands, p.recoMC, p.ev.mc, p.ev.mcParents, p.ev.mcDaughters);
  const int j = signalIndex(truth);
  REQUIRE(j >= 0);
  const DTFRunner runner(b, {421}, true);

  const auto inputs = runner.inputs(cands, j, p.ev.recos, p.ts, p.ev.Bz);
  REQUIRE(inputs.size() == 4);
  for (const auto& in : inputs) {
    CHECK(in.hasTrack());
    CHECK(in.track->Bz == p.ev.Bz);
  }
  CHECK_THAT(inputs[0].mass(), WithinAbs(ParticleTable::mass(321), 1e-6));

  // primary vertex at the true position, beam-spot-like errors
  edm4hep::VertexData pv{};
  pv.position = {float(p.chain.pv.X()), float(p.chain.pv.Y()), float(p.chain.pv.Z())};
  pv.covMatrix[0] = 2e-5f;
  pv.covMatrix[2] = 1e-8f;
  pv.covMatrix[5] = 0.09f;
  const auto fits = runner(cands, p.ev.recos, p.ts, p.ev.Bz, pv);
  REQUIRE(fits.size() == cands.leaves.size());
  REQUIRE(fits[j].converged());
  CHECK(fitStatus(fits)[j] == 0);
  CHECK(fits[j].nodes.size() == 7);
  // 20 helix + 3 PV + 1 mass - 4 (B vertex, flight) - 1 (D0 flight) - 12 momenta
  CHECK(fitNdf(fits)[j] == 7);
  CHECK(fitChi2(fits)[j] < 1e-3);
  CHECK_THAT(fitNodeMass(fits, 2)[j], WithinAbs(ParticleTable::mass(421), 1e-5));
  CHECK_THAT(fitNodeCtau(fits, 2)[j], WithinAbs(0.123, 1e-4));
  CHECK_THAT(fitNodeCtau(fits, 0)[j], WithinAbs(0.455, 1e-4));
  CHECK_THAT(fitNodeDecayVertex(fits, 2, 0)[j], WithinAbs(p.chain.dVertex.X(), 1e-4));
  CHECK_THAT(fitNodeDecayVertex(fits, 0, 2)[j], WithinAbs(p.chain.bVertex.Z(), 1e-4));
  // the D* is a resonance: it decays at the B vertex, no lifetime
  CHECK_THAT(fitNodeDecayVertex(fits, 1, 1)[j], WithinAbs(p.chain.bVertex.Y(), 1e-4));
  CHECK(fitNodeCtau(fits, 1)[j] == kInvalid);
}
