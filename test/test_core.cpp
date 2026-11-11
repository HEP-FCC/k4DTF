// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include "k4DTF/core/DecayNodeInput.h"
#include "k4DTF/core/ParticleMeasurement.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/FitResult.h"
#include "k4DTF/ParticleTable.h"
#include "k4DTF/TreeInput.h"

using Catch::Matchers::WithinAbs;
using namespace k4DTF;
using namespace k4DTF;

// ---------------------------------------------------------------------------
// DTFParticle EDM4hep conversions
// ---------------------------------------------------------------------------

TEST_CASE("DTFParticle::fromEDM4hep maps momentum and basic fields", "[dtfparticle]") {
  edm4hep::ReconstructedParticleData rp{};
  rp.momentum = {1.0f, 2.0f, 3.0f};
  rp.energy   = 4.0f;
  rp.charge   = -1.0f;
  rp.PDG      = -211;

  const auto dp = DTFParticle::fromEDM4hep(rp);

  CHECK_THAT(dp.p4.Px(), WithinAbs(1.0, 1e-6));
  CHECK_THAT(dp.p4.Py(), WithinAbs(2.0, 1e-6));
  CHECK_THAT(dp.p4.Pz(), WithinAbs(3.0, 1e-6));
  CHECK_THAT(dp.p4.E(),  WithinAbs(4.0, 1e-6));
  CHECK(dp.charge == -1.0);
  CHECK(dp.pdg    == -211);
}

TEST_CASE("DTFParticle::fromEDM4hep copies all 10 covariance elements", "[dtfparticle]") {
  edm4hep::ReconstructedParticleData rp{};
  for (int i = 0; i < 10; ++i) rp.covMatrix[i] = static_cast<float>(i + 1);

  const auto dp = DTFParticle::fromEDM4hep(rp);

  for (int i = 0; i < 10; ++i)
    CHECK_THAT(dp.cov[i], WithinAbs(static_cast<double>(i + 1), 1e-6));
}

TEST_CASE("DTFParticle::toEDM4hep maps momentum and basic fields", "[dtfparticle]") {
  DTFParticle dp;
  dp.p4.SetPxPyPzE(1.0, 2.0, 3.0, 4.0);
  dp.charge = +1.0;
  dp.pdg    = 321;

  const auto rp = dp.toEDM4hep();

  CHECK_THAT(rp.momentum.x, WithinAbs(1.0f, 1e-5f));
  CHECK_THAT(rp.momentum.y, WithinAbs(2.0f, 1e-5f));
  CHECK_THAT(rp.momentum.z, WithinAbs(3.0f, 1e-5f));
  CHECK_THAT(rp.energy,     WithinAbs(4.0f, 1e-5f));
  CHECK(rp.charge == 1.0f);
  CHECK(rp.PDG    == 321);
}

TEST_CASE("DTFParticle::toEDM4hep fills mass from 4-momentum", "[dtfparticle]") {
  DTFParticle dp;
  const double mK = 0.4937, pz = 2.0;
  dp.p4.SetPxPyPzE(0.0, 0.0, pz, std::sqrt(mK * mK + pz * pz));

  const auto rp = dp.toEDM4hep();

  CHECK_THAT(rp.mass, WithinAbs(static_cast<float>(mK), 1e-4f));
}

TEST_CASE("DTFParticle::toEDM4hep copies all 10 covariance elements", "[dtfparticle]") {
  DTFParticle dp;
  for (int i = 0; i < 10; ++i) dp.cov[i] = static_cast<double>(i + 1);

  const auto rp = dp.toEDM4hep();

  for (int i = 0; i < 10; ++i)
    CHECK_THAT(rp.covMatrix[i], WithinAbs(static_cast<float>(i + 1), 1e-5f));
}

TEST_CASE("DTFParticle round-trip fromEDM4hep->toEDM4hep preserves all fields", "[dtfparticle]") {
  edm4hep::ReconstructedParticleData orig{};
  orig.momentum = {1.5f, -0.3f, 2.1f};
  orig.energy   = 3.0f;
  orig.charge   = +1.0f;
  orig.PDG      = 321;
  for (int i = 0; i < 10; ++i) orig.covMatrix[i] = static_cast<float>(i) * 0.1f;

  const auto result = DTFParticle::fromEDM4hep(orig).toEDM4hep();

  CHECK_THAT(result.momentum.x, WithinAbs(orig.momentum.x, 1e-5f));
  CHECK_THAT(result.momentum.y, WithinAbs(orig.momentum.y, 1e-5f));
  CHECK_THAT(result.momentum.z, WithinAbs(orig.momentum.z, 1e-5f));
  CHECK_THAT(result.energy,     WithinAbs(orig.energy,     1e-5f));
  CHECK(result.charge == orig.charge);
  CHECK(result.PDG    == orig.PDG);
  for (int i = 0; i < 10; ++i)
    CHECK_THAT(result.covMatrix[i], WithinAbs(orig.covMatrix[i], 1e-5f));
}

// ---------------------------------------------------------------------------
// DecayDescriptor
// ---------------------------------------------------------------------------

TEST_CASE("DecayDescriptor parses simple two-body decay", "[descriptor]") {
  const auto desc = DecayDescriptor::fromString("B0 -> K+ pi-");
  const auto pdgs = desc.finalStatePdgs();
  REQUIRE(pdgs.size() == 2);
  CHECK(pdgs[0] == 321);   // K+
  CHECK(pdgs[1] == -211);  // pi-
}

TEST_CASE("DecayDescriptor round-trips through toString", "[descriptor]") {
  const std::string input = "B0 -> K+ pi-";
  const auto desc = DecayDescriptor::fromString(input);
  CHECK_FALSE(desc.toString().empty());
}

TEST_CASE("DecayDescriptor rejects missing arrow", "[descriptor]") {
  CHECK_THROWS_AS(DecayDescriptor::fromString("B0 K+ pi-"),
                  std::invalid_argument);
}

TEST_CASE("DecayDescriptor parses nested decays", "[descriptor]") {
  const auto desc = DecayDescriptor::fromString("B- -> (D0 -> K- pi+) K-");
  const auto pdgs = desc.finalStatePdgs();
  REQUIRE(pdgs.size() == 3);
  CHECK(pdgs[0] == -321);
  CHECK(pdgs[1] == 211);
  CHECK(pdgs[2] == -321);

  const auto& root = desc.root();
  CHECK(root.pdg == -521);
  REQUIRE(root.daughters.size() == 2);
  CHECK(root.daughters[0].pdg == 421);
  CHECK(root.daughters[0].daughters.size() == 2);
  CHECK(root.daughters[1].isLeaf());
}

TEST_CASE("DecayDescriptor flattens in pre-order with contiguous leaf ranges", "[descriptor]") {
  const auto desc = DecayDescriptor::fromString("B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+");
  const auto& f = desc.flatten();
  REQUIRE(f.size() == 7);
  // 0 B0, 1 D*-, 2 D~0, 3 K+, 4 pi-, 5 pi- (soft), 6 pi+
  const int pdgs[7] = {511, -413, -421, 321, -211, -211, 211};
  for (int i = 0; i < 7; ++i) CHECK(f[i].pdg == pdgs[i]);
  CHECK(f[0].parent == -1);
  CHECK(f[1].parent == 0);
  CHECK(f[2].parent == 1);
  CHECK(f[5].parent == 1);
  CHECK(f[6].parent == 0);
  CHECK(f[0].nLeaves == 4);
  CHECK(f[1].firstLeaf == 0);
  CHECK(f[1].nLeaves == 3);
  CHECK(f[2].nLeaves == 2);
  CHECK(f[6].leafIndex == 3);
  CHECK(f[5].leafIndex == 2);
  CHECK(desc.nLeaves() == 4);
  CHECK(desc.leafNodes() == std::vector<int>{3, 4, 5, 6});
  CHECK_THAT(f[2].mass, WithinAbs(1.86484, 1e-5));
}

TEST_CASE("DecayDescriptor toString round-trips nested decays", "[descriptor]") {
  const auto a = DecayDescriptor::fromString("B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+");
  const auto b = DecayDescriptor::fromString(a.toString());
  CHECK(a.toString() == "B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+");
  CHECK(b.finalStatePdgs() == a.finalStatePdgs());
  CHECK(b.nNodes() == a.nNodes());
}

TEST_CASE("DecayDescriptor charge conjugation", "[descriptor]") {
  const auto d = DecayDescriptor::fromString("B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+");
  const auto cc = d.chargeConjugate();
  CHECK(cc.finalStatePdgs() == std::vector<int>{-321, 211, 211, -211});
  CHECK(cc.root().pdg == -511);
  CHECK_FALSE(d.isSelfConjugate());
  CHECK(DecayDescriptor::fromString("B0 -> (J/psi(1S) -> mu+ mu-) (KS0 -> pi+ pi-)").isSelfConjugate());
  // the CC string parses back
  CHECK(DecayDescriptor::fromString(cc.toString()).finalStatePdgs() == cc.finalStatePdgs());
}

TEST_CASE("DecayDescriptor syntax errors", "[descriptor]") {
  CHECK_THROWS_AS(DecayDescriptor::fromString("B0 -> (D0 -> K- pi+ pi+"), std::invalid_argument);
  CHECK_THROWS_AS(DecayDescriptor::fromString("B0 -> D0 -> K- pi+"), std::invalid_argument);
  CHECK_THROWS_AS(DecayDescriptor::fromString("B0 ->"), std::invalid_argument);
  CHECK_THROWS_AS(DecayDescriptor::fromString("B0 -> NotAParticle"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// ParticleTable
// ---------------------------------------------------------------------------

TEST_CASE("ParticleTable understands LHCb, EvtGen and other spellings", "[particles]") {
  using namespace ParticleTable;
  CHECK(pdgFromName("D~0") == -421);
  CHECK(pdgFromName("anti-D0") == -421);
  CHECK(pdgFromName("D0bar") == -421);
  CHECK(pdgFromName("D*(2010)-") == -413);
  CHECK(pdgFromName("D*-") == -413);
  CHECK(pdgFromName("KS0") == 310);
  CHECK(pdgFromName("K_S0") == 310);
  CHECK(pdgFromName("p~-") == -2212);
  CHECK(pdgFromName("anti-p-") == -2212);
  CHECK(pdgFromName("Lambda_c~-") == -4122);
  CHECK(pdgFromName("e+") == -11);
  CHECK(pdgFromName("mu-") == 13);
  CHECK(pdgFromName("J/psi(1S)") == 443);
  CHECK(pdgFromName("B_s~0") == -531);
  CHECK(pdgFromName("-421") == -421);
  CHECK_THROWS_AS(pdgFromName("Dzero"), std::invalid_argument);
  CHECK(chargeConjugate(310) == 310);
  CHECK(chargeConjugate(421) == -421);
  CHECK(charge(-211) == -1.0);
  CHECK(charge(11) == -1.0);
  CHECK(nameFromPdg(-421) == "D~0");
  CHECK(nameFromPdg(-413) == "D*(2010)-");
  CHECK_THAT(mass(-511), WithinAbs(5.27958, 1e-6));
}

// ---------------------------------------------------------------------------
// DTFTrack helpers
// ---------------------------------------------------------------------------

TEST_CASE("DTFTrack 4-momentum covariance matches numerical propagation", "[track]") {
  auto t = DTFTrack::fromPositionMomentum({0.1, 0.2, 0.3}, {2.0, -1.0, 0.5}, +1.0, 2.0);
  // correlated (phi, omega, tanLambda) covariance
  t.setCov(1, 1, 1e-6);
  t.setCov(2, 2, 1e-10);
  t.setCov(4, 4, 4e-6);
  t.setCov(2, 1, 2e-9);
  t.setCov(4, 2, -1e-9);
  const double m = 0.493677;
  const auto analytic = t.p4Covariance(m);

  auto p4of = [&](const DTFTrack& tr) {
    const auto p = tr.momentum();
    return std::array<double, 4>{p.X(), p.Y(), p.Z(), std::sqrt(p.Mag2() + m * m)};
  };
  const int params[3] = {1, 2, 4};
  double J[4][3];
  for (int k = 0; k < 3; ++k) {
    const double h = 1e-6 * std::max(1e-6, std::abs(t.par[params[k]]));
    auto up = t, dn = t;
    up.par[params[k]] += h;
    dn.par[params[k]] -= h;
    const auto a = p4of(up), b = p4of(dn);
    for (int i = 0; i < 4; ++i) J[i][k] = (a[i] - b[i]) / (2 * h);
  }
  int idx = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j <= i; ++j, ++idx) {
      double s = 0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) s += J[i][k] * t.covAt(params[k], params[l]) * J[j][l];
      CHECK_THAT(analytic[idx], WithinAbs(s, 1e-6 * std::abs(s) + 1e-14));
    }
  }
}

TEST_CASE("DTFParticle from track fills an empty Delphes covariance", "[track]") {
  const auto t = DTFTrack::fromPositionMomentum({0, 0, 0}, {1.0, 2.0, 3.0}, -1.0, 2.0);
  auto ts = t.toEDM4hep();
  for (int i = 0; i < 5; ++i) ts.covMatrix[i * (i + 1) / 2 + i] = 1e-6f;
  edm4hep::ReconstructedParticleData rp{};
  rp.momentum = {1.0f, 2.0f, 3.0f};
  rp.charge = -1.0f;
  rp.mass = 0.13957f;
  rp.energy = std::sqrt(14.0f + rp.mass * rp.mass);

  const auto dp = DTFParticle::fromEDM4hep(rp, ts, 2.0, 0.493677);
  REQUIRE(dp.hasTrack());
  CHECK_THAT(dp.mass(), WithinAbs(0.493677, 1e-6));
  CHECK(dp.cov[0] > 0);  // sigma^2(px)
  CHECK(dp.cov[9] > 0);  // sigma^2(E)
}

TEST_CASE("ParticleTable resonances", "[particles]") {
  CHECK(ParticleTable::isResonance(-413));
  CHECK(ParticleTable::isResonance(443));
  CHECK_FALSE(ParticleTable::isResonance(421));
  CHECK_FALSE(ParticleTable::isResonance(310));
  CHECK_FALSE(ParticleTable::isResonance(511));
}

TEST_CASE("makeDecayNodeInput builds the tree in descriptor order", "[treeinput]") {
  const auto desc = DecayDescriptor::fromString("B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+");
  std::vector<DTFParticle> leaves;
  for (int i = 0; i < 4; ++i) {
    auto t = DTFTrack::fromPositionMomentum({0, 0, 0}, {1.0 + i, 0.5, 0.2}, i == 0 || i == 3 ? 1 : -1, 2.0);
    for (int k = 0; k < 5; ++k) t.setCov(k, k, 1e-6);
    auto ts = t.toEDM4hep();
    edm4hep::ReconstructedParticleData rp{};
    rp.momentum = {float(1.0 + i), 0.5f, 0.2f};
    rp.charge = static_cast<float>(t.charge());
    rp.mass = 0.13957f;
    leaves.push_back(DTFParticle::fromEDM4hep(rp, ts, 2.0));
  }
  const auto in = makeDecayNodeInput(desc, leaves, {421});
  CHECK(in.pdgCode == 511);
  CHECK_FALSE(in.isResonance);
  REQUIRE(in.daughters.size() == 2);
  const auto& dst = in.daughters[0];
  CHECK(dst.pdgCode == -413);
  CHECK(dst.isResonance);
  CHECK_FALSE(dst.applyMassConstraint);
  const auto& d0 = dst.daughters[0];
  CHECK(d0.pdgCode == -421);
  CHECK(d0.applyMassConstraint);  // 421 also constrains its charge conjugate
  CHECK_THAT(d0.mass, WithinAbs(1.86484, 1e-5));
  const auto& k = d0.daughters[0];
  CHECK(k.pdgCode == 321);
  CHECK_THAT(k.mass, WithinAbs(0.493677, 1e-6));  // mass hypothesis from the descriptor
  REQUIRE(std::holds_alternative<TrackMeasurement>(k.measurement));
  CHECK(std::get<TrackMeasurement>(k.measurement).helix.has_value());

  CHECK_THROWS_AS(makeDecayNodeInput(desc, {leaves[0]}), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// DecayTreeFitter — core API tests
// ---------------------------------------------------------------------------

namespace {

DecayNodeInput makeTrack(const std::string& name, int pdg, double px, double py,
                          double pz, double sigma) {
  DecayNodeInput input;
  input.name = name;
  input.pdgCode = pdg;
  input.mass = std::abs(pdg) == 321 ? 0.493677 : 0.13957;

  TrackMeasurement meas;
  meas.momentum = Eigen::Vector3d(px, py, pz);
  meas.cov = Eigen::Matrix3d::Zero();
  meas.cov(0, 0) = sigma * sigma;
  meas.cov(1, 1) = sigma * sigma;
  meas.cov(2, 2) = sigma * sigma;
  meas.charge = (pdg > 0) ? 1.0 : -1.0;
  input.measurement = meas;

  return input;
}

DecayNodeInput makeTwoBodyDecay() {
  DecayNodeInput mother;
  mother.name = "B0";
  mother.pdgCode = 511;
  mother.daughters.push_back(
      makeTrack("K+", 321, 0.0, 0.0, 2.0, 0.01));
  mother.daughters.push_back(
      makeTrack("pi-", -211, 0.0, 0.0, -1.5, 0.01));

  ResonanceMeasurement motherMeas;
  motherMeas.momentum = Eigen::Vector4d(0.0, 0.0, 0.5, 3.6);
  motherMeas.cov = Eigen::Matrix4d::Identity() * 0.01;
  mother.measurement = motherMeas;

  return mother;
}

}  // namespace

TEST_CASE("DecayTreeFitter fits two-body decay without throwing", "[fitter]") {
  DecayTreeFitter dtf;
  auto input = makeTwoBodyDecay();
  FitResult result = dtf.fit(input);
  CHECK(result.nIter > 0);
  CHECK(result.chi2 >= 0.0);
}

TEST_CASE("DecayTreeFitter produces valid chi2", "[fitter]") {
  DecayTreeFitter dtf;
  auto input = makeTwoBodyDecay();
  FitResult result = dtf.fit(input);
  CHECK(result.chi2 >= 0.0);
}

TEST_CASE("DecayTreeFitter with interaction point", "[fitter]") {
  DecayTreeFitter dtf;
  auto input = makeTwoBodyDecay();

  InteractionPointInput ip;
  ip.position = Eigen::Vector3d(0.0, 0.0, 0.0);
  ip.cov = Eigen::Matrix3d::Identity() * 0.01;

  FitResult result = dtf.fit(input, ip);
  CHECK(failure(result.status) == false);
}

TEST_CASE("DecayTreeFitter needs a mass hypothesis for tracks", "[fitter]") {
  DecayTreeFitter dtf;
  auto input = makeTwoBodyDecay();
  input.daughters[0].mass = -1;
  CHECK(dtf.fit(input).status == ErrCode::badsetup);
}

TEST_CASE("DecayTreeFitter applies mass constraints", "[fitter]") {
  // momentum-only tracks: 6 measurements, 6 track parameters, no vertex
  // information, so the only constraint is the mass of the head
  DecayTreeFitter dtf;
  auto input = makeTwoBodyDecay();
  input.daughters[0] = makeTrack("K+", 321, 0.3, 0.1, 2.0, 0.01);
  input.daughters[1] = makeTrack("pi-", -211, -0.2, 0.05, -1.5, 0.01);
  const auto free = dtf.fit(input);
  REQUIRE(free.converged());
  CHECK_THAT(free.chi2, WithinAbs(0.0, 1e-8));

  input.applyMassConstraint = true;
  input.massConstraintValue = 3.5;
  const auto constrained = dtf.fit(input);
  REQUIRE(constrained.converged());
  CHECK(constrained.chi2 > 1.0);
  REQUIRE_FALSE(constrained.nodes.empty());
  CHECK_THAT(constrained.nodes[0].mass, WithinAbs(3.5, 1e-6));
  // no position measurement: a kinematic fit with one constraint, no vertex
  CHECK(constrained.ndf == 1);
  CHECK_FALSE(constrained.nodes[0].hasDecayVertex);
  // daughter masses are the hypotheses
  CHECK_THAT(constrained.nodes[1].mass, WithinAbs(0.493677, 1e-6));
  CHECK_THAT(constrained.nodes[2].mass, WithinAbs(0.13957, 1e-6));
}

TEST_CASE("FitResult stays valid after fit() returns", "[result]") {
  DecayTreeFitter dtf;
  FitResult result = dtf.fit(makeTwoBodyDecay());
  REQUIRE(result.converged());
  const FitResult copy = result;
  REQUIRE(copy.fitParams() != nullptr);
  REQUIRE(copy.nodes.size() == 3);
  CHECK(copy.nodes[0].children == std::vector<int>{1, 2});
  CHECK(copy.nodes[1].parent == 0);
  CHECK(copy.nodes[1].name == "K+");
}

// ---------------------------------------------------------------------------
// FitResult
// ---------------------------------------------------------------------------

TEST_CASE("FitResult chi2 probability", "[result]") {
  FitResult r;
  r.chi2 = 2.0;
  r.ndf = 2;
  CHECK_THAT(r.prob(), WithinAbs(std::exp(-1.0), 1e-9));
}

TEST_CASE("FitResult chi2 probability is 0 when ndf <= 0", "[result]") {
  FitResult r;
  r.ndf = 0;
  CHECK(r.prob() == 0.0);
}

TEST_CASE("FitResult default status is success", "[result]") {
  FitResult r;
  CHECK(failure(r.status) == false);
}
