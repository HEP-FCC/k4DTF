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

using Catch::Matchers::WithinAbs;
using namespace k4dtf;
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
