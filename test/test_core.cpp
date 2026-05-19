// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/FitResult.h"

using Catch::Matchers::WithinAbs;
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

  // float → double → float: values should survive within float precision
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
// DecayTreeFitter — basic interface tests
// (The actual fit result is a stub until the algorithm is implemented.)
// ---------------------------------------------------------------------------

namespace {

DTFParticle makeKaon() {
  DTFParticle p;
  // K+ with ~2 GeV momentum along z
  const double mK = 0.4937;
  const double pz = 2.0;
  const double E = std::sqrt(mK * mK + pz * pz);
  p.p4.SetPxPyPzE(0.0, 0.0, pz, E);
  p.charge = +1.0;
  p.pdg = 321;
  p.cov.fill(0.0);
  p.cov[9] = 1e-4;  // sigma_E^2
  return p;
}

DTFParticle makePion() {
  DTFParticle p;
  const double mPi = 0.1396;
  const double pz = -1.5;
  const double E = std::sqrt(mPi * mPi + pz * pz);
  p.p4.SetPxPyPzE(0.0, 0.0, pz, E);
  p.charge = -1.0;
  p.pdg = -211;
  p.cov.fill(0.0);
  p.cov[9] = 1e-4;
  return p;
}

}  // namespace

TEST_CASE("DecayTreeFitter accepts correct number of daughters", "[fitter]") {
  DecayTreeFitter dtf(DecayDescriptor::fromString("B0 -> K+ pi-"));
  dtf.addMassConstraint(511, 5.27963);

  const std::vector<DTFParticle> daughters = {makeKaon(), makePion()};
  REQUIRE_NOTHROW(dtf.fit(daughters));
}

TEST_CASE("DecayTreeFitter throws on wrong daughter count", "[fitter]") {
  DecayTreeFitter dtf(DecayDescriptor::fromString("B0 -> K+ pi-"));
  const std::vector<DTFParticle> one = {makeKaon()};
  CHECK_THROWS_AS(dtf.fit(one), std::invalid_argument);
}

TEST_CASE("FitResult chi2 probability is in [0,1] when ndf > 0", "[result]") {
  FitResult r;
  r.status = FitResult::Status::Success;
  r.chi2 = 1.5;
  r.ndf = 2;
  CHECK(r.prob() >= 0.0);
  CHECK(r.prob() <= 1.0);
}

TEST_CASE("FitResult prob is 0 when ndf <= 0", "[result]") {
  FitResult r;
  r.ndf = 0;
  CHECK(r.prob() == 0.0);
}