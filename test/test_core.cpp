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