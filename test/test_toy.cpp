// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

// Helix conventions and a toy-MC validation of the fitter: decays with
// exact vertices, tracks as helices smeared with a known correlated
// covariance. A correct fit gives unit pulls and chi2 with mean ndf.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <map>
#include <string>

#include "ToyDecays.h"
#include "k4DTF/CalorimeterResolution.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/core/Helix.h"

using Catch::Matchers::WithinAbs;
using namespace k4DTF;

// ---------------------------------------------------------------------------
// Helix
// ---------------------------------------------------------------------------

TEST_CASE("Helix from position and momentum is a consistent helix", "[helix]") {
  const double Bz = 2.0;
  // soft track (R ~ 0.3 m) produced far from the origin, so that the sagitta
  // distinguishes the two senses of rotation
  const Eigen::Vector3d x(12.0, -7.0, 30.0);
  const Eigen::Vector3d p(0.15, 0.08, -0.06);
  for (double q : {+1.0, -1.0}) {
    const auto par = helix::fromPositionMomentum(x, p, q, Bz);
    CHECK(helix::charge(par, Bz) == q);
    const Eigen::Vector3d pc = helix::momentum(par, Bz);
    CHECK_THAT(std::hypot(pc.x(), pc.y()), WithinAbs(std::hypot(p.x(), p.y()), 1e-12));
    CHECK_THAT(pc.z(), WithinAbs(p.z(), 1e-12));

    // a positive charge in Bz > 0 (omega > 0) turns clockwise seen from +z:
    // the centre is to the right of the momentum
    const double R = 1.0 / std::abs(par(helix::Omega));
    const Eigen::Vector3d pca = helix::pca(par);
    const double sgn = par(helix::Omega) > 0 ? 1.0 : -1.0;
    const double ux = pc.x() / std::hypot(pc.x(), pc.y()), uy = pc.y() / std::hypot(pc.x(), pc.y());
    const double cx = pca.x() + sgn * R * uy, cy = pca.y() - sgn * R * ux;
    CHECK_THAT(std::hypot(x.x() - cx, x.y() - cy), WithinAbs(R, 1e-6));

    // through the origin: D0 = Z0 = 0
    const auto par0 = helix::fromPositionMomentum(Eigen::Vector3d::Zero(), p, q, Bz);
    CHECK_THAT(par0(helix::D0), WithinAbs(0.0, 1e-12));
    CHECK_THAT(par0(helix::Z0), WithinAbs(0.0, 1e-12));
    CHECK_THAT(par0(helix::Phi), WithinAbs(std::atan2(p.y(), p.x()), 1e-12));

    // reference point: shifting both the point and the reference changes nothing
    const Eigen::Vector3d ref(0.3, -0.2, 1.0);
    const auto a = helix::fromPositionMomentum(x + ref, p, q, Bz, ref);
    const auto b = helix::fromPositionMomentum(x, p, q, Bz);
    CHECK((a - b).norm() < 1e-9);
  }
}

TEST_CASE("Helix momentum covariance matches numerical propagation", "[helix]") {
  const double Bz = 2.0;
  const auto par = helix::fromPositionMomentum({0.1, 0.2, 0.3}, {2.0, -1.0, 0.5}, +1.0, Bz);
  const auto cov = toy::Generator::covariance({2.0, -1.0, 0.5}, Bz);
  const Eigen::Matrix3d analytic = helix::momentumCovariance(par, cov, Bz);
  Eigen::Matrix<double, 3, 5> J;
  for (int k = 0; k < 5; ++k) {
    const double h = 1e-7 * std::max(1e-6, std::abs(par(k)));
    auto up = par, dn = par;
    up(k) += h;
    dn(k) -= h;
    J.col(k) = (helix::momentum(up, Bz) - helix::momentum(dn, Bz)) / (2 * h);
  }
  const Eigen::Matrix3d numeric = J * cov * J.transpose();
  CHECK((analytic - numeric).norm() < 1e-6 * numeric.norm());
}

TEST_CASE("Analytic helix Jacobian matches numerical derivatives", "[helix]") {
  const double Bz = 2.0;
  const Eigen::Vector3d ref(0.1, -0.05, 0.2);
  // prompt, displaced (B/D), far displaced (KS), soft and backward tracks
  const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> cases = {
      {{0.0, 0.0, 0.0}, {1.5, -0.7, 2.0}},     {{0.8, -1.2, 3.0}, {12.0, 5.0, -8.0}},
      {{250.0, 310.0, -120.0}, {0.9, 1.1, 0.4}}, {{3.0, 2.0, -1.0}, {0.12, -0.08, 0.05}},
      {{-5.0, 1.0, 40.0}, {-2.0, -3.0, -6.0}}};
  for (double q : {+1.0, -1.0}) {
    for (const auto& [x, p] : cases) {
      const auto analytic = helix::jacobian(x, p, q, Bz, ref);
      Eigen::Matrix<double, 5, 6> numeric;
      Eigen::Matrix<double, 6, 1> v;
      v << x, p;
      for (int k = 0; k < 6; ++k) {
        const double h = 1e-6 * std::max(1.0, std::abs(v(k)));
        auto up = v, dn = v;
        up(k) += h;
        dn(k) -= h;
        helix::Vector5d d = helix::fromPositionMomentum(up.head<3>(), up.tail<3>(), q, Bz, ref) -
                            helix::fromPositionMomentum(dn.head<3>(), dn.tail<3>(), q, Bz, ref);
        d(helix::Phi) = helix::wrapPhi(d(helix::Phi));
        numeric.col(k) = d / (2 * h);
      }
      for (int i = 0; i < 5; ++i)
        for (int k = 0; k < 6; ++k) {
          INFO("q " << q << " x " << x.transpose() << " d(par " << i << ")/d(" << k << ")");
          CHECK_THAT(analytic(i, k), WithinAbs(numeric(i, k), 1e-6 * (1.0 + std::abs(numeric(i, k)))));
        }
    }
  }
}

TEST_CASE("Helix positions and momenta along the track", "[helix]") {
  const double Bz = 2.0;
  const Eigen::Vector3d ref(0.1, -0.05, 0.2);
  for (double q : {+1.0, -1.0}) {
    const Eigen::Vector3d x(310.0, -150.0, 80.0), p(0.8, 0.3, -0.4);
    const auto par = helix::fromPositionMomentum(x, p, q, Bz, ref);
    const double s = helix::arcLengthNear(par, x, ref);
    CHECK((helix::positionAt(par, s, ref) - x).norm() < 1e-9);
    CHECK((helix::momentumAt(par, Bz, s) - p).norm() < 1e-12);
    CHECK((helix::positionAt(par, 0.0, ref) - helix::pca(par, ref)).norm() < 1e-12);
    // the track moves outwards from x: the PCA to the reference is behind it
    CHECK(s > 0);
  }
}

TEST_CASE("Two-helix vertex seed", "[helix]") {
  const double Bz = 2.0;
  auto seed = [&](const Eigen::Vector3d& v, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
                  double tolerance = 1e-3) {
    Eigen::Vector3d out;
    REQUIRE(helix::intersect(helix::fromPositionMomentum(v, p1, +1, Bz), Eigen::Vector3d::Zero(),
                             helix::fromPositionMomentum(v, p2, -1, Bz), Eigen::Vector3d::Zero(), out, tolerance));
    return out;
  };
  // D0 a few mm from the beam, KS at 40 cm, and nearly parallel daughters
  // (crossing angle ~2 deg) of a prompt-ish decay and of a boosted KS
  const Eigen::Vector3d d0(1.2, -0.8, 2.0), ks(280.0, 290.0, -100.0);
  CHECK((seed(d0, {1.5, 0.4, 0.8}, {0.3, -1.1, 2.0}) - d0).norm() < 1e-6);
  CHECK((seed(ks, {0.9, 0.7, -0.2}, {0.4, 1.2, -0.3}) - ks).norm() < 1e-6);
  CHECK((seed(d0, {8.0, 0.0, 1.0}, {7.0, 0.25, 1.1}) - d0).norm() < 1e-6);
  CHECK((seed(ks, {5.0, 5.2, -1.8}, {4.5, 4.9, -1.6}) - ks).norm() < 1e-6);
  // Nearly parallel daughters seen with errors: the helices come within
  // 61 um of each other at the point nearest the beam, a 1.4 mm displaced
  // vertex is not resolved with a 0.1 mm tolerance and the seed is the
  // point nearest the beam. The fit takes it from there.
  CHECK((seed(d0, {8.0, 0.0, 1.0}, {7.0, 0.25, 1.1}, 0.1) - Eigen::Vector3d(0.0150, -0.8209, 1.8331)).norm() < 1e-3);
  // A boosted KS whose daughters also cross exactly 12 cm behind the beam
  // line: the crossing ahead of the beam is the vertex.
  const Eigen::Vector3d ks2(21.2, 390.3, 309.1);
  CHECK((seed(ks2, {0.2826107255, 2.383302578, 1.893404036}, {-0.04807012709, 1.899391814, 1.49834698}, 0.2) - ks2)
            .norm() < 1e-4);
}

TEST_CASE("Phi wrapping", "[helix]") {
  CHECK_THAT(helix::wrapPhi(3 * toy::kPi), WithinAbs(toy::kPi, 1e-12));
  CHECK_THAT(helix::wrapPhi(-toy::kPi), WithinAbs(toy::kPi, 1e-12));
  CHECK_THAT(helix::wrapPhi(0.5), WithinAbs(0.5, 1e-15));
  CHECK(!std::isfinite(helix::wrapPhi(INFINITY)));  // returns, no endless loop
}

// ---------------------------------------------------------------------------
// Toy MC
// ---------------------------------------------------------------------------

namespace {

struct Stat {
  double n = 0, s = 0, s2 = 0;
  void add(double x) {
    n += 1;
    s += x;
    s2 += x * x;
  }
  double mean() const { return s / n; }
  double rms() const { return std::sqrt(s2 / n - mean() * mean()); }
};

struct ToyResult {
  int nFit = 0, nConverged = 0, ndf = -1;
  Stat chi2, probBelow5;
  std::map<std::string, Stat> pulls;
};

/// Fit nEvents toy decays and collect pulls by "<node index>_<quantity>".
ToyResult runToy(toy::Particle (*make)(), bool withPV, int nEvents, unsigned seed,
                 const std::vector<int>& massConstrained, double maxDecayRadius = 1e9,
                 const CalorimeterResolution* calo = nullptr) {
  toy::Generator gen(seed);
  if (calo) gen.useDelphesCalorimeter(*calo);
  std::normal_distribution<double> n01(0, 1);
  const Eigen::Vector3d bs(4.5e-3, 2e-5, 0.3);  // luminous region [mm]
  ToyResult out;
  DecayTreeFitter dtf;
  for (int ev = 0; ev < nEvents; ++ev) {
    auto head = make();
    const Eigen::Vector3d pv(bs(0) * n01(gen.rng()), bs(1) * n01(gen.rng()), bs(2) * n01(gen.rng()));
    gen.generate(head, pv);
    // decays outside the tracker would not be reconstructed
    while (toy::maxDecayRadius(head) > maxDecayRadius) gen.generate(head, pv);
    std::optional<InteractionPointInput> ip;
    if (withPV) ip = InteractionPointInput{Eigen::Vector3d::Zero(), Eigen::Vector3d(bs.cwiseProduct(bs)).asDiagonal()};
    const auto r = dtf.fit(gen.input(head, massConstrained), ip);
    ++out.nFit;
    if (!r.converged()) continue;
    ++out.nConverged;
    out.ndf = r.ndf;
    out.chi2.add(r.chi2);
    out.probBelow5.add(r.prob() < 0.05 ? 1.0 : 0.0);

    std::vector<const toy::Particle*> truth;
    toy::flatten(head, truth);
    REQUIRE(truth.size() == r.nodes.size());
    for (std::size_t i = 0; i < truth.size(); ++i) {
      const auto& n = r.nodes[i];
      const auto& t = *truth[i];
      const std::string k = std::to_string(i) + "_";
      for (int c = 0; c < 3; ++c) out.pulls[k + "p" + "xyz"[c]].add((n.p4(c) - t.p4(c)) / std::sqrt(n.p4Cov(c, c)));
      if (n.hasDecayVertex && !t.resonance)
        for (int c = 0; c < 3; ++c)
          out.pulls[k + "vtx" + "xyz"[c]].add((n.decayVertex(c) - t.decay(c)) / std::sqrt(n.decayVertexCov(c, c)));
      if (n.hasLifetime) out.pulls[k + "ctau"].add((n.ctau - t.trueCtau) / n.ctauErr);
    }
  }
  return out;
}

void checkToy(const ToyResult& r, int expectedNdf) {
  INFO("converged " << r.nConverged << "/" << r.nFit);
  CHECK(r.nConverged >= 0.99 * r.nFit);
  CHECK(r.ndf == expectedNdf);
  // <chi2> = ndf, with sqrt(2 ndf / N) statistical error
  CHECK_THAT(r.chi2.mean(), WithinAbs(expectedNdf, 5 * std::sqrt(2.0 * expectedNdf / r.nConverged)));
  CHECK_THAT(r.probBelow5.mean(), WithinAbs(0.05, 0.03));
  for (const auto& [name, s] : r.pulls) {
    INFO("pull " << name << ": mean " << s.mean() << " rms " << s.rms());
    CHECK(std::abs(s.mean()) < 0.2);
    CHECK(std::abs(s.rms() - 1.0) < 0.12);
  }
}

}  // namespace

TEST_CASE("Toy B- -> (D0 -> K- pi+) K-: pulls and chi2", "[toy]") {
  // 3 helices (15) + D0 mass (1) - B vertex (3) - D0 flight (1) - momenta (9) = 3
  checkToy(runToy(toy::BuToD0K, false, 600, 1, {421}), 3);
}

TEST_CASE("Toy B- -> (D0 -> K- pi+) K- with primary vertex", "[toy]") {
  // + PV position (3) - B flight (1) = 5
  checkToy(runToy(toy::BuToD0K, true, 600, 2, {421}), 5);
}

TEST_CASE("Toy B0 -> (D*- -> (D~0 -> K+ pi-) pi-) pi+ with primary vertex", "[toy]") {
  // 4 helices (20) + PV (3) + D0 mass (1) - B vertex+flight (4) - D0 flight (1) - momenta (12) = 7
  checkToy(runToy(toy::BdToDstPi, true, 600, 3, {-421}), 7);
}

TEST_CASE("Toy B0 -> (K*0 -> K+ pi-) gamma with primary vertex: photons in all directions", "[toy][photon]") {
  // 2 helices (10) + photon (3) + PV (3) - B vertex+flight (4) - momenta (9) = 3
  // Photons cover |cos theta| < 0.8 around the B direction, so many are in the
  // barrel (pz ~ 0), where a px/pz, py/pz parametrisation would diverge.
  checkToy(runToy(toy::BdToKstGamma, true, 600, 4, {}), 3);
}

TEST_CASE("Toy B- -> (D0 -> K- pi+ (pi0 -> gamma gamma)) pi- with Delphes IDEA photons", "[toy][photon]") {
  // 3 helices (15) + 2 photons (6) + PV (3) + D0 and pi0 masses (2)
  //   - B vertex+flight (4) - D0 flight (1) - momenta (15) = 6
  // Photons are placed as Delphes does: uniformly within their 0.02 x 0.02
  // tower, with the errors of CalorimeterResolution.
  const CalorimeterResolution calo;
  checkToy(runToy(toy::BuToD0Pi_KPiPi0, true, 600, 6, {421, 111}, 1e9, &calo), 6);
}

TEST_CASE("Calorimeter resolution model", "[photon]") {
  CalorimeterResolution calo;
  CHECK_THAT(calo.angularResolution(), WithinAbs(0.02 / std::sqrt(6.0), 1e-12));
  calo.smearTowerCenter = false;
  CHECK_THAT(calo.angularResolution(), WithinAbs(0.02 / std::sqrt(12.0), 1e-12));
  calo.smearTowerCenter = true;
  CHECK_THAT(std::sqrt(calo.energyVariance(1.0)), WithinAbs(std::sqrt(0.005 * 0.005 + 0.03 * 0.03 + 0.002 * 0.002), 1e-12));

  // front face: barrel at 45 deg, endcap at small angle
  CHECK((calo.frontFace({1.0, 0.0, 1.0}) - Eigen::Vector3d(2250, 0, 2250)).norm() < 1e-9);
  CHECK((calo.frontFace({0.1, 0.0, -1.0}) - Eigen::Vector3d(250, 0, -2500)).norm() < 1e-9);

  // position covariance: nothing along the line of sight, r sigma(theta)
  // along theta and r sin(theta) sigma(phi) along phi
  const Eigen::Vector3d pos(2250.0, 0.0, 1000.0);
  const auto cov = calo.positionCovariance(pos);
  const Eigen::Vector3d n = pos.normalized();
  CHECK(std::abs(n.dot(cov * n)) < 1e-9);
  const double r = pos.norm(), st = 2250.0 / r, s = calo.angularResolution();
  CHECK_THAT(cov(1, 1), WithinAbs(std::pow(2250.0 * s, 2), 1e-6));
  const Eigen::Vector3d uTheta(1000.0 / r, 0.0, -st);
  CHECK_THAT(uTheta.dot(cov * uTheta), WithinAbs(std::pow(r * s * st, 2), 1e-6));
}

TEST_CASE("Toy B0 -> (J/psi -> mu+ mu-) (KS0 -> pi+ pi-) with primary vertex: KS up to 60 cm", "[toy][seed]") {
  // 4 helices (20) + PV (3) - B vertex+flight (4) - KS flight (1) - momenta (12) = 6
  // Without vertex seeding about 1% of these fits fail and others converge
  // to wrong minima (<chi2> ~ 4000).
  checkToy(runToy(toy::BdToJpsiKS, true, 600, 5, {}, 600.0), 6);
}
