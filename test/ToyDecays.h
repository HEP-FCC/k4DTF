// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

/// Toy decays for testing the fitter core: exact kinematics and vertices,
/// tracks as helices smeared with a known (correlated) covariance.
///
/// Only Eigen and k4DTF/core are used, so the tests run without ROOT I/O.

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "k4DTF/core/DecayNodeInput.h"
#include "k4DTF/core/Helix.h"

namespace k4dtf::toy {

constexpr double kPi = 3.14159265358979323846;

struct Particle {
  int pdg{0};
  std::string name;
  double mass{0};
  double ctau{0};  ///< [mm]; 0 for resonances and stable particles
  double charge{0};
  bool resonance{false};
  std::vector<Particle> daughters;

  // filled by generate()
  Eigen::Vector4d p4{Eigen::Vector4d::Zero()};
  Eigen::Vector3d production{Eigen::Vector3d::Zero()};
  Eigen::Vector3d decay{Eigen::Vector3d::Zero()};
  double trueCtau{0};
};

inline Particle leaf(int pdg, const std::string& name, double mass, double charge) {
  Particle p;
  p.pdg = pdg;
  p.name = name;
  p.mass = mass;
  p.charge = charge;
  return p;
}

inline Particle composite(int pdg, const std::string& name, double mass, double ctau, std::vector<Particle> dau,
                          bool resonance = false) {
  Particle p;
  p.pdg = pdg;
  p.name = name;
  p.mass = mass;
  p.ctau = ctau;
  p.resonance = resonance;
  p.daughters = std::move(dau);
  return p;
}

/// B- -> (D0 -> K- pi+) K-
inline Particle BuToD0K() {
  const double mK = 0.493677, mPi = 0.13957, mD0 = 1.86484, mB = 5.27926;
  return composite(-521, "B-", mB, 0.491,
                   {composite(421, "D0", mD0, 0.1229, {leaf(-321, "K-", mK, -1), leaf(211, "pi+", mPi, +1)}),
                    leaf(-321, "K-", mK, -1)});
}

/// B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+
inline Particle BdToDstPi() {
  const double mK = 0.493677, mPi = 0.13957, mD0 = 1.86484, mDst = 2.01026, mB = 5.27958;
  auto D0 = composite(-421, "D~0", mD0, 0.1229, {leaf(321, "K+", mK, +1), leaf(-211, "pi-", mPi, -1)});
  auto Dst = composite(-413, "D*(2010)-", mDst, 0.0, {D0, leaf(-211, "pi-", mPi, -1)}, true);
  return composite(511, "B0", mB, 0.455, {Dst, leaf(211, "pi+", mPi, +1)});
}

/// B0 -> (K*0 -> K+ pi-) gamma
inline Particle BdToKstGamma() {
  const double mK = 0.493677, mPi = 0.13957, mKst = 0.89555, mB = 5.27958;
  auto Kst = composite(313, "K*0", mKst, 0.0, {leaf(321, "K+", mK, +1), leaf(-211, "pi-", mPi, -1)}, true);
  return composite(511, "B0", mB, 0.455, {Kst, leaf(22, "gamma", 0.0, 0)});
}

/// B0 -> (J/psi -> mu+ mu-) (KS0 -> pi+ pi-): a KS flying tens of cm
inline Particle BdToJpsiKS() {
  const double mMu = 0.1056584, mPi = 0.13957, mJpsi = 3.096916, mKS = 0.497614, mB = 5.27958;
  auto Jpsi = composite(443, "J/psi", mJpsi, 0.0, {leaf(-13, "mu+", mMu, +1), leaf(13, "mu-", mMu, -1)}, true);
  auto KS = composite(310, "KS0", mKS, 26.84, {leaf(211, "pi+", mPi, +1), leaf(-211, "pi-", mPi, -1)});
  return composite(511, "B0", mB, 0.455, {Jpsi, KS});
}

/// Largest transverse decay radius [mm] of any composite in the tree.
inline double maxDecayRadius(const Particle& p) {
  double r = p.daughters.empty() ? 0.0 : std::hypot(p.decay.x(), p.decay.y());
  for (const auto& d : p.daughters) r = std::max(r, maxDecayRadius(d));
  return r;
}

class Generator {
public:
  explicit Generator(unsigned seed, double Bz = 2.0) : m_rng(seed), m_Bz(Bz) {}

  double Bz() const { return m_Bz; }
  std::mt19937_64& rng() { return m_rng; }

  /// Decay `head` produced at `pv` with momentum |p| in [pmin, pmax], isotropic.
  void generate(Particle& head, const Eigen::Vector3d& pv, double pmin = 10, double pmax = 40) {
    std::uniform_real_distribution<double> u(0, 1);
    const double pm = pmin + (pmax - pmin) * u(m_rng);
    // keep away from the beam axis (|cos theta| < 0.8)
    const double ct = -0.8 + 1.6 * u(m_rng), phi = 2 * kPi * u(m_rng), st = std::sqrt(1 - ct * ct);
    const Eigen::Vector3d p(pm * st * std::cos(phi), pm * st * std::sin(phi), pm * ct);
    head.p4 << p, std::sqrt(p.squaredNorm() + head.mass * head.mass);
    decay(head, pv);
  }

  /// Smeared helix of a final-state particle with covariance `cov`.
  HelixMeasurement measure(const Particle& p, const helix::Matrix5d& cov) {
    HelixMeasurement h;
    h.Bz = m_Bz;
    h.cov = cov;
    h.par = helix::fromPositionMomentum(p.production, p.p4.head<3>(), p.charge, m_Bz);
    const Eigen::LLT<helix::Matrix5d> llt(cov);
    std::normal_distribution<double> n(0, 1);
    helix::Vector5d z;
    for (int i = 0; i < 5; ++i) z(i) = n(m_rng);
    h.par += llt.matrixL() * z;
    return h;
  }

  /// A Delphes-like helix covariance for a track of transverse momentum pt:
  /// 10 um impact parameters, correlated D0-phi and omega-phi.
  static helix::Matrix5d covariance(const Eigen::Vector3d& p, double Bz) {
    const double pt = std::hypot(p.x(), p.y());
    const double omega = helix::kC * Bz / pt;
    const double sD0 = 0.010, sZ0 = 0.012, sPhi = 1e-4, sOm = 2e-3 * omega, sTl = 1.5e-4;
    helix::Matrix5d c = helix::Matrix5d::Zero();
    c(0, 0) = sD0 * sD0;
    c(1, 1) = sPhi * sPhi;
    c(2, 2) = sOm * sOm;
    c(3, 3) = sZ0 * sZ0;
    c(4, 4) = sTl * sTl;
    c(0, 1) = c(1, 0) = -0.4 * sD0 * sPhi;
    c(1, 2) = c(2, 1) = 0.3 * sPhi * sOm;
    c(3, 4) = c(4, 3) = -0.3 * sZ0 * sTl;
    return c;
  }

  /// Calorimeter of the toy: a cylinder of radius kCaloR closed by endcaps
  /// at |z| = kCaloZ [mm], roughly the IDEA calorimeter front face.
  static constexpr double kCaloR = 2100, kCaloZ = 2500;

  /// Smeared cluster of a photon: the point where it reaches the calorimeter,
  /// with transverse resolution 3 mm, 20 mm along the shower (depth), and
  /// sigma(E)/E = 11%/sqrt(E) + 1% (dual-readout-like).
  PhotonMeasurement measurePhoton(const Particle& p) {
    const Eigen::Vector3d n = p.p4.head<3>().normalized();
    const Eigen::Vector3d x0 = p.production;
    // path length to the barrel and to the endcap; the smaller one is hit
    const double a = n.x() * n.x() + n.y() * n.y();
    const double b = x0.x() * n.x() + x0.y() * n.y();
    const double c = x0.x() * x0.x() + x0.y() * x0.y() - kCaloR * kCaloR;
    double t = a > 0 ? (-b + std::sqrt(b * b - a * c)) / a : 1e30;
    if (n.z() != 0) t = std::min(t, ((n.z() > 0 ? kCaloZ : -kCaloZ) - x0.z()) / n.z());
    const Eigen::Vector3d hit = x0 + t * n;

    const double E = p.p4(3);
    const double sT = 3.0, sD = 20.0, sE = E * std::hypot(0.11 / std::sqrt(E), 0.01);
    PhotonMeasurement m;
    m.cov.setZero();
    m.cov.topLeftCorner<3, 3>() =
        sT * sT * (Eigen::Matrix3d::Identity() - n * n.transpose()) + sD * sD * n * n.transpose();
    m.cov(3, 3) = sE * sE;
    const Eigen::LLT<Eigen::Matrix4d> llt(m.cov);
    std::normal_distribution<double> n01(0, 1);
    Eigen::Vector4d z;
    for (int i = 0; i < 4; ++i) z(i) = n01(m_rng);
    const Eigen::Vector4d smear = llt.matrixL() * z;
    m.clusterPos = hit + smear.head<3>();
    m.energy = E + smear(3);
    return m;
  }

  /// Build the fitter input; helix measurements for charged leaves,
  /// calorimeter clusters for photons.
  DecayNodeInput input(const Particle& p, const std::vector<int>& massConstrained = {}) {
    DecayNodeInput in;
    in.pdgCode = p.pdg;
    in.name = p.name;
    in.mass = p.mass;
    if (p.daughters.empty() && p.charge == 0) {
      in.measurement = measurePhoton(p);
      return in;
    }
    if (p.daughters.empty()) {
      TrackMeasurement t;
      t.charge = p.charge;
      t.helix = measure(p, covariance(p.p4.head<3>(), m_Bz));
      t.momentum = helix::momentum(t.helix->par, m_Bz);
      t.cov = helix::momentumCovariance(t.helix->par, t.helix->cov, m_Bz);
      in.measurement = t;
      return in;
    }
    in.isResonance = p.resonance;
    for (int pdg : massConstrained)
      if (pdg == p.pdg) in.applyMassConstraint = true;
    for (const auto& d : p.daughters) in.daughters.push_back(input(d, massConstrained));
    return in;
  }

private:
  std::mt19937_64 m_rng;
  double m_Bz;

  static Eigen::Vector4d boost(const Eigen::Vector4d& q, const Eigen::Vector4d& frame) {
    const double m = std::sqrt(frame(3) * frame(3) - frame.head<3>().squaredNorm());
    const Eigen::Vector3d b = frame.head<3>() / frame(3);
    const double b2 = b.squaredNorm(), g = frame(3) / m;
    const double bp = b.dot(q.head<3>());
    const double g2 = b2 > 0 ? (g - 1) / b2 : 0;
    Eigen::Vector4d out;
    out << q.head<3>() + (g2 * bp + g * q(3)) * b, g * (q(3) + bp);
    return out;
  }

  void decay(Particle& p, const Eigen::Vector3d& production) {
    p.production = production;
    if (p.daughters.empty()) return;
    std::exponential_distribution<double> e(1.0);
    p.trueCtau = p.ctau > 0 ? p.ctau * e(m_rng) : 0.0;
    const Eigen::Vector3d dir = p.p4.head<3>().normalized();
    const double L = p.trueCtau * p.p4.head<3>().norm() / p.mass;
    p.decay = production + L * dir;

    // sequential two- or more-body decay: flat phase space by acceptance-rejection
    // is overkill for tests; use isotropic two-body decays (all test modes are two-body)
    std::uniform_real_distribution<double> u(0, 1);
    auto& a = p.daughters[0];
    auto& b = p.daughters[1];
    const double M = p.mass;
    const double q = std::sqrt((M * M - std::pow(a.mass + b.mass, 2)) * (M * M - std::pow(a.mass - b.mass, 2))) / (2 * M);
    const double ct = 2 * u(m_rng) - 1, phi = 2 * kPi * u(m_rng), st = std::sqrt(1 - ct * ct);
    const Eigen::Vector3d qv(q * st * std::cos(phi), q * st * std::sin(phi), q * ct);
    Eigen::Vector4d qa, qb;
    qa << qv, std::sqrt(q * q + a.mass * a.mass);
    qb << -qv, std::sqrt(q * q + b.mass * b.mass);
    a.p4 = boost(qa, p.p4);
    b.p4 = boost(qb, p.p4);
    decay(a, p.decay);
    decay(b, p.decay);
  }
};

/// All particles of a tree in pre-order (matches FitResult::nodes).
inline void flatten(const Particle& p, std::vector<const Particle*>& out) {
  out.push_back(&p);
  for (const auto& d : p.daughters) flatten(d, out);
}

}  // namespace k4dtf::toy
