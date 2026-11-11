// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DTFParticle.h"

#include <algorithm>
#include <cmath>

#include "k4DTF/core/Helix.h"

namespace k4DTF {


namespace {

/// Index into a lower-triangular row-major array.
constexpr int lowerTri(int i, int j) {
  if (i < j) std::swap(i, j);
  return i * (i + 1) / 2 + j;
}

}  // namespace

// ---------------------------------------------------------------------------
// DTFTrack
// ---------------------------------------------------------------------------

double DTFTrack::covAt(int i, int j) const { return cov[lowerTri(i, j)]; }
void DTFTrack::setCov(int i, int j, double v) { cov[lowerTri(i, j)] = v; }

namespace {
helix::Vector5d toEigen(const std::array<double, 5>& par) {
  helix::Vector5d v;
  for (int i = 0; i < 5; ++i) v(i) = par[i];
  return v;
}
}  // namespace

double DTFTrack::charge() const { return helix::charge(toEigen(par), Bz); }

ROOT::Math::XYZVector DTFTrack::momentum() const {
  const auto p = helix::momentum(toEigen(par), Bz);
  return {p.x(), p.y(), p.z()};
}

ROOT::Math::XYZPoint DTFTrack::pca() const {
  const auto x = helix::pca(toEigen(par), {referencePoint.X(), referencePoint.Y(), referencePoint.Z()});
  return {x.x(), x.y(), x.z()};
}

std::array<double, 10> DTFTrack::p4Covariance(double mass) const {
  std::array<double, 10> out{};
  if (par[Omega] == 0.0) return out;
  const auto pv = toEigen(par);
  const Eigen::Vector3d p = helix::momentum(pv, Bz);
  const double E = std::sqrt(p.squaredNorm() + mass * mass);
  Eigen::Matrix<double, 4, 5> J;
  J.topRows<3>() = helix::momentumJacobian(pv, Bz);
  J.row(3) = p.transpose() * J.topRows<3>() / E;
  helix::Matrix5d C;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) C(i, j) = covAt(i, j);
  const Eigen::Matrix4d c4 = J * C * J.transpose();
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b <= a; ++b) out[lowerTri(a, b)] = c4(a, b);
  return out;
}

DTFTrack DTFTrack::fromEDM4hep(const edm4hep::TrackState& ts, double Bz) {
  DTFTrack t;
  t.par = {ts.D0, ts.phi, ts.omega, ts.Z0, ts.tanLambda};
  for (int i = 0; i < 15; ++i) t.cov[i] = ts.covMatrix[i];
  t.referencePoint = {ts.referencePoint.x, ts.referencePoint.y, ts.referencePoint.z};
  t.Bz = Bz;
  return t;
}

edm4hep::TrackState DTFTrack::toEDM4hep() const {
  edm4hep::TrackState ts{};
  ts.location = edm4hep::TrackState::AtOther;
  ts.D0 = static_cast<float>(par[D0]);
  ts.phi = static_cast<float>(par[Phi]);
  ts.omega = static_cast<float>(par[Omega]);
  ts.Z0 = static_cast<float>(par[Z0]);
  ts.tanLambda = static_cast<float>(par[TanLambda]);
  ts.referencePoint = {static_cast<float>(referencePoint.X()), static_cast<float>(referencePoint.Y()),
                       static_cast<float>(referencePoint.Z())};
  for (int i = 0; i < 15; ++i) ts.covMatrix[i] = static_cast<float>(cov[i]);
  return ts;
}

k4DTF::HelixMeasurement DTFTrack::toHelixMeasurement() const {
  k4DTF::HelixMeasurement h;
  h.par = toEigen(par);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) h.cov(i, j) = covAt(i, j);
  h.referencePoint = {referencePoint.X(), referencePoint.Y(), referencePoint.Z()};
  h.Bz = Bz;
  return h;
}

DTFTrack DTFTrack::fromPositionMomentum(const ROOT::Math::XYZPoint& x, const ROOT::Math::XYZVector& p,
                                        double charge, double Bz) {
  DTFTrack t;
  t.Bz = Bz;
  const auto h = helix::fromPositionMomentum({x.X(), x.Y(), x.Z()}, {p.X(), p.Y(), p.Z()}, charge, Bz);
  for (int i = 0; i < 5; ++i) t.par[i] = h(i);
  return t;
}

// ---------------------------------------------------------------------------
// DTFVertex
// ---------------------------------------------------------------------------

DTFVertex DTFVertex::fromEDM4hep(const edm4hep::VertexData& v) {
  DTFVertex out;
  out.position = {v.position.x, v.position.y, v.position.z};
  for (int i = 0; i < 6; ++i) out.cov[i] = v.covMatrix[i];
  return out;
}

// ---------------------------------------------------------------------------
// DTFParticle
// ---------------------------------------------------------------------------

void DTFParticle::setMassHypothesis(double m) {
  const double p2 = p4.Px() * p4.Px() + p4.Py() * p4.Py() + p4.Pz() * p4.Pz();
  p4.SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), std::sqrt(p2 + m * m));
  if (track) cov = track->p4Covariance(m);
}

DTFParticle DTFParticle::fromEDM4hep(const edm4hep::ReconstructedParticleData& rp) {
  DTFParticle dp;
  dp.p4.SetPxPyPzE(rp.momentum.x, rp.momentum.y, rp.momentum.z, rp.energy);
  dp.charge = rp.charge;
  dp.pdg    = rp.PDG;
  for (int i = 0; i < 10; ++i) dp.cov[i] = rp.covMatrix[i];
  return dp;
}

DTFParticle DTFParticle::fromEDM4hep(const edm4hep::ReconstructedParticleData& rp, const edm4hep::TrackState& ts,
                                     double Bz, double massHypothesis) {
  DTFParticle dp = fromEDM4hep(rp);
  dp.track = DTFTrack::fromEDM4hep(ts, Bz);
  const double m = massHypothesis > 0 ? massHypothesis : rp.mass;
  if (massHypothesis > 0) dp.setMassHypothesis(massHypothesis);
  const bool noCov = std::all_of(dp.cov.begin(), dp.cov.end(), [](double c) { return c == 0.0; });
  if (noCov) dp.cov = dp.track->p4Covariance(m);
  return dp;
}

edm4hep::ReconstructedParticleData DTFParticle::toEDM4hep() const {
  edm4hep::ReconstructedParticleData rp{};
  rp.momentum  = {static_cast<float>(p4.Px()),
                  static_cast<float>(p4.Py()),
                  static_cast<float>(p4.Pz())};
  rp.energy    = static_cast<float>(p4.E());
  rp.mass      = static_cast<float>(mass());
  rp.charge    = static_cast<float>(charge);
  rp.PDG       = pdg;
  for (int i = 0; i < 10; ++i) rp.covMatrix[i] = static_cast<float>(cov[i]);
  return rp;
}

}  // namespace k4DTF
