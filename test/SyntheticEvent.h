// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

/// A hand-made event in the EDM4hep POD layout of k4SimDelphes files, for
/// testing candidate building and truth matching without input files.
///
/// Columns are filled like the corresponding branches:
///   mc, mcParents, mcDaughters   <-> Particle, _Particle_parents/_daughters.index
///   recos, rpTracks              <-> ReconstructedParticles, _ReconstructedParticles_tracks.index
///   tracks, trackStates          <-> EFlowTrack, _EFlowTrack_trackStates
///   linkFrom, linkTo             <-> _RecoMCLink_from/_to.index
/// Charged reco particles get exact momenta and the exact helix of their MC
/// particle (DTFTrack::fromPositionMomentum) with a diagonal covariance, and
/// the pion mass/PDG like Delphes.

#include <cmath>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "k4DTF/DTFParticle.h"
#include "k4DTF/ParticleTable.h"

namespace k4DTF::test {

struct SyntheticEvent {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc;
  ROOT::VecOps::RVec<int> mcParents, mcDaughters;
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recos;
  ROOT::VecOps::RVec<int> rpTracks;
  ROOT::VecOps::RVec<edm4hep::TrackData> tracks;
  ROOT::VecOps::RVec<edm4hep::TrackState> trackStates;
  ROOT::VecOps::RVec<int> linkFrom, linkTo;
  double Bz{2.0};

  // relation bookkeeping while building
  std::vector<std::vector<int>> parentsOf, daughtersOf;

  /// Add an MC particle; returns its index.
  int addMC(int pdg, const ROOT::Math::PxPyPzEVector& p4, const ROOT::Math::XYZPoint& vertex, int parent = -1) {
    edm4hep::MCParticleData p{};
    p.PDG = pdg;
    p.charge = static_cast<float>(ParticleTable::charge(pdg));
    p.mass = p4.M();
    p.momentum = {p4.Px(), p4.Py(), p4.Pz()};
    p.vertex = {vertex.X(), vertex.Y(), vertex.Z()};
    p.generatorStatus = 1;
    mc.push_back(p);
    parentsOf.emplace_back();
    daughtersOf.emplace_back();
    const int idx = static_cast<int>(mc.size()) - 1;
    if (parent >= 0) {
      parentsOf[idx].push_back(parent);
      daughtersOf[parent].push_back(idx);
      mc[parent].generatorStatus = 2;
    }
    return idx;
  }

  /// Reconstruct MC particle m as a charged track (Delphes-like) or a neutral.
  int addReco(int m) {
    const auto& p = mc.at(m);
    edm4hep::ReconstructedParticleData rp{};
    rp.momentum = {static_cast<float>(p.momentum.x), static_cast<float>(p.momentum.y),
                   static_cast<float>(p.momentum.z)};
    rp.charge = p.charge;
    const double mass = p.charge != 0 ? ParticleTable::mass(211) : 0.0;
    const double pmag = std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y +
                                  p.momentum.z * p.momentum.z);
    rp.mass = static_cast<float>(mass);
    rp.energy = static_cast<float>(std::sqrt(pmag * pmag + mass * mass));
    rp.PDG = p.charge > 0 ? 211 : (p.charge < 0 ? -211 : 22);
    if (p.charge != 0) {
      auto trk = DTFTrack::fromPositionMomentum({p.vertex.x, p.vertex.y, p.vertex.z},
                                                {p.momentum.x, p.momentum.y, p.momentum.z}, p.charge, Bz);
      const double sigma[5] = {0.010, 1e-4, 1e-3 * std::abs(trk.par[DTFTrack::Omega]), 0.015, 1e-4};
      for (int i = 0; i < 5; ++i) trk.setCov(i, i, sigma[i] * sigma[i]);
      auto ts = trk.toEDM4hep();
      ts.location = edm4hep::TrackState::AtIP;
      edm4hep::TrackData td{};
      td.trackStates_begin = static_cast<unsigned>(trackStates.size());
      trackStates.push_back(ts);
      td.trackStates_end = static_cast<unsigned>(trackStates.size());
      rp.tracks_begin = static_cast<unsigned>(rpTracks.size());
      rpTracks.push_back(static_cast<int>(tracks.size()));
      rp.tracks_end = static_cast<unsigned>(rpTracks.size());
      tracks.push_back(td);
    } else {
      rp.tracks_begin = rp.tracks_end = static_cast<unsigned>(rpTracks.size());
    }
    recos.push_back(rp);
    linkFrom.push_back(static_cast<int>(recos.size()) - 1);
    linkTo.push_back(m);
    return static_cast<int>(recos.size()) - 1;
  }

  /// Fill the begin/end ranges and the flat relation columns.
  void finalise() {
    mcParents.clear();
    mcDaughters.clear();
    for (std::size_t i = 0; i < mc.size(); ++i) {
      mc[i].parents_begin = static_cast<unsigned>(mcParents.size());
      for (int q : parentsOf[i]) mcParents.push_back(q);
      mc[i].parents_end = static_cast<unsigned>(mcParents.size());
      mc[i].daughters_begin = static_cast<unsigned>(mcDaughters.size());
      for (int d : daughtersOf[i]) mcDaughters.push_back(d);
      mc[i].daughters_end = static_cast<unsigned>(mcDaughters.size());
    }
  }
};

// ---------------------------------------------------------------------------
// Kinematics helpers
// ---------------------------------------------------------------------------

/// Boost p (given in the rest frame of `frame`) to the frame where `frame` has its 4-momentum.
inline ROOT::Math::PxPyPzEVector boost(const ROOT::Math::PxPyPzEVector& p, const ROOT::Math::PxPyPzEVector& frame) {
  const double m = frame.M();
  const double bx = frame.Px() / frame.E(), by = frame.Py() / frame.E(), bz = frame.Pz() / frame.E();
  const double b2 = bx * bx + by * by + bz * bz;
  const double g = frame.E() / m;
  const double bp = bx * p.Px() + by * p.Py() + bz * p.Pz();
  const double g2 = b2 > 0 ? (g - 1.0) / b2 : 0.0;
  return {p.Px() + g2 * bp * bx + g * bx * p.E(), p.Py() + g2 * bp * by + g * by * p.E(),
          p.Pz() + g2 * bp * bz + g * bz * p.E(), g * (p.E() + bp)};
}

/// Two-body decay of `mother` into masses m1, m2 at polar/azimuthal angles in the rest frame.
inline std::pair<ROOT::Math::PxPyPzEVector, ROOT::Math::PxPyPzEVector> twoBody(
    const ROOT::Math::PxPyPzEVector& mother, double m1, double m2, double cosTheta, double phi) {
  const double M = mother.M();
  const double q = std::sqrt((M * M - (m1 + m2) * (m1 + m2)) * (M * M - (m1 - m2) * (m1 - m2))) / (2 * M);
  const double s = std::sqrt(1 - cosTheta * cosTheta);
  const double qx = q * s * std::cos(phi), qy = q * s * std::sin(phi), qz = q * cosTheta;
  const ROOT::Math::PxPyPzEVector a(qx, qy, qz, std::sqrt(q * q + m1 * m1));
  const ROOT::Math::PxPyPzEVector b(-qx, -qy, -qz, std::sqrt(q * q + m2 * m2));
  return {boost(a, mother), boost(b, mother)};
}

/// Decay point after a proper decay length ctau [mm].
inline ROOT::Math::XYZPoint flight(const ROOT::Math::XYZPoint& from, const ROOT::Math::PxPyPzEVector& p4,
                                   double ctau) {
  const double s = ctau / p4.M();  // L = ctau * p / m along p
  return {from.X() + s * p4.Px(), from.Y() + s * p4.Py(), from.Z() + s * p4.Pz()};
}

/// Indices of the MC particles of one B0 -> (D*- -> (D~0 -> K+ pi-) pi-) pi+ chain
/// (or its charge conjugate).
struct BdChain {
  int B, Dst, D0, K, piD0, piSoft, piB;
  int Bproduced;  ///< = B unless mixed
  ROOT::Math::XYZPoint pv, bVertex, dVertex;
};

/// Add B0 -> (D*- -> (D~0 -> K+ pi-) pi-) pi+ (cc = charge conjugate) to the event.
/// With mixed, the B is produced as its antiparticle and oscillates (MC record
/// B~0 -> B0 -> ...). ctauB and ctauD are the proper decay lengths [mm].
inline BdChain addBdChain(SyntheticEvent& ev, int zIndex, const ROOT::Math::XYZPoint& pv, bool cc, bool mixed,
                          double ctauB = 0.455, double ctauD = 0.123) {
  const int s = cc ? -1 : 1;
  const double mB = ParticleTable::mass(511), mDst = ParticleTable::mass(413), mD0 = ParticleTable::mass(421);
  const double mK = ParticleTable::mass(321), mPi = ParticleTable::mass(211);

  const double pB[3] = {12.0, -7.0, 25.0};
  const ROOT::Math::PxPyPzEVector b4(pB[0], pB[1], pB[2],
                                     std::sqrt(pB[0] * pB[0] + pB[1] * pB[1] + pB[2] * pB[2] + mB * mB));
  BdChain c{};
  c.pv = pv;
  int parent = zIndex;
  if (mixed) {
    c.Bproduced = ev.addMC(-s * 511, b4, pv, zIndex);
    parent = c.Bproduced;
  }
  c.B = ev.addMC(s * 511, b4, pv, parent);
  if (!mixed) c.Bproduced = c.B;

  c.bVertex = flight(pv, b4, ctauB);
  const auto [dst4, piB4] = twoBody(b4, mDst, mPi, 0.3, 1.1);
  c.Dst = ev.addMC(-s * 413, dst4, c.bVertex, c.B);
  c.piB = ev.addMC(s * 211, piB4, c.bVertex, c.B);

  const auto [d04, piS4] = twoBody(dst4, mD0, mPi, -0.5, 2.0);
  c.D0 = ev.addMC(-s * 421, d04, c.bVertex, c.Dst);  // strong decay: no flight
  c.piSoft = ev.addMC(-s * 211, piS4, c.bVertex, c.Dst);

  c.dVertex = flight(c.bVertex, d04, ctauD);
  const auto [k4, piD4] = twoBody(d04, mK, mPi, 0.7, -0.4);
  c.K = ev.addMC(s * 321, k4, c.dVertex, c.D0);
  c.piD0 = ev.addMC(-s * 211, piD4, c.dVertex, c.D0);
  return c;
}

/// e+ e- -> Z with the Z at the primary vertex; returns the Z index.
inline int addZ(SyntheticEvent& ev, const ROOT::Math::XYZPoint& pv) {
  const int em = ev.addMC(11, {0, 0, 45.6, 45.6}, pv);
  const int ep = ev.addMC(-11, {0, 0, -45.6, 45.6}, pv);
  const int z = ev.addMC(23, {0, 0, 0, 91.2}, pv, em);
  ev.parentsOf[z].push_back(ep);
  ev.daughtersOf[ep].push_back(z);
  return z;
}

}  // namespace k4DTF::test
