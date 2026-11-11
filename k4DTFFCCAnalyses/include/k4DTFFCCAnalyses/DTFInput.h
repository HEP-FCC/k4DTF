// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

/// Building decay-tree candidates and their MC truth in FCCAnalyses / RDataFrame.
///
/// FCCAnalyses reads EDM4hep files at the POD level, where relations are
/// index columns, so everything here works on indices:
///
///   * a candidate is the list of ReconstructedParticles indices of its final-
///     state particles, in the depth-first leaf order of the decay descriptor
///     (the order DecayTreeFitter::fit() expects);
///   * nodes of the tree are numbered as in DecayDescriptor::flatten()
///     (pre-order, node 0 = head).
///
/// Column names for EDM4hep v1 files from k4SimDelphes (check with podio-dump):
///   ReconstructedParticles                     RVec<ReconstructedParticleData>
///   _ReconstructedParticles_tracks.index       reco -> EFlowTrack
///   _ReconstructedParticles_clusters.index     reco -> EFlowPhoton / EFlowNeutralHadron
///   EFlowPhoton                                RVec<ClusterData> (photon clusters)
///   EFlowTrack, _EFlowTrack_trackStates        RVec<TrackData>, RVec<TrackState>
///   _RecoMCLink_from.index, _RecoMCLink_to.index   reco <-> MC links
///   Particle, _Particle_parents.index, _Particle_daughters.index
///   magFieldBz                                 RVec<float> (one entry)

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ClusterData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"
#include "k4DTF/CalorimeterResolution.h"
#include "k4DTF/DecayDescriptor.h"

namespace k4DTF {

using Vec_rp = ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>;
using Vec_mc = ROOT::VecOps::RVec<edm4hep::MCParticleData>;
using Vec_ts = ROOT::VecOps::RVec<edm4hep::TrackState>;
using Vec_trk = ROOT::VecOps::RVec<edm4hep::TrackData>;
using Vec_cl = ROOT::VecOps::RVec<edm4hep::ClusterData>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_cand = ROOT::VecOps::RVec<Vec_i>;

/// Value written for quantities that are not available (unmatched, not fitted).
constexpr float kInvalid = -999.f;

/// TrackState::location given to placeholder states of particles without a track.
constexpr int kNoTrack = -1;

/// ClusterData::type given to placeholder clusters of particles that are not photons.
constexpr int kNoCluster = -1;

// ---------------------------------------------------------------------------
// Per reconstructed particle
// ---------------------------------------------------------------------------

/// MC index for every reconstructed particle from the RecoMCLink columns
/// (-1 if not linked; the first link is used if there are several).
Vec_i recoToMC(const Vec_rp& recos, const Vec_i& linkFrom, const Vec_i& linkTo);

/// PDG code of the linked MC particle for every reconstructed particle
/// (0 if not linked). Used as perfect PID in the test samples.
Vec_i truthPDG(const Vec_i& recoMC, const Vec_mc& mc);

/// Track state for every reconstructed particle, aligned with
/// ReconstructedParticles (unlike FCCAnalyses' getRP2TRK, which skips
/// neutrals). Particles without a track get a zero state with
/// location = kNoTrack. Uses the first track state of the first track.
Vec_ts trackStatePerReco(const Vec_rp& recos, const Vec_i& rpTrackIdx, const Vec_trk& tracks,
                         const Vec_ts& trackStates);

/// Photon cluster for every reconstructed particle, aligned with
/// ReconstructedParticles, ready for the fit: position, energy and their
/// errors. Particles that are not photons get a placeholder with
/// type = kNoCluster and zero energy.
///
/// A neutral particle is a photon if its first cluster is in photonClusters
/// (EFlowPhoton for Delphes). The cluster relation holds an index without
/// the collection, and neutral hadrons point into EFlowNeutralHadron, so the
/// cluster is accepted only if its energy equals the particle's (Delphes
/// writes the same tower energy to both). Neutral particles with PDG 22 are
/// photons without a matching cluster too.
///
/// What is missing is filled in, so that files from every k4SimDelphes
/// version work (before v00-09, clusters have neither a position nor errors):
///   - position: the cluster's; else the particle's reference point (set to
///     the cluster position by v00-09); else where the particle's direction
///     meets the calorimeter front face (calo.frontFace());
///   - position covariance: the cluster's; else calo.positionCovariance();
///   - energy: the particle's; its error: the cluster's, else from calo.
/// With modelErrors the errors always come from calo, also for files that
/// have them: k4SimDelphes v00-09 writes tower size / sqrt(12), while the
/// error of a position drawn uniformly within the tower (Delphes
/// SmearTowerCenter, as in the IDEA card) is tower size / sqrt(6).
Vec_cl photonClusterPerReco(const Vec_rp& recos, const Vec_i& rpClusterIdx, const Vec_cl& photonClusters,
                            const CalorimeterResolution& calo = CalorimeterResolution(), bool modelErrors = false);

/// First element of a one-entry column such as magFieldBz, or def if empty.
float firstOr(const Vec_f& v, float def);

/// pid with 0 for every reconstructed particle that does not come from a
/// decay of a particle with |PDG| = |headPdg| (its MC particle has no such
/// ancestor, or it is not linked). Passed to CandidateBuilder, this keeps
/// only particles that can be part of a signal candidate, which avoids most
/// of the combinatorics in busy events; select the truth-matched candidates
/// afterwards with selectCandidates(..., TruthMatch::isSignal).
Vec_i signalPID(const Vec_i& pid, const Vec_i& recoMC, const Vec_mc& mc, const Vec_i& mcParents, int headPdg);

// ---------------------------------------------------------------------------
// Candidates
// ---------------------------------------------------------------------------

struct Candidates {
  Vec_cand leaves;  ///< per candidate: reco indices of the leaves, descriptor leaf order
  Vec_i isCC;       ///< per candidate: 1 if built for the charge-conjugate descriptor
};

/// Builds all candidates for a decay descriptor, node by node.
///
/// Leaves are taken from the reconstructed particles whose PID hypothesis
/// (pid, one entry per reco particle) equals the leaf's PDG code and whose
/// charge matches. Composite nodes combine their daughters' candidates
/// without reusing a particle; identical sibling subtrees (e.g. the two pi-
/// in D- -> K+ pi- pi-) are combined once. If a mass window is given for a
/// node's |PDG|, its invariant mass (leaf masses from the descriptor) must lie
/// inside. With includeCC the charge-conjugate descriptor is built as well,
/// unless it has the same final state.
class CandidateBuilder {
public:
  using MassWindows = std::map<int, std::pair<double, double>>;  ///< |PDG| -> (min, max) [GeV]

  CandidateBuilder(const std::string& descriptor, MassWindows massWindows = {}, bool includeCC = true,
                   std::size_t maxCandidates = 10000);

  Candidates operator()(const Vec_rp& recos, const Vec_i& pid) const;

  const DecayDescriptor& descriptor() const { return m_desc; }
  const DecayDescriptor& ccDescriptor() const { return m_cc; }

private:
  DecayDescriptor m_desc;
  DecayDescriptor m_cc;
  MassWindows m_windows;
  bool m_includeCC;
  std::size_t m_max;

  Vec_cand build(const DecayDescriptor& d, const Vec_rp& recos, const Vec_i& pid) const;
};

/// Number of candidates.
int nCandidates(const Candidates& c);

/// Candidates j with mask[j] != 0, e.g. mask = TruthMatch::isSignal.
Candidates selectCandidates(const Candidates& c, const Vec_i& mask);

/// Reco index of leaf `leaf` for every candidate.
Vec_i leafRecoIndex(const Candidates& c, int leaf);

/// Leaf `leaf` of every candidate with the descriptor's mass hypothesis and
/// PDG code (charge-conjugated for CC candidates); the 4-momentum covariance
/// is propagated from the track if the reconstructed particle has none.
Vec_rp leafParticles(const Candidates& c, const Vec_rp& recos, const Vec_ts& tsPerReco,
                     const CandidateBuilder& b, int leaf, float Bz);

/// Track state of leaf `leaf` for every candidate.
Vec_ts leafTrackStates(const Candidates& c, const Vec_ts& tsPerReco, int leaf);

/// Cluster (from photonClusterPerReco) of leaf `leaf` for every candidate.
Vec_cl leafClusters(const Candidates& c, const Vec_cl& clusterPerReco, int leaf);

/// Invariant mass [GeV] of flat node `node` for every candidate, with the
/// descriptor's leaf mass hypotheses.
Vec_f nodeMass(const Candidates& c, const Vec_rp& recos, const CandidateBuilder& b, int node);

/// Number of the candidate's leaf tracks that were used in a primary-vertex
/// fit (pvTracks, e.g. VertexFitterSimple::get_PrimaryTracks).
Vec_i nLeavesInPV(const Candidates& c, const Vec_ts& tsPerReco, const Vec_ts& pvTracks);

// ---------------------------------------------------------------------------
// MC truth
// ---------------------------------------------------------------------------

struct TruthMatch {
  Vec_cand nodeMC;     ///< per candidate, per flat node: MC index, -1 if not matched
  Vec_i isSignal;      ///< 1 if every node is matched to the right particle
  Vec_i nMatched;      ///< number of matched nodes
  Vec_i headOrigin;    ///< MC index of the head before B0/Bs mixing (= head if unmixed), -1 if unmatched
};

/// Matches candidates to the MC decay chain.
///
/// A leaf is matched if its reco particle links to an MC particle with the
/// leaf's PDG code. A composite node is matched if all its daughters are
/// matched, share the same MC mother with the node's PDG code, and that
/// mother has no other daughters (photons excepted, for PHOTOS). For the
/// head, a mixing chain (B0 -> anti-B0 -> ...) is followed back to the
/// produced particle, stored in headOrigin.
class TruthMatcher {
public:
  explicit TruthMatcher(const CandidateBuilder& builder);

  TruthMatch operator()(const Candidates& c, const Vec_i& recoMC, const Vec_mc& mc, const Vec_i& mcParents,
                        const Vec_i& mcDaughters) const;

private:
  DecayDescriptor m_desc;
  DecayDescriptor m_cc;
};

/// Truth of the candidates j with mask[j] != 0 (same order as
/// selectCandidates(c, mask)).
TruthMatch selectCandidates(const TruthMatch& t, const Vec_i& mask);

/// Matched MC particle of node `node` for every candidate (default-constructed,
/// PDG 0, if not matched). For the head this is the decaying particle.
Vec_mc nodeMCParticle(const TruthMatch& t, const Vec_mc& mc, int node);

/// True production vertex coordinate (0 = x, 1 = y, 2 = z) [mm] of node `node`.
/// For the head the vertex of headOrigin is used.
Vec_d nodeTrueProductionVertex(const TruthMatch& t, const Vec_mc& mc, int node, int coord);

/// True decay vertex coordinate [mm] of composite node `node`: production
/// vertex of its first MC daughter (k4SimDelphes does not fill endpoint).
Vec_d nodeTrueDecayVertex(const TruthMatch& t, const Vec_mc& mc, const Vec_i& mcDaughters, int node, int coord);

/// True proper decay length c*tau [mm] of composite node `node`:
/// |decay - production vertex| * m / p.
Vec_d nodeTrueCtau(const TruthMatch& t, const Vec_mc& mc, const Vec_i& mcDaughters, int node);

/// Exact helix parameters of the MC particle matched to leaf `leaf`
/// (DTFTrack::fromPositionMomentum), as an edm4hep::TrackState with zero
/// covariance; location = kNoTrack if not matched. Compare with
/// leafTrackStates() for track-parameter pulls.
Vec_ts leafTrueTrackState(const TruthMatch& t, const Vec_mc& mc, const CandidateBuilder& b, int leaf, float Bz);

/// True primary vertex coordinate [mm]: production vertex of the Z (PDG 23),
/// or of the first MC particle without parents.
double truePrimaryVertex(const Vec_mc& mc, int coord);

}  // namespace k4DTF
