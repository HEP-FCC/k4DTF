// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <optional>
#include <string>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "edm4hep/ClusterData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"
#include "k4DTF/core/ParticleMeasurement.h"

namespace k4DTF {

/// Helix parameters of a charged track at a reference point, with ROOT types
/// for analysis code. The conventions and formulae are those of
/// k4DTF::helix (k4DTF/core/Helix.h), used by the fitter.
///
/// Parameters and units follow edm4hep::TrackState as written by k4SimDelphes
/// (Delphes IDEA card):
///   par = (D0 [mm], phi [rad], omega [1/mm], Z0 [mm], tanLambda)
/// with
///   - phi       azimuth of the momentum at the point of closest approach (PCA),
///   - omega     signed curvature, omega = -2 C where C is Delphes' half-curvature,
///               so sign(omega) = sign(charge * Bz),
///   - D0, Z0    such that the PCA is reference + (-D0 sin(phi), D0 cos(phi), Z0)
///               (Delphes TrkUtil::ParToX convention),
///   - reference point (0, 0, 0) for Delphes tracks.
/// cov is the lower triangle of the 5x5 covariance in EDM4hep order
/// (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) ...; the 6th (time) row of the EDM4hep
/// matrix is not used (Delphes leaves it at zero).
struct DTFTrack {
  enum Param { D0 = 0, Phi = 1, Omega = 2, Z0 = 3, TanLambda = 4 };

  std::array<double, 5> par{};
  std::array<double, 15> cov{};
  ROOT::Math::XYZPoint referencePoint{};  // [mm]
  double Bz{0.0};                         // [T]

  /// Covariance element (i, j), symmetric access.
  double covAt(int i, int j) const;
  void setCov(int i, int j, double v);

  /// Charge from the sign of omega and Bz.
  double charge() const;

  /// Momentum [GeV] at the PCA.
  ROOT::Math::XYZVector momentum() const;

  /// Point of closest approach [mm].
  ROOT::Math::XYZPoint pca() const;

  /// (px, py, pz, E) covariance at the PCA for a mass hypothesis, propagated
  /// from (phi, omega, tanLambda); 10-element lower triangle as in DTFParticle::cov.
  std::array<double, 10> p4Covariance(double mass) const;

  static DTFTrack fromEDM4hep(const edm4hep::TrackState& ts, double Bz);
  edm4hep::TrackState toEDM4hep() const;

  /// Measurement for the fitter (k4DTF::TrackMeasurement::helix).
  k4DTF::HelixMeasurement toHelixMeasurement() const;

  /// Exact helix of a particle produced at x [mm] with momentum p [GeV] and
  /// the given charge, reference point (0, 0, 0). Same formulae as Delphes
  /// TrkUtil::XPtoPar, so a perfectly measured Delphes track reproduces it.
  /// Use it to compute truth track parameters from MC particles.
  static DTFTrack fromPositionMomentum(const ROOT::Math::XYZPoint& x, const ROOT::Math::XYZVector& p,
                                       double charge, double Bz);
};

/// A photon's calorimeter cluster: position [mm] with covariance, energy
/// [GeV] with its error, as in edm4hep::Cluster. The fitter assumes the
/// photon flies from its production vertex to the cluster position.
struct DTFCluster {
  ROOT::Math::XYZPoint position{};    // [mm]
  std::array<double, 6> positionCov{};  // lower triangle (xx, yx, yy, zx, zy, zz) [mm^2]
  double energy{0.0};                 // [GeV]
  double energyError{0.0};            // [GeV]

  /// Momentum of a massless particle from the origin to the cluster.
  ROOT::Math::XYZVector momentum() const;

  /// (px, py, pz, E) covariance of momentum(); 10-element lower triangle.
  std::array<double, 10> p4Covariance() const;

  static DTFCluster fromEDM4hep(const edm4hep::ClusterData& c);
  edm4hep::ClusterData toEDM4hep() const;

  /// Measurement for the fitter.
  k4DTF::PhotonMeasurement toPhotonMeasurement() const;
};

/// A vertex with position covariance, e.g. the primary vertex.
struct DTFVertex {
  ROOT::Math::XYZPoint position{};  // [mm]
  std::array<double, 6> cov{};      // lower triangle (xx, yx, yy, zx, zy, zz) [mm^2]

  static DTFVertex fromEDM4hep(const edm4hep::VertexData& v);
};

/// Internal particle representation used by the fitter.
///
/// Covariance matrix convention: lower triangle of the 4×4 (px,py,pz,E)
/// covariance stored row-major → 10 elements, matching the layout of
/// edm4hep::ReconstructedParticleData::covMatrix.
///
/// Charged final-state particles should also carry their track (helix
/// parameters and covariance); a vertex-constrained decay-tree fit needs it,
/// and Delphes reconstructed particles have an empty 4-momentum covariance.
struct DTFParticle {
  ROOT::Math::PxPyPzEVector p4;
  std::array<double, 10> cov{};  // (px,py,pz,E) covariance
  double charge{0.0};
  int pdg{0};
  std::string name;
  std::optional<DTFTrack> track;  ///< set for charged final-state particles
  std::optional<DTFCluster> cluster;  ///< set for photons

  /// Invariant mass from the 4-momentum.
  double mass() const { return p4.M(); }

  bool hasTrack() const { return track.has_value(); }
  bool hasCluster() const { return cluster.has_value(); }

  /// Keep the 3-momentum, set E = sqrt(p^2 + m^2). If a track is attached,
  /// cov is recomputed for the new mass.
  void setMassHypothesis(double m);

  /// Construct from an EDM4hep reconstructed particle data struct.
  static DTFParticle fromEDM4hep(const edm4hep::ReconstructedParticleData& rp);

  /// Construct a charged final-state particle from its reconstructed particle
  /// and track state. With massHypothesis > 0 the energy is recomputed for
  /// that mass. If rp has no 4-momentum covariance (Delphes), cov is
  /// propagated from the track covariance.
  static DTFParticle fromEDM4hep(const edm4hep::ReconstructedParticleData& rp, const edm4hep::TrackState& ts,
                                 double Bz, double massHypothesis = -1.0);

  /// Construct a photon from its reconstructed particle and cluster (see
  /// photonClusterPerReco() in k4DTFFCCAnalyses for Delphes files). The
  /// 4-momentum points from the origin to the cluster, with the cluster
  /// energy; its covariance is propagated from the cluster if rp has none.
  static DTFParticle fromEDM4hep(const edm4hep::ReconstructedParticleData& rp, const edm4hep::ClusterData& cl);

  /// Convert back to an EDM4hep reconstructed particle data struct.
  edm4hep::ReconstructedParticleData toEDM4hep() const;
};

}  // namespace k4DTF
