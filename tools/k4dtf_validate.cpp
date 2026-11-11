// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

/// k4dtf-validate: run k4DTF on the candidates of a stage-1 ntuple
/// (examples/dtf_test_samples/analysis/stage1_dtf_inputs.py) and write one
/// row per candidate with fit results, MC truth and pulls.
///
///   k4dtf-validate --input dtf_inputs_Bd2DstPi.root --output dtf_validation.root
///       --descriptor "B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+"
///       --labels B,Dst,D0,K,piD0,piSoft,piB [--mass-constraint 421] [--pv]
///       [--signal-only] [--max-events N]
///
/// The fitter is linked directly, so after changing the fitter a rebuild and
/// re-run is all that is needed; no EDM4hep input or FCCAnalyses. --pv uses
/// PV_reco as interaction-point constraint; --mass-constraint applies to the
/// listed PDG codes and their charge conjugates.
///
/// Output tree "candidates" (kInvalid = -999 where not available):
///   event, cand, nCand, isCC, isSignal
///   status, chi2, ndf, nIter, prob                    fit summary
///   pv_{true,reco,pull}_{x,y,z}                       primary vertex input vs MC
///   <node>_trueM, _trueP, _fitM, _fitMErr, _fitP, _pullP    every node
///   <node>_preM                                       composite: mass before the fit
///   <node>_trueVtx_{x,y,z}, _fitVtx_{x,y,z}, _fitVtxErr_{x,y,z}, _pullVtx_{x,y,z}
///   <node>_trueCtau, _fitCtau, _fitCtauErr, _pullCtau composite nodes
///   <leaf>_pull{D0,Phi,Omega,Z0,TanLambda}            input track vs MC helix (fit independent)

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/TreeInput.h"
#include "k4DTF/DTFParticle.h"

namespace {

constexpr float kInvalid = -999.f;
constexpr int kNoTrack = -1;  // TrackState::location of placeholder states (stage 1)

struct Options {
  std::string input, output = "dtf_validation.root", descriptor;
  std::vector<std::string> labels;
  std::vector<int> massConstraints;
  bool pv = false, signalOnly = false;
  long maxEvents = -1;
};

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  for (std::string item; std::getline(ss, item, sep);)
    if (!item.empty()) out.push_back(item);
  return out;
}

void usage() {
  std::cerr << "Usage: k4dtf-validate --input FILE --descriptor DESC --labels L0,L1,... [--output FILE]\n"
               "                      [--mass-constraint PDG[,PDG...]] [--pv] [--signal-only] [--max-events N]\n";
}

bool parse(int argc, char** argv, Options& o) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
      return argv[++i];
    };
    if (a == "--input") o.input = next();
    else if (a == "--output") o.output = next();
    else if (a == "--descriptor") o.descriptor = next();
    else if (a == "--labels") o.labels = split(next(), ',');
    else if (a == "--mass-constraint")
      for (const auto& p : split(next(), ',')) o.massConstraints.push_back(std::stoi(p));
    else if (a == "--pv") o.pv = true;
    else if (a == "--signal-only") o.signalOnly = true;
    else if (a == "--max-events") o.maxEvents = std::stol(next());
    else if (a == "-h" || a == "--help") return false;
    else throw std::invalid_argument("unknown option " + a);
  }
  return !o.input.empty() && !o.descriptor.empty() && !o.labels.empty();
}

double norm(const edm4hep::Vector3d& v) { return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

double wrapPi(double a) {
  constexpr double pi = 3.14159265358979323846;
  while (a > pi) a -= 2 * pi;
  while (a < -pi) a += 2 * pi;
  return a;
}

float pull(double fit, double truth, double err) {
  return (err > 0 && fit != kInvalid && truth != kInvalid) ? static_cast<float>((fit - truth) / err) : kInvalid;
}

/// Output branches: floats and ints addressed by name.
struct Row {
  std::map<std::string, float> f;
  std::map<std::string, int> i;
  Long64_t event = 0;

  void book(TTree& t) {
    t.Branch("event", &event);
    for (auto& [name, v] : i) t.Branch(name.c_str(), &v);
    for (auto& [name, v] : f) t.Branch(name.c_str(), &v);
  }
  void reset() {
    for (auto& [name, v] : f) v = kInvalid;
    for (auto& [name, v] : i) v = -1;
  }
};

}  // namespace

int main(int argc, char** argv) {
  Options opt;
  try {
    if (!parse(argc, argv, opt)) {
      usage();
      return 1;
    }
  } catch (const std::exception& e) {
    std::cerr << "k4dtf-validate: " << e.what() << "\n";
    usage();
    return 1;
  }

  // --- fitter configuration --------------------------------------------------
  const auto desc = k4DTF::DecayDescriptor::fromString(opt.descriptor);
  const auto ccDesc = desc.chargeConjugate();
  const auto& flat = desc.flatten();
  if (opt.labels.size() != flat.size()) {
    std::cerr << "k4dtf-validate: " << opt.labels.size() << " labels for " << flat.size() << " nodes of '"
              << desc.toString() << "'\n";
    return 1;
  }
  const k4DTF::DecayTreeFitter fitter;
  const auto leafNodes = desc.leafNodes();
  const auto nLeaves = leafNodes.size();

  // --- input -------------------------------------------------------------------
  std::unique_ptr<TFile> in(TFile::Open(opt.input.c_str()));
  if (!in || in->IsZombie()) {
    std::cerr << "k4dtf-validate: cannot open " << opt.input << "\n";
    return 1;
  }
  TTreeReader reader("events", in.get());
  TTreeReaderValue<float> Bz(reader, "Bz");
  TTreeReaderValue<edm4hep::VertexData> pvReco(reader, "PV_reco");
  std::unique_ptr<TTreeReaderValue<double>> pvTrue[3];
  for (int k = 0; k < 3; ++k)
    pvTrue[k] = std::make_unique<TTreeReaderValue<double>>(reader, (std::string("PV_true_") + "xyz"[k]).c_str());
  TTreeReaderArray<int> isCC(reader, "cand_isCC");
  TTreeReaderArray<int> isSignal(reader, "cand_isSignal");

  using RecoArr = TTreeReaderArray<edm4hep::ReconstructedParticleData>;
  using TrackArr = TTreeReaderArray<edm4hep::TrackState>;
  using MCArr = TTreeReaderArray<edm4hep::MCParticleData>;
  using DArr = TTreeReaderArray<double>;
  using FArr = TTreeReaderArray<float>;
  std::vector<std::unique_ptr<RecoArr>> leafReco;
  std::vector<std::unique_ptr<TrackArr>> leafTrack, leafTrue;
  std::vector<std::unique_ptr<MCArr>> nodeMC;
  std::map<int, std::unique_ptr<FArr>> preMass;
  std::map<int, std::unique_ptr<DArr>> trueCtau;
  std::map<int, std::array<std::unique_ptr<DArr>, 3>> trueVtx;
  const char* xyz[3] = {"x", "y", "z"};

  for (std::size_t n = 0; n < flat.size(); ++n) {
    const auto& L = opt.labels[n];
    nodeMC.push_back(std::make_unique<MCArr>(reader, (L + "_mc").c_str()));
    if (flat[n].isLeaf()) continue;
    preMass[n] = std::make_unique<FArr>(reader, (L + "_M").c_str());
    trueCtau[n] = std::make_unique<DArr>(reader, (L + "_trueCtau").c_str());
    for (int k = 0; k < 3; ++k)
      trueVtx[n][k] = std::make_unique<DArr>(reader, (L + "_trueDecayVtx_" + xyz[k]).c_str());
  }
  for (int n : leafNodes) {
    const auto& L = opt.labels[n];
    leafReco.push_back(std::make_unique<RecoArr>(reader, (L + "_reco").c_str()));
    leafTrack.push_back(std::make_unique<TrackArr>(reader, (L + "_track").c_str()));
    leafTrue.push_back(std::make_unique<TrackArr>(reader, (L + "_trueTrack").c_str()));
  }

  // --- output --------------------------------------------------------------------
  Row row;
  for (const char* name : {"cand", "nCand", "isCC", "isSignal", "status", "ndf", "nIter"}) row.i[name] = -1;
  for (const char* name : {"chi2", "prob"}) row.f[name] = kInvalid;
  for (const char* c : {"x", "y", "z"})
    for (const char* q : {"pv_true_", "pv_reco_", "pv_pull_"}) row.f[std::string(q) + c] = kInvalid;
  const char* helix[5] = {"D0", "Phi", "Omega", "Z0", "TanLambda"};
  for (std::size_t n = 0; n < flat.size(); ++n) {
    const auto& L = opt.labels[n];
    for (const char* q : {"_trueM", "_trueP", "_fitM", "_fitMErr", "_fitP", "_pullP"}) row.f[L + q] = kInvalid;
    if (flat[n].isLeaf()) {
      for (const char* h : helix) row.f[L + "_pull" + h] = kInvalid;
      continue;
    }
    for (const char* q : {"_preM", "_trueCtau", "_fitCtau", "_fitCtauErr", "_pullCtau"}) row.f[L + q] = kInvalid;
    for (const char* c : xyz)
      for (const char* q : {"_trueVtx_", "_fitVtx_", "_fitVtxErr_", "_pullVtx_"}) row.f[L + q + c] = kInvalid;
  }
  std::unique_ptr<TFile> out(TFile::Open(opt.output.c_str(), "RECREATE"));
  TTree tree("candidates", ("k4DTF validation: " + desc.toString()).c_str());
  row.book(tree);

  // --- loop ------------------------------------------------------------------------
  long nEvents = 0, nCands = 0, nSignal = 0;
  std::map<int, long> statusCount;
  const int vtxDiag[3] = {0, 2, 5};

  while (reader.Next()) {
    if (opt.maxEvents >= 0 && nEvents >= opt.maxEvents) break;
    ++nEvents;
    const auto nCand = isCC.GetSize();
    const auto pv = k4DTF::DTFVertex::fromEDM4hep(*pvReco);
    std::optional<k4DTF::InteractionPointInput> ip;
    if (opt.pv) ip = k4DTF::makeInteractionPoint(pv);

    for (std::size_t j = 0; j < nCand; ++j) {
      if (opt.signalOnly && !isSignal[j]) continue;
      row.reset();
      row.event = reader.GetCurrentEntry();
      row.i["cand"] = static_cast<int>(j);
      row.i["nCand"] = static_cast<int>(nCand);
      row.i["isCC"] = isCC[j];
      row.i["isSignal"] = isSignal[j];
      ++nCands;
      nSignal += isSignal[j];
      for (int k = 0; k < 3; ++k) {
        const double t = **pvTrue[k];
        const double rec = k == 0 ? pv.position.X() : (k == 1 ? pv.position.Y() : pv.position.Z());
        const double err = std::sqrt(std::max(0.0, pv.cov[vtxDiag[k]]));
        row.f[std::string("pv_true_") + xyz[k]] = static_cast<float>(t);
        row.f[std::string("pv_reco_") + xyz[k]] = static_cast<float>(rec);
        row.f[std::string("pv_pull_") + xyz[k]] = pull(rec, t, err);
      }

      // inputs: leaves with mass hypotheses (set in stage 1) and tracks
      std::vector<k4DTF::DTFParticle> leaves;
      for (std::size_t l = 0; l < nLeaves; ++l) {
        const auto& rp = (*leafReco[l])[j];
        const auto& ts = (*leafTrack[l])[j];
        auto p = ts.location != kNoTrack ? k4DTF::DTFParticle::fromEDM4hep(rp, ts, *Bz, rp.mass)
                                         : k4DTF::DTFParticle::fromEDM4hep(rp);
        leaves.push_back(std::move(p));

        // track-parameter pulls against the MC helix (independent of the fit)
        const auto& tt = (*leafTrue[l])[j];
        const auto& Lb = opt.labels[leafNodes[l]];
        if (ts.location != kNoTrack && tt.location != kNoTrack) {
          const double reco[5] = {ts.D0, ts.phi, ts.omega, ts.Z0, ts.tanLambda};
          const double tru[5] = {tt.D0, tt.phi, tt.omega, tt.Z0, tt.tanLambda};
          for (int k = 0; k < 5; ++k) {
            const double var = ts.covMatrix[k * (k + 1) / 2 + k];
            const double d = k == 1 ? wrapPi(reco[k] - tru[k]) : reco[k] - tru[k];
            row.f[Lb + "_pull" + helix[k]] = var > 0 ? static_cast<float>(d / std::sqrt(var)) : kInvalid;
          }
        }
      }

      k4DTF::FitResult r;
      try {
        r = fitter.fit(k4DTF::makeDecayNodeInput(isCC[j] ? ccDesc : desc, leaves, opt.massConstraints), ip);
      } catch (const std::exception&) {
        r.status = k4DTF::ErrCode::badsetup;
      }
      ++statusCount[static_cast<int>(r.status)];
      row.i["status"] = static_cast<int>(r.status);
      row.i["ndf"] = r.ndf;
      row.i["nIter"] = r.nIter;
      if (r.converged()) {
        row.f["chi2"] = static_cast<float>(r.chi2);
        row.f["prob"] = static_cast<float>(r.prob());
      }

      for (std::size_t n = 0; n < flat.size(); ++n) {
        const auto& L = opt.labels[n];
        const auto& mc = (*nodeMC[n])[j];
        const bool matched = mc.PDG != 0;
        if (matched) {
          row.f[L + "_trueM"] = static_cast<float>(mc.mass);
          row.f[L + "_trueP"] = static_cast<float>(norm(mc.momentum));
        }
        const bool fitted = r.converged() && n < r.nodes.size();
        if (fitted) {
          const auto& node = r.nodes[n];
          const Eigen::Vector3d p = node.p4.head<3>();
          row.f[L + "_fitM"] = static_cast<float>(node.mass);
          row.f[L + "_fitP"] = static_cast<float>(p.norm());
          if (node.massErr > 0) row.f[L + "_fitMErr"] = static_cast<float>(node.massErr);
          if (matched && p.norm() > 0) {
            const Eigen::Vector3d u = p / p.norm();
            const double pErr = std::sqrt(std::max(0.0, double(u.transpose() * node.p4Cov.topLeftCorner<3, 3>() * u)));
            row.f[L + "_pullP"] = pull(p.norm(), norm(mc.momentum), pErr);
          }
        }
        if (flat[n].isLeaf()) continue;

        row.f[L + "_preM"] = (*preMass[n])[j];
        const double tCtau = (*trueCtau[n])[j];
        row.f[L + "_trueCtau"] = static_cast<float>(tCtau);
        double tv[3];
        for (int k = 0; k < 3; ++k) {
          tv[k] = (*trueVtx[n][k])[j];
          row.f[L + "_trueVtx_" + xyz[k]] = static_cast<float>(tv[k]);
        }
        if (fitted && r.nodes[n].hasDecayVertex) {
          const auto& node = r.nodes[n];
          for (int k = 0; k < 3; ++k) {
            const double err = std::sqrt(std::max(0.0, node.decayVertexCov(k, k)));
            row.f[L + "_fitVtx_" + xyz[k]] = static_cast<float>(node.decayVertex(k));
            row.f[L + "_fitVtxErr_" + xyz[k]] = static_cast<float>(err);
            row.f[L + "_pullVtx_" + xyz[k]] = pull(node.decayVertex(k), tv[k], err);
          }
        }
        if (fitted && r.nodes[n].hasLifetime) {
          const auto& node = r.nodes[n];
          row.f[L + "_fitCtau"] = static_cast<float>(node.ctau);
          row.f[L + "_fitCtauErr"] = static_cast<float>(node.ctauErr);
          row.f[L + "_pullCtau"] = pull(node.ctau, tCtau, node.ctauErr);
        }
      }
      tree.Fill();
    }
  }

  out->cd();
  tree.Write();
  out->Close();

  std::cout << "k4dtf-validate: " << nEvents << " events, " << nCands << " candidates (" << nSignal
            << " truth-matched) -> " << opt.output << "\n  fit status:";
  for (const auto& [s, c] : statusCount) std::cout << " [" << static_cast<k4DTF::ErrCode>(s) << "] " << c;
  std::cout << "\n";
  return 0;
}
