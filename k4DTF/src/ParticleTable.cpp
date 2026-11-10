// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/ParticleTable.h"

#include <cstdlib>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

namespace k4DTF::ParticleTable {

namespace {

struct Entry {
  int pdg;
  double mass;    // GeV
  double charge;  // e, of the particle with positive PDG code
  bool selfConjugate;
  std::vector<std::string> names;  // first one is canonical
};

// Masses from EvtGen evt.pdl (FCC-config).
const std::vector<Entry>& entries() {
  static const std::vector<Entry> table = {
      {11, 0.0005110, -1, false, {"e-"}},
      {13, 0.1056584, -1, false, {"mu-"}},
      {15, 1.7768200, -1, false, {"tau-"}},
      {22, 0.0, 0, true, {"gamma"}},
      {111, 0.1349766, 0, true, {"pi0"}},
      {211, 0.1395702, +1, false, {"pi+"}},
      {321, 0.4936770, +1, false, {"K+"}},
      {310, 0.4976140, 0, true, {"KS0", "K_S0", "KS", "Ks0"}},
      {130, 0.4976140, 0, true, {"KL0", "K_L0", "KL"}},
      {311, 0.4976140, 0, false, {"K0"}},
      {313, 0.8958100, 0, false, {"K*(892)0", "K*0"}},
      {323, 0.8916600, +1, false, {"K*(892)+", "K*+"}},
      {113, 0.7752600, 0, true, {"rho(770)0", "rho0"}},
      {213, 0.7752600, +1, false, {"rho(770)+", "rho+"}},
      {221, 0.5478620, 0, true, {"eta"}},
      {331, 0.9577800, 0, true, {"eta_prime", "eta'"}},
      {223, 0.7826500, 0, true, {"omega(782)", "omega"}},
      {333, 1.0194610, 0, true, {"phi(1020)", "phi"}},
      {443, 3.0969160, 0, true, {"J/psi(1S)", "J/psi", "Jpsi"}},
      {100443, 3.6861090, 0, true, {"psi(2S)"}},
      {2212, 0.9382721, +1, false, {"p+", "p", "proton"}},
      {2112, 0.9395654, 0, false, {"n0", "n"}},
      {3122, 1.1156830, 0, false, {"Lambda0", "Lambda", "L0"}},
      {411, 1.8696100, +1, false, {"D+"}},
      {421, 1.8648400, 0, false, {"D0"}},
      {431, 1.9683000, +1, false, {"D_s+", "Ds+"}},
      {413, 2.0102600, +1, false, {"D*(2010)+", "D*+"}},
      {423, 2.0069600, 0, false, {"D*(2007)0", "D*0"}},
      {433, 2.1121000, +1, false, {"D*_s+", "D_s*+", "Ds*+"}},
      {4122, 2.2864600, +1, false, {"Lambda_c+", "Lc+"}},
      {4232, 2.4678000, +1, false, {"Xi_c+"}},
      {4132, 2.4708800, 0, false, {"Xi_c0"}},
      {521, 5.2792600, +1, false, {"B+"}},
      {511, 5.2795800, 0, false, {"B0", "Bd0", "B_d0"}},
      {531, 5.3667700, 0, false, {"B_s0", "Bs0"}},
      {541, 6.2756000, +1, false, {"B_c+", "Bc+"}},
      {5122, 5.6195000, 0, false, {"Lambda_b0", "Lb0", "Lambda_b"}},
  };
  return table;
}

char lastChar(const std::string& s) { return s.empty() ? '\0' : s.back(); }
bool hasChargeSuffix(const std::string& s) { return lastChar(s) == '+' || lastChar(s) == '-'; }

std::string flipChargeSuffix(std::string s) {
  if (lastChar(s) == '+')
    s.back() = '-';
  else if (lastChar(s) == '-')
    s.back() = '+';
  return s;
}

const std::unordered_map<std::string, int>& nameMap() {
  static const auto map = [] {
    std::unordered_map<std::string, int> m;
    for (const auto& e : entries()) {
      for (const auto& n : e.names) {
        m.emplace(n, e.pdg);
        // Charged particles: the opposite-charge spelling is the antiparticle.
        if (hasChargeSuffix(n)) m.emplace(flipChargeSuffix(n), -e.pdg);
      }
    }
    return m;
  }();
  return map;
}

const std::unordered_map<int, const Entry*>& pdgMap() {
  static const auto map = [] {
    std::unordered_map<int, const Entry*> m;
    for (const auto& e : entries()) m.emplace(e.pdg, &e);
    return m;
  }();
  return map;
}

bool endsWith(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Resolve antiparticle spellings by stripping the marker and conjugating.
// For charged names the charge suffix already states the antiparticle
// ("p~-", "anti-Lambda_c-"), so the stripped name resolves directly.
bool resolve(const std::string& name, int& pdg) {
  if (const auto it = nameMap().find(name); it != nameMap().end()) {
    pdg = it->second;
    return true;
  }
  auto conjugated = [&](const std::string& rest) {
    int p = 0;
    if (rest.empty() || !resolve(rest, p)) return false;
    pdg = hasChargeSuffix(rest) ? p : chargeConjugate(p);
    return true;
  };
  if (name.rfind("anti-", 0) == 0) return conjugated(name.substr(5));
  if (const auto pos = name.find('~'); pos != std::string::npos) {
    return conjugated(name.substr(0, pos) + name.substr(pos + 1));
  }
  if (endsWith(name, "_bar")) return conjugated(name.substr(0, name.size() - 4));
  if (endsWith(name, "bar")) return conjugated(name.substr(0, name.size() - 3));
  return false;
}

}  // namespace

int pdgFromName(const std::string& name) {
  int pdg = 0;
  if (resolve(name, pdg)) return pdg;

  if (auto* p = TDatabasePDG::Instance()->GetParticle(name.c_str())) return p->PdgCode();

  // Integer PDG code, e.g. "-421"
  char* end = nullptr;
  const long code = std::strtol(name.c_str(), &end, 10);
  if (!name.empty() && end != nullptr && *end == '\0' && code != 0) return static_cast<int>(code);

  throw std::invalid_argument("Unknown particle name or PDG code: '" + name + "'");
}

std::string nameFromPdg(int pdg) {
  if (const auto it = pdgMap().find(pdg); it != pdgMap().end()) return it->second->names.front();
  if (const auto it = pdgMap().find(-pdg); it != pdgMap().end()) {
    const auto& n = it->second->names.front();
    if (hasChargeSuffix(n)) return flipChargeSuffix(n);
    // neutral antiparticle, LHCb style: "D0" -> "D~0", "K*(892)0" -> "K*(892)~0"
    return n.substr(0, n.size() - 1) + "~" + n.back();
  }
  if (auto* p = TDatabasePDG::Instance()->GetParticle(pdg)) return p->GetName();
  return std::to_string(pdg);
}

double mass(int pdg) {
  if (const auto it = pdgMap().find(std::abs(pdg)); it != pdgMap().end()) return it->second->mass;
  if (auto* p = TDatabasePDG::Instance()->GetParticle(pdg)) return p->Mass();
  return -1.0;
}

double charge(int pdg) {
  if (const auto it = pdgMap().find(std::abs(pdg)); it != pdgMap().end())
    return pdg > 0 ? it->second->charge : -it->second->charge;
  if (auto* p = TDatabasePDG::Instance()->GetParticle(pdg)) return p->Charge() / 3.0;
  return 0.0;
}

bool isSelfConjugate(int pdg) {
  if (const auto it = pdgMap().find(std::abs(pdg)); it != pdgMap().end()) return it->second->selfConjugate;
  return TDatabasePDG::Instance()->GetParticle(-pdg) == nullptr;
}

bool isResonance(int pdg) {
  switch (std::abs(pdg)) {
    case 111: case 113: case 213: case 221: case 223: case 313: case 323: case 331: case 333:  // light
    case 413: case 423: case 433:                                                              // D*, D*s
    case 443: case 100443:                                                                     // charmonium
    case 513: case 523: case 533:                                                              // B*
      return true;
    default:
      return false;
  }
}

}  // namespace k4DTF::ParticleTable
