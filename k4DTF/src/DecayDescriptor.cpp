// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DecayDescriptor.h"

#include <algorithm>
#include <functional>
#include <sstream>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

namespace k4DTF {

namespace {

int pdgFromName(const std::string& name) {
  auto* db = TDatabasePDG::Instance();
  auto* p = db->GetParticle(name.c_str());
  if (!p) {
    // Try parsing as integer PDG code
    try {
      return std::stoi(name);
    } catch (...) {
      throw std::invalid_argument("Unknown particle name or PDG code: " + name);
    }
  }
  return static_cast<int>(p->PdgCode());
}

std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\n\r");
  const auto e = s.find_last_not_of(" \t\n\r");
  return (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
}

}  // namespace

// ---------------------------------------------------------------------------

DecayDescriptor DecayDescriptor::fromString(const std::string& desc) {
  const auto arrowPos = desc.find("->");
  if (arrowPos == std::string::npos)
    throw std::invalid_argument("DecayDescriptor: missing '->' in: " + desc);

  const std::string motherStr = trim(desc.substr(0, arrowPos));
  const std::string daughtersStr = trim(desc.substr(arrowPos + 2));

  Node root;
  root.name = motherStr;
  root.pdg = pdgFromName(motherStr);
  root.daughters = parseDaughters(daughtersStr);

  return DecayDescriptor(std::move(root));
}

DecayDescriptor::DecayDescriptor(Node root) : m_root(std::move(root)) {}

// ---------------------------------------------------------------------------

std::vector<DecayDescriptor::Node> DecayDescriptor::parseDaughters(
    const std::string& s) {
  std::vector<Node> daughters;
  // TODO: implement nested "(X -> a b)" syntax
  // For now, split on whitespace at the top level
  std::istringstream ss(s);
  std::string token;
  while (ss >> token) {
    Node n;
    n.name = token;
    n.pdg = pdgFromName(token);
    daughters.push_back(std::move(n));
  }
  return daughters;
}

// ---------------------------------------------------------------------------

std::vector<int> DecayDescriptor::finalStatePdgs() const {
  std::vector<int> pdgs;
  std::function<void(const Node&)> collect = [&](const Node& n) {
    if (n.isLeaf()) {
      pdgs.push_back(n.pdg);
    } else {
      for (const auto& d : n.daughters) collect(d);
    }
  };
  collect(m_root);
  return pdgs;
}

std::string DecayDescriptor::toString() const {
  std::function<std::string(const Node&)> fmt = [&](const Node& n) -> std::string {
    if (n.isLeaf()) return n.name;
    std::string s = n.name + " -> ";
    for (const auto& d : n.daughters) s += fmt(d) + " ";
    return s;
  };
  return fmt(m_root);
}

}  // namespace k4DTF
