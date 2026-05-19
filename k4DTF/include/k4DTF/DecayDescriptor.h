// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace k4DTF {

/// Describes a decay topology for the fitter.
///
/// Simple format: "B0 -> K+ pi-" or nested "B0 -> (D0 -> K- pi+) pi+"
/// PDG codes or particle names may be used; the fitter queries TDatabasePDG.
class DecayDescriptor {
public:
  struct Node {
    int pdg{0};
    std::string name;
    std::vector<Node> daughters;

    bool isLeaf() const { return daughters.empty(); }
  };

  /// Build from a decay-descriptor string, e.g. "B0 -> K+ pi-".
  static DecayDescriptor fromString(const std::string& desc);

  /// Build directly from a pre-constructed node tree.
  explicit DecayDescriptor(Node root);

  const Node& root() const { return m_root; }

  /// Flat list of final-state (leaf) PDG codes in tree order.
  std::vector<int> finalStatePdgs() const;

  std::string toString() const;

private:
  Node m_root;

  static Node parseNode(const std::string& token);
  static std::vector<Node> parseDaughters(const std::string& s);
};

}  // namespace k4DTF