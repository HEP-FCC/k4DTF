// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace k4DTF {

/// Describes a decay topology for the fitter.
///
/// Grammar:  decay := mother "->" daughter+ ;  daughter := name | "(" decay ")"
///
///   "B0 -> K+ pi-"
///   "B- -> (D0 -> K- pi+) K-"
///   "B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+"
///
/// Names are resolved by k4DTF::ParticleTable (LHCb, EvtGen and ROOT
/// spellings, or integer PDG codes).
///
/// Node numbering: flatten() lists all nodes in depth-first pre-order, so
/// node 0 is the head and every composite node is followed by its subtree.
/// Leaves (final-state particles) are numbered separately, also depth-first;
/// this is the order in which DecayTreeFitter::fit() expects its inputs. The
/// leaves of any composite node form a contiguous range [firstLeaf, firstLeaf + nLeaves).
class DecayDescriptor {
public:
  struct Node {
    int pdg{0};
    std::string name;
    double mass{-1.0};  ///< nominal mass [GeV] from ParticleTable (-1 if unknown)
    std::vector<Node> daughters;

    bool isLeaf() const { return daughters.empty(); }
  };

  /// Flat view of one node, see flatten().
  struct FlatNode {
    int pdg{0};
    std::string name;
    double mass{-1.0};
    int parent{-1};             ///< index of the mother node, -1 for the head
    std::vector<int> children;  ///< indices of the daughter nodes
    int leafIndex{-1};          ///< position among the leaves (-1 for composites)
    int firstLeaf{0};           ///< first leaf below this node
    int nLeaves{0};             ///< number of leaves below this node (1 for a leaf)
    int depth{0};               ///< 0 for the head

    bool isLeaf() const { return children.empty(); }
  };

  /// Build from a decay-descriptor string; throws std::invalid_argument on
  /// syntax errors or unknown particle names.
  static DecayDescriptor fromString(const std::string& desc);

  /// Build directly from a pre-constructed node tree.
  explicit DecayDescriptor(Node root);

  const Node& root() const { return m_root; }

  /// Flat list of final-state (leaf) PDG codes in tree order.
  std::vector<int> finalStatePdgs() const;

  /// All nodes in depth-first pre-order (node 0 = head), see class comment.
  const std::vector<FlatNode>& flatten() const { return m_flat; }

  std::size_t nNodes() const { return m_flat.size(); }
  std::size_t nLeaves() const;

  /// Indices (into flatten()) of the leaf nodes, in leaf order.
  std::vector<int> leafNodes() const;

  /// The same topology with every particle replaced by its antiparticle.
  DecayDescriptor chargeConjugate() const;

  /// True if the charge-conjugate descriptor has the same final state,
  /// e.g. "B0 -> (J/psi -> mu+ mu-) (KS0 -> pi+ pi-)".
  bool isSelfConjugate() const;

  /// Canonical string with nested decays in parentheses; parses back to the
  /// same tree.
  std::string toString() const;

  /// Canonical string of the subtree below one flat node.
  std::string toString(int flatIndex) const;

private:
  Node m_root;
  std::vector<FlatNode> m_flat;

  void buildFlat();
  static Node parseDecay(const std::string& s);
  static std::vector<Node> parseDaughters(const std::string& s);
  static Node makeParticle(const std::string& name);
};

}  // namespace k4DTF
