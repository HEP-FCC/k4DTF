// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/DecayDescriptor.h"

#include <algorithm>
#include <cctype>
#include <functional>

#include "k4DTF/ParticleTable.h"

namespace k4DTF {

namespace {

std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\n\r");
  const auto e = s.find_last_not_of(" \t\n\r");
  return (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
}

[[noreturn]] void syntaxError(const std::string& what, const std::string& in) {
  throw std::invalid_argument("DecayDescriptor: " + what + " in: '" + in + "'");
}

/// Index of the ')' matching the '(' at position open, or npos.
std::size_t matchingParen(const std::string& s, std::size_t open) {
  int depth = 0;
  for (std::size_t i = open; i < s.size(); ++i) {
    if (s[i] == '(') ++depth;
    if (s[i] == ')' && --depth == 0) return i;
  }
  return std::string::npos;
}

/// Position of the "->" that is not inside parentheses, or npos.
std::size_t topLevelArrow(const std::string& s) {
  int depth = 0;
  for (std::size_t i = 0; i + 1 < s.size(); ++i) {
    if (s[i] == '(') ++depth;
    if (s[i] == ')') --depth;
    if (depth == 0 && s[i] == '-' && s[i + 1] == '>') return i;
  }
  return std::string::npos;
}

bool isSpace(char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; }

}  // namespace

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

DecayDescriptor DecayDescriptor::fromString(const std::string& desc) {
  return DecayDescriptor(parseDecay(trim(desc)));
}

DecayDescriptor::DecayDescriptor(Node root) : m_root(std::move(root)) { buildFlat(); }

DecayDescriptor::Node DecayDescriptor::makeParticle(const std::string& name) {
  Node n;
  n.name = name;
  n.pdg = ParticleTable::pdgFromName(name);
  n.mass = ParticleTable::mass(n.pdg);
  return n;
}

DecayDescriptor::Node DecayDescriptor::parseDecay(const std::string& s) {
  const auto arrow = topLevelArrow(s);
  if (arrow == std::string::npos) syntaxError("missing '->'", s);

  const std::string mother = trim(s.substr(0, arrow));
  if (mother.empty()) syntaxError("missing mother particle", s);
  if (std::any_of(mother.begin(), mother.end(), isSpace))
    syntaxError("mother must be a single particle name", s);

  Node root = makeParticle(mother);
  root.daughters = parseDaughters(s.substr(arrow + 2));
  if (root.daughters.empty()) syntaxError("no daughters", s);
  return root;
}

std::vector<DecayDescriptor::Node> DecayDescriptor::parseDaughters(const std::string& s) {
  std::vector<Node> daughters;
  std::size_t i = 0;
  while (i < s.size()) {
    if (isSpace(s[i])) {
      ++i;
      continue;
    }
    if (s[i] == ')') syntaxError("unmatched ')'", s);
    if (s[i] == '(') {
      // nested decay "(X -> a b)"
      const auto close = matchingParen(s, i);
      if (close == std::string::npos) syntaxError("unmatched '('", s);
      daughters.push_back(parseDecay(trim(s.substr(i + 1, close - i - 1))));
      i = close + 1;
      continue;
    }
    // particle name; parentheses inside a name belong to it: "D*(2010)-", "psi(2S)"
    std::string token;
    while (i < s.size() && !isSpace(s[i])) {
      if (s[i] == ')') syntaxError("unmatched ')'", s);
      if (s[i] == '(') {
        const auto close = matchingParen(s, i);
        if (close == std::string::npos) syntaxError("unmatched '('", s);
        token += s.substr(i, close - i + 1);
        i = close + 1;
        continue;
      }
      token += s[i++];
    }
    if (token == "->") syntaxError("unexpected '->' (nested decays need parentheses)", s);
    daughters.push_back(makeParticle(token));
  }
  return daughters;
}

// ---------------------------------------------------------------------------
// Flat view
// ---------------------------------------------------------------------------

void DecayDescriptor::buildFlat() {
  m_flat.clear();
  int nextLeaf = 0;
  std::function<int(const Node&, int, int)> visit = [&](const Node& n, int parent, int depth) -> int {
    const int idx = static_cast<int>(m_flat.size());
    FlatNode f;
    f.pdg = n.pdg;
    f.name = n.name;
    f.mass = n.mass;
    f.parent = parent;
    f.depth = depth;
    f.firstLeaf = nextLeaf;
    m_flat.push_back(f);
    if (n.isLeaf()) {
      m_flat[idx].leafIndex = nextLeaf++;
    } else {
      for (const auto& d : n.daughters) {
        const int child = visit(d, idx, depth + 1);
        m_flat[idx].children.push_back(child);
      }
    }
    m_flat[idx].nLeaves = nextLeaf - m_flat[idx].firstLeaf;
    return idx;
  };
  visit(m_root, -1, 0);
}

std::vector<int> DecayDescriptor::finalStatePdgs() const {
  std::vector<int> pdgs;
  for (const auto& f : m_flat)
    if (f.isLeaf()) pdgs.push_back(f.pdg);
  return pdgs;
}

std::size_t DecayDescriptor::nLeaves() const { return m_flat.empty() ? 0 : m_flat.front().nLeaves; }

std::vector<int> DecayDescriptor::leafNodes() const {
  std::vector<int> out;
  for (std::size_t i = 0; i < m_flat.size(); ++i)
    if (m_flat[i].isLeaf()) out.push_back(static_cast<int>(i));
  return out;
}

// ---------------------------------------------------------------------------
// Charge conjugation
// ---------------------------------------------------------------------------

DecayDescriptor DecayDescriptor::chargeConjugate() const {
  std::function<Node(const Node&)> conj = [&](const Node& n) {
    Node c;
    c.pdg = ParticleTable::chargeConjugate(n.pdg);
    c.name = (c.pdg == n.pdg) ? n.name : ParticleTable::nameFromPdg(c.pdg);
    c.mass = n.mass;
    for (const auto& d : n.daughters) c.daughters.push_back(conj(d));
    return c;
  };
  return DecayDescriptor(conj(m_root));
}

bool DecayDescriptor::isSelfConjugate() const {
  auto a = finalStatePdgs();
  auto b = chargeConjugate().finalStatePdgs();
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  return a == b;
}

// ---------------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------------

std::string DecayDescriptor::toString(int flatIndex) const {
  const auto& f = m_flat.at(flatIndex);
  if (f.isLeaf()) return f.name;
  std::string s = f.name + " ->";
  for (int c : f.children) {
    const auto sub = toString(c);
    s += " " + (m_flat[c].isLeaf() ? sub : "(" + sub + ")");
  }
  return s;
}

std::string DecayDescriptor::toString() const { return m_flat.empty() ? "" : toString(0); }

}  // namespace k4DTF
