// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

namespace k4DTF {

/// Particle names, PDG codes and nominal masses used by the decay descriptor.
///
/// A built-in table covers the hadrons and leptons common in flavour decay
/// chains, accepting the usual spellings:
///   LHCb       "D~0", "B~0", "D*(2010)-", "KS0", "J/psi(1S)", "Lambda_c+"
///   EvtGen     "anti-D0", "anti-B0", "D*-", "K_S0", "J/psi"
///   other      "D0bar", "D0_bar", "Lc+", "Lb0"
/// Unknown names fall back to TDatabasePDG, then to an integer PDG code.
/// Masses [GeV] in the table are those of EvtGen's evt.pdl as shipped in
/// FCC-config, i.e. the masses the FCC samples are generated with, so mass
/// constraints in validation fits are unbiased. Masses of TDatabasePDG
/// fallbacks come from ROOT's table.
namespace ParticleTable {

/// PDG code for a particle name; throws std::invalid_argument if unknown.
int pdgFromName(const std::string& name);

/// Canonical (LHCb-style) name for a PDG code, or the code as a string.
std::string nameFromPdg(int pdg);

/// Nominal mass [GeV]; negative if the particle is unknown.
double mass(int pdg);

/// Charge [e]; 0 for unknown particles.
double charge(int pdg);

/// True for particles that are their own antiparticle (gamma, pi0, KS0, J/psi, ...).
bool isSelfConjugate(int pdg);

/// True for particles that decay strongly or electromagnetically (D*, K*, rho,
/// phi, J/psi, pi0, eta, ...): in a decay-tree fit they have no flight length.
/// Weakly decaying hadrons and leptons (B, D, Lambda_c, KS0, tau, ...) are not.
bool isResonance(int pdg);

/// PDG code of the charge-conjugate particle.
inline int chargeConjugate(int pdg) { return isSelfConjugate(pdg) ? pdg : -pdg; }

}  // namespace ParticleTable
}  // namespace k4DTF
