// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include "k4DTF/DTFParticle.h"
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/core/DecayNodeInput.h"

namespace k4DTF {

/// Builds the input of k4dtf::DecayTreeFitter::fit() from a decay descriptor
/// and the measured final-state particles.
///
/// - leaves: one DTFParticle per leaf of the descriptor, in depth-first leaf
///   order (DecayDescriptor::leafNodes()). Charged leaves with a track become
///   helix measurements; charged leaves without one use their 3-momentum and
///   the (px, py, pz) block of DTFParticle::cov. Neutral leaves are not
///   supported yet (std::invalid_argument).
/// - Masses: every node gets the nominal mass of its PDG code
///   (ParticleTable); for leaves that is the mass hypothesis of the track.
/// - Composites decaying strongly or electromagnetically
///   (ParticleTable::isResonance) are resonances, the others get a vertex.
/// - massConstraints: PDG codes of composites to mass-constrain to their
///   nominal mass; a code also applies to its charge conjugate.
///
/// FitResult::nodes is then in DecayDescriptor::flatten() order.
k4dtf::DecayNodeInput makeDecayNodeInput(const DecayDescriptor& desc, const std::vector<DTFParticle>& leaves,
                                         const std::vector<int>& massConstraints = {});

/// Primary vertex for the fitter.
k4dtf::InteractionPointInput makeInteractionPoint(const DTFVertex& v);

}  // namespace k4DTF
