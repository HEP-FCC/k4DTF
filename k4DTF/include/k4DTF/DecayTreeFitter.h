// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include "k4DTF/ErrCode.h"
#include "k4DTF/FitResult.h"
#include "k4DTF/core/DecayNodeInput.h"

namespace k4dtf {

class DecayTreeFitter {
public:
  struct Configuration {
    int maxIter{30};
    double convergenceCriterion{1e-6};
    int verbose{0};
  };

  DecayTreeFitter() = default;
  explicit DecayTreeFitter(const Configuration& config);

  FitResult fit(const DecayNodeInput& input,
                std::optional<InteractionPointInput> pv = std::nullopt);

  void setConfiguration(const Configuration& config) { m_config = config; }
  const Configuration& configuration() const { return m_config; }

private:
  Configuration m_config;
};

} // namespace k4dtf
