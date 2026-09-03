// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#include "k4DTF/FitResult.h"

#include <cmath>

namespace k4dtf {

double FitResult::prob() const {
  if (ndf <= 0) return 0.0;
  return 0.0;
}

} // namespace k4dtf
