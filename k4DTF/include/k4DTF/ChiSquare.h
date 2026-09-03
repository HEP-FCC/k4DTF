// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace k4dtf {

struct ChiSquare {
  double value{0};
  int    ndf{0};

  ChiSquare() = default;
  ChiSquare(double chi, int n) : value(chi), ndf(n) {}

  ChiSquare& operator+=(const ChiSquare& rhs) {
    value += rhs.value;
    ndf += rhs.ndf;
    return *this;
  }
};

} // namespace k4dtf
