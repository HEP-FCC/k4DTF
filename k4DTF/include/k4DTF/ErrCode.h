// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <iostream>

namespace k4dtf {

enum class ErrCode : unsigned int {
  success             = 0,
  pocafailure         = 1,
  baddistance         = 2,
  inversionerror      = 4,
  badsetup            = 8,
  divergingconstraint = 16,
  slowdivergingfit    = 32,
  fastdivergingfit    = 64,
  filtererror         = 128
};

inline ErrCode& operator|=(ErrCode& a, ErrCode b) {
  a = static_cast<ErrCode>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
  return a;
}

inline ErrCode operator|(ErrCode a, ErrCode b) { return a |= b; }

inline bool failure(ErrCode e) { return static_cast<unsigned int>(e) != 0; }

inline std::ostream& operator<<(std::ostream& os, ErrCode ec) {
  auto flag = static_cast<unsigned int>(ec);
  os << "flag=" << flag << " ";
  if (flag) {
    if (flag & static_cast<unsigned int>(ErrCode::pocafailure)) os << "pocafailure ";
    if (flag & static_cast<unsigned int>(ErrCode::baddistance)) os << "baddistance ";
    if (flag & static_cast<unsigned int>(ErrCode::inversionerror)) os << "inversionerror ";
    if (flag & static_cast<unsigned int>(ErrCode::badsetup)) os << "badsetup ";
    if (flag & static_cast<unsigned int>(ErrCode::divergingconstraint)) os << "divergingconstraint ";
    if (flag & static_cast<unsigned int>(ErrCode::slowdivergingfit)) os << "slowdivergingfit ";
    if (flag & static_cast<unsigned int>(ErrCode::fastdivergingfit)) os << "fastdivergingfit ";
    if (flag & static_cast<unsigned int>(ErrCode::filtererror)) os << "filtererror ";
  } else {
    os << "success ";
  }
  return os;
}

} // namespace k4dtf
