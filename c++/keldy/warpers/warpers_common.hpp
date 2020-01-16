//******************************************************************************
//
// keldy
//
// Copyright (C) 2019, The Simons Foundation
// authors: Philipp Dumitrescu
//
// keldy is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// keldy is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// keldy. If not, see <http://www.gnu.org/licenses/>.
//
//******************************************************************************

#pragma once

#include "../common.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <triqs/gfs.hpp>

namespace keldy {

// Coordinate Transforms (ui <-> vi):
inline std::vector<double> vi_from_ui(double t_max, std::vector<double> const &u_times) {
  std::vector<double> v_times = u_times;
  std::sort(v_times.begin(), v_times.end(), std::greater<>());
  int n = v_times.size();
  for (int i = n - 1; i > 0; i--) {
    v_times[i] = (v_times[i - 1] - v_times[i]);
  }
  v_times[0] = t_max - v_times[0];
  return v_times;
}

inline std::vector<double> ui_from_vi(double t_max, std::vector<double> const &v_times) {
  std::vector<double> u_times = v_times;
  u_times[0] = t_max - u_times[0];
  std::partial_sum(u_times.begin(), u_times.end(), u_times.begin(), std::minus<>());
  return u_times;
}

// Some Simple Functions:

struct CPP2PY_IGNORE identity_function {
  double operator()([[maybe_unused]] double t) { return 1.0; }
};

struct CPP2PY_IGNORE linear_function {
  double operator()(double t) { return t; }
};

} // namespace keldy
