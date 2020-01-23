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
#include "warpers_common.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>

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

// **********************

class warper_plasma_uv_t {
 private:
  double t_max = 0.0; // default constructor

 public:
  warper_plasma_uv_t(double t_max_) : t_max{t_max_} {}

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    return ui_from_vi(t_max, li_vec);
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    return vi_from_ui(t_max, ui_vec);
  }

  [[nodiscard]] double jacobian([[maybe_unused]] std::vector<double> const &li_vec) const { return 1.0; }

  [[nodiscard]] double operator()([[maybe_unused]] std::vector<double> const &ui_vec) const { return 1.0; }
};

} // namespace keldy
