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

#include <complex>

namespace keldy {

// Parameters / definitions used in both real and imag time calculations
using dcomplex = std::complex<double>;
using time_real_t = double;
using time_imag_t = double;
using orbital_t = int;

enum spin_t { up = 0, down = 1 };
enum keldysh_idx_t { forward = 0, backward = 1 };

/// Point of the Keldysh Contour (time, keldysh_idx, timesplit)
class contour_pt_t {
 public:
  time_real_t time{};
  keldysh_idx_t k_idx = forward;
  int timesplit_n = 0; // time-spliting order to dinstiguish vertices at equal times

  bool operator==(const contour_pt_t &other) const {
    return (time == other.time) && (k_idx == other.k_idx) && (timesplit_n == other.timesplit_n);
  }
};

int compare_3way(const contour_pt_t &a, const contour_pt_t &b);

} // namespace keldy
