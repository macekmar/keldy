/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019 The Simons Foundation
 *   authors: Philipp Dumitrescu
 *
 * keldy is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * keldy is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * keldy. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "common.hpp"
#include <cmath>

namespace keldy {

/// Time ordering along the Keldysh Basis:
///
/// 3-way compare (spaceship): a <=> b:
/// (a <=> b) < 0  if a < b
/// (a <=> b) > 0  if a > b
/// (a <=> b) == 0 if a == b
///
/// Mapping takes care of correct keldysh ordering [0 = forward contour; 1 = backward contour]
///  a    b    (a.time > b.time)   L/G      a <=> b
///  0    0           1             G       a  >  b
///  0    0           0             L       a  <  b
///  1    1           1             L       a  >  b
///  1    1           0             G       a  >  b
///
///  0    1           *             L       a  <  b
///  1    0           *             G       a  >  b
///
/// At equal times and Keldysh index we need to point-split times according to external integer (time-ordering).
/// For self contractions (no external ordering), use $g^< \sim c^\dag c$ since this is normal ordering defined by V
int compare_3way(const contour_pt_t &a, const contour_pt_t &b) {
  int k_idx_3way = a.k_idx - b.k_idx;
  if (k_idx_3way != 0) {
    return k_idx_3way;
  }
  auto time_3way = a.time - b.time;
  if (time_3way != 0.0) {
    return int(1 - 2 * int(a.k_idx)) * (2 * int(!std::signbit(time_3way)) - 1);
  }
  return (1 - 2 * int(a.k_idx)) * (a.timesplit_n - b.timesplit_n);
}

} // namespace keldy
