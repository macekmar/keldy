/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019 The Simons Foundation
 *   authors: Philipp Dumitrescu, Corentin Bertrand
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

#pragma once

#include "../common.hpp"
#include "../binner.hpp"
#include "model.hpp"

#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <triqs/arrays/math_functions.hpp>
#include <triqs/arrays/mpi.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/first_include.hpp>
#include <algorithm>
#include <vector>
#include <tuple>

namespace keldy::impurity_oneband {

using namespace triqs::arrays;
using namespace triqs::gfs;

class integrand_g_kernel_single_omega {
  g0_keldysh_contour_t g0;
  gf_index_t g_idx_X; // Fixed Point in Kernal
  double const omega;

  // bool expand_col = true; // expand_row = false
  // double condition_numebr_tol;

 public:
  /// Returns integrand for the specified times
  using result_t = dcomplex;
  [[nodiscard]] std::pair<result_t, int> operator()(std::vector<double> const &times,
                                                    bool const keep_u_hypercube = true) const;

  integrand_g_kernel_single_omega(g0_keldysh_contour_t g0_, gf_index_t g_idx_X_, double omega_)
     : g0(std::move(g0_)), g_idx_X(g_idx_X_), omega(omega_){};
};

} // namespace keldy::impurity_oneband
