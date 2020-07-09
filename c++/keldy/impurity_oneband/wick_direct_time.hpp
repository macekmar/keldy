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
#include "../binner.hpp"
#include "model.hpp"

#include <triqs/utility/first_include.hpp>
#include <triqs/utility/variant.hpp>
#include <vector>

namespace keldy::impurity_oneband {

class integrand_g_direct_time {
  g0_keldysh_contour_t g0;
  gf_index_t external_A;
  gf_index_t external_B;
  double cutoff;

 public:
  // Specify return type of function to easily check type compability
  using result_t = binner::sparse_binner_t<1>;

  /// Returns integrand for the specified times
  [[nodiscard]] std::pair<result_t, int> operator()(std::vector<double> const &times,
                                                    bool const keep_u_hypercube = true) const;

  integrand_g_direct_time(g0_keldysh_contour_t g0_, gf_index_t external_A_, gf_index_t external_B_, double cutoff_ = 0.)
     : g0(std::move(g0_)), external_A(external_A_), external_B(external_B_), cutoff(cutoff_) {

    if ((external_A.contour.time < 0.) || (external_B.contour.time < 0.)) {
      TRIQS_RUNTIME_ERROR << "An external point has negative time.";
    }
    double const t_max = g0.get_time_max();
    if ((external_A.contour.time > t_max) || (external_B.contour.time > t_max)) {
      TRIQS_RUNTIME_ERROR << "An external point is out of g0 time window.";
    }
    if ((external_A.orbital != 0 || external_B.orbital != 0) && (!g0.model.make_dot_lead)) {
      TRIQS_RUNTIME_ERROR << "Need g0_model with make_dot_lead == true to calculate off_diagonal gf.";
    }
  };
};

} // namespace keldy::impurity_oneband
