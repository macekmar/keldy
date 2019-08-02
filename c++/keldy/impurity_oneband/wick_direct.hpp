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

#include "model.hpp"
#include "keldy/common.hpp"

#include <triqs/utility/first_include.hpp>
#include <triqs/utility/variant.hpp>
#include <triqs/gfs.hpp>
#include <vector>


using namespace triqs::gfs;

namespace keldy::impurity_oneband {

// struct measure_operator_t {
//   std::vector<gf_index_t> creation_operators;
//   std::vector<gf_index_t> annihilation_operators;
// };

class integrand_g_t1t2_direct {
  public:
  g0_keldysh_contour_t g0;
  gf_index_t external_A;
  gf_index_t external_B;
  /// Returns integrand for the specified times
  dcomplex operator()(std::vector<double> const &times) const;
  integrand_g_t1t2_direct(g0_keldysh_contour_t g0_, gf_index_t external_A_, gf_index_t external_B_) :
  g0(std::move(g0_)), external_A(std::move(external_A_)), external_B(std::move(external_B_)) {};
};

// struct integrand_weight_kernel {
//   /// Returns integrand for the specified configuratiton
//   dcomplex operator()(configuration_t const &config) const;
// };

} // namespace anderson_re
