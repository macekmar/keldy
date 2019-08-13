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

#include "keldy/common.hpp"

#include <triqs/gfs.hpp>
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/variant.hpp>
#include <vector>

using namespace triqs::gfs;

namespace keldy::impurity_oneband {

class CPP2PY_IGNORE model_param_t {
 public:
  double beta = 1.0;
  double bias_V = 0.0;
  double eps_d = 0.0;
  double Gamma = 1.0;
  double time_max = +100.0; // (time_limit_min = - time_limit_max)
  int nr_time_points_gf = 1000;
  double alpha = 0.0;
  std::string bath_type = "flatband";
};

// fake function to get cpp2py to create adaptor for model_param_t
CPP2PY_ARG_AS_DICT inline void fake(model_param_t const &temp){};

/// Point of the Contour Keldysh Green Function (time, spin, keldysh_idx)
class gf_index_t {
 public:
  // time / contour related:
  time_real_t time;
  keldysh_idx_t k_idx;
  int timesplit_n = 0; // time-spliting order to dinstiguish vertices at equal times

  spin_t spin;

  /// Constructor: (time, spin, keldysh_idx)
  gf_index_t() = default;
  gf_index_t(time_real_t time_, int spin_, int k_idx_)
     : time(time_), k_idx(keldysh_idx_t(k_idx_)), spin(spin_t(spin_)) {}
  gf_index_t(time_real_t time_, int spin_, int k_idx_, int timesplit_n_)
     : time(time_), k_idx(keldysh_idx_t(k_idx_)), timesplit_n(timesplit_n_), spin(spin_t(spin_)) {}
};

/// Defines model throuh non-interacting Green function g_lesser / g_greater
class g0_model {
 public:
  /// Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, scalar_valued> g0_lesser;

  /// Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, scalar_valued> g0_greater;

  explicit g0_model(model_param_t const &parameters);

  void make_semicircular_model();
  void make_flat_band();

  model_param_t param_; // g0_keldysh_contour_t will need access to alpha
};

/// Adapt g0_lesser and g0_greater into Green function on Keldysh contour
struct g0_keldysh_contour_t {
  g0_model model;
  /// Evalutate G, passing two Keldysh contour points
  dcomplex operator()(gf_index_t const &a, gf_index_t const &b, bool internal_point = true) const;
  g0_keldysh_contour_t(g0_model model_) : model(std::move(model_)){};
};

// class g0_keldysh_LO_t {
//   // g_ai_t g0_retarded()
//   // g_ai_t g0_advanced()
//   // g_ai_t g0_keldysh()
// };

} // namespace keldy::impurity_oneband
